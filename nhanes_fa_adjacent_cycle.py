#!/usr/bin/env python3
"""Evaluate fixed NHANES block choices with the frozen within-cycle RBF recipe."""
from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from hashlib import file_digest, sha256
import importlib.util
import json
from pathlib import Path
import platform
import re
import sys

import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata, t


BENCHMARK = 0.11828095881172637
BOOT_SEED = 2026111002
BOOT_DRAWS = 5000
CHOICES = ("LBXTSA", "LBXEPA", "LBXDTA", "LBXEPA", "LBXDTA")
PROTOCOL = "NHANES_FA_ADJACENT_CYCLE_PROTOCOL_v0.1.md"
PINNED = {
    "nhanes_fa_run.py": "f85d5235cd2fd609145be9ba2b444cd2c4b394a8b3ed6c761c942e765a730297",
    "nhanes_fa/frozen/p6_n2.py": "681ce852609e6ab4fcaf0f0ff603477c64e2351e10bd119811d45f008d4bbe8e",
    "nhanes_fa/frozen/p5_gate3a.py": "4600e2f3e3bb652affadcc6f6935380d8b9a77a4237e7cb306f056f849478130",
    "nhanes_fa_blocks_v0.4.json": "e3e1d8ac636afd9690e4b5b00b884e21ccf79dad77993e4a75d8776ef3a14338",
    "nhanes_fa_stab_out/dev_fit.json": "e9d894aa2e0b549ce4c02712cbbbcb2fd679ffddbac85ce9539420125dc46669",
    "nhanes_fa_out/FAS_G_development_v04_result.json": "f7f9f5cad1c25edcd9e395ae4198b5997880d80c58552fee4aa77b91121c2101",
}
DEV_INPUTS = {
    "FAS_G.xpt": "5a8a907c9797b00d4f9d4dfbf0613e68d14d833d332f2aa3d9dfbedb3525063b",
    "DEMO_G.xpt": "eaf0525d1952626885af3e935415a1f66ad62c18698080e7354789c125af252d",
}


class Unassessable(ValueError):
    """A mandatory protocol requirement failed."""


def require(condition, message):
    if not condition:
        raise Unassessable(message)


def digest(path):
    with path.open("rb") as stream:
        return file_digest(stream, "sha256").hexdigest()


def read(path):
    return json.loads(path.read_text())


def write_new(path, value):
    payload = json.dumps(value, indent=2, allow_nan=False) + "\n"
    with path.open("x") as stream:
        stream.write(payload)


def utc():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def timestamp(text):
    value = datetime.fromisoformat(text)
    require(value.utcoffset() is not None, "opening times must include a time zone")
    return value


def verify_opening(path, root):
    """Check the separately authorised, independently timestamped opening record."""
    a = read(path)
    require(a["cycle"] == "2013–14", "opening authorisation is for the wrong cycle")
    require(a["protocol_sha256"] == digest(root / PROTOCOL), "protocol hash mismatch")
    require(a["runner_sha256"] == digest(Path(__file__)), "runner hash mismatch")
    require(re.fullmatch(r"[0-9a-f]{40}", a["public_commit"]), "missing public commit")
    require(re.fullmatch(r"[0-9a-f]{64}", a["timestamp_receipt_sha256"]),
            "missing verified timestamp receipt hash")
    stamped = timestamp(a["timestamp_utc"])
    authorised = timestamp(a["authorised_utc"])
    require(stamped <= authorised <= datetime.now(timezone.utc),
            "timestamp must precede authorisation")
    records = a["downloads"]
    require(set(records) == {"FAS_H.xpt", "DEMO_H.xpt"}, "wrong download inventory")
    for record in records.values():
        require(stamped < timestamp(record["download_started_utc"]) <= datetime.now(timezone.utc),
                "download must start after the independent timestamp")
        require(re.fullmatch(r"[0-9a-f]{64}", record["sha256"]), "missing download hash")
    return a


def load_design(wrapper, fas_path, demo_path, columns, cycle, dev_record):
    fas = pd.read_sas(fas_path, format="xport")
    demo = pd.read_sas(demo_path, format="xport")
    cycle_code = 7 if cycle == "G" else 8
    require("SDDSRVYR" in demo and (demo.SDDSRVYR == cycle_code).all(),
            "demographic file cycle code does not match requested cycle")
    for frame in (fas, demo):
        ids = frame["SEQN"].to_numpy(float)
        require(np.isfinite(ids).all() and (ids > 0).all() and (ids == np.floor(ids)).all(),
                "SEQN must be finite positive integers")
    d = wrapper.load_cycle(fas_path, demo_path, {"O": columns, "P": [], "T": []},
                           _frames=(fas, demo))
    age = d.RIDAGEYR.to_numpy(float)
    weight = d.WTFAS2YR.to_numpy(float)
    require(np.isfinite(age).all() and (age >= 0).all(), "invalid age in joined file")
    require(np.isfinite(weight).all() and (weight >= 0).all(), "invalid subsample weight")
    eligible = d.eligible.to_numpy()
    require(eligible.sum() >= 10, "fewer than ten eligible rows for five outer and inner folds")
    require(np.isfinite(d.loc[eligible, columns].to_numpy(float)).all(),
            "eligible analytes must be finite")
    ids = d.loc[eligible, "SEQN"].astype("int64").astype(str).tolist()
    if cycle != "G":
        dev_ids = set(dev_record["seqn_sha256_list"])
        require(not any(sha256(str(int(i)).encode()).hexdigest() in dev_ids for i in fas.SEQN),
                "FAS identifiers overlap the eligible development participants")
    return d, weight, eligible, ids


class Survey:
    def __init__(self, d, weight, eligible, wrapper):
        # Keep all positive-weight PSUs, including those with zero domain contribution.
        strata = d.SDMVSTRA.to_numpy(float)
        psu = d.SDMVPSU.to_numpy(float)
        for codes in (strata, psu):
            require(np.isfinite(codes).all() and (codes > 0).all() and
                    (codes == np.floor(codes)).all(), "invalid masked design codes")
        self.strata = strata[weight > 0]
        self.psu = psu[weight > 0]
        self.domain = eligible[weight > 0]
        self.weight = weight[eligible]
        self.wrapper = wrapper
        self.groups = []
        counts = []
        for h in np.unique(self.strata):
            mh = self.strata == h
            full = np.unique(self.psu[mh])
            inside = np.unique(self.psu[mh & self.domain])
            counts.append({"stratum": int(h), "full_psus": len(full), "domain_psus": len(inside)})
            if len(inside):
                require(len(full) >= 2, "one full-design PSU in a stratum with eligible rows")
                self.groups.append([np.flatnonzero((mh & (self.psu == p))[self.domain])
                                    for p in full])
        n_strata = sum(row["domain_psus"] > 0 for row in counts)
        n_psus = sum(row["domain_psus"] for row in counts)
        self.df = n_psus - n_strata
        require(self.df > 0, "eligible domain has no positive design degrees of freedom")
        self.counts = {"eligible": int(eligible.sum()), "full_design_rows": len(self.strata),
                       "domain_strata": n_strata, "domain_psus": n_psus, "df": self.df,
                       "psu_support_by_stratum": counts}
        self.critical = float(t.ppf(.975, self.df))

    def stat(self, num, den):
        require(np.isfinite(num).all() and np.isfinite(den).all() and
                np.isfinite(self.weight @ den) and self.weight @ den > 0,
                "nonfinite contrast or nonpositive weighted base loss")
        y, b, w = (np.zeros(len(self.domain)) for _ in range(3))
        y[self.domain], b[self.domain], w[self.domain] = num, den, self.weight
        value, se, _, _, _ = self.wrapper.ratio_linearised(y, b, w, self.strata, self.psu)
        half = self.critical * se
        require(np.isfinite([value, se, half]).all(), "nonfinite survey result")
        return {"estimate": value, "se": se, "ci95": [value-half, value+half]}

    def bootstrap_weights(self):
        rng = np.random.default_rng(BOOT_SEED)
        draws = []
        for group in self.groups:
            n = len(group)
            counts = rng.multinomial(n-1, np.full(n, 1/n), size=BOOT_DRAWS)
            draws.append(counts * n/(n-1))
        return np.concatenate(draws, axis=1)

    def psu_sums(self, mat):
        weighted = mat * self.weight[:, None]
        return np.stack([weighted[rows].sum(axis=0)
                         for group in self.groups for rows in group])


def block_losses(n2, gate, d, block, scale, fold, index):
    e = d[d.eligible]
    O, X, Y = (e[block[key]].to_numpy(float) for key in ("O", "P", "T"))
    Y = Y / scale
    sq, bd = (np.full((len(Y), X.shape[1]+1), np.nan) for _ in range(2))
    picks = []
    for k in range(n2.K_OUT):
        te = fold == k
        tr = ~te
        picks.append(gate.g1_pick(O[tr], Y[tr], X[tr]))
        _, best = n2.inner_scores("rbf", O[tr], X[tr], Y[tr], n2.SEED+k)
        for j in range(X.shape[1]+1):
            args = (X[tr, j-1], X[te, j-1]) if j else ()
            pred = n2.predict_grid("rbf", O[tr], O[te], Y[tr], *args)[best[j]]
            errors = (pred-Y[te])**2
            sq[te, j] = errors.sum(axis=1)
            bd[te, j] = np.minimum(errors, n2.CAP).sum(axis=1)
        print(f"block {index+1}/5, fold {k+1}/5 complete", file=sys.stderr, flush=True)
    require(np.isfinite(sq).all() and np.isfinite(bd).all(), "nonfinite or missing held-out losses")
    return sq, bd, picks


def spearman(scores, realised):
    a = rankdata(scores)
    b = rankdata(realised, axis=-1)
    a = a - a.mean()
    b = b - b.mean(axis=-1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        return (b @ a) / np.sqrt((a*a).sum() * (b*b).sum(axis=-1))


def finite_or_null(value):
    return float(value) if np.isfinite(value) else None


def summaries(fitted, survey, blocks, choices, scores):
    base_sq, chosen_sq, base_bd, chosen_bd = (np.zeros(len(survey.weight)) for _ in range(4))
    gains, ranks, rank_draws, swaps, disagreements = [], [], [], [], []
    boot = survey.bootstrap_weights()
    for i, ((sq, bd, picks), block, j, s) in enumerate(zip(fitted, blocks, choices, scores), 1):
        base_sq += sq[:, 0]
        chosen_sq += sq[:, j+1]
        base_bd += bd[:, 0]
        chosen_bd += bd[:, j+1]
        base = sq[:, 0]
        aug = sq[:, 1:]
        gains.append({"block": i, "choice": block["P"][j],
                      **survey.stat(base-aug[:, j], base)})
        num, den = survey.psu_sums(base[:, None]-aug), survey.psu_sums(base[:, None])
        realised = num.sum(axis=0) / den.sum()
        ranks.append({"block": i, "spearman": finite_or_null(spearman(s, realised)),
                      "candidate_gains": [{"candidate": p, "gain": float(r)}
                                          for p, r in zip(block["P"], realised)]})
        with np.errstate(divide="ignore", invalid="ignore"):
            replicate = (boot @ num)/(boot @ den)
        valid = np.isfinite(replicate).all(axis=1)
        draws = np.full(BOOT_DRAWS, np.nan)
        draws[valid] = spearman(s, replicate[valid])
        rank_draws.append(draws)
        order = np.argsort(-np.asarray(s), kind="stable")
        pairs = {"adjacent": zip(order[:-1], order[1:]),
                 "winner_vs_others": ((j, k) for k in order if k != j)}
        swaps.append({"block": i, **{
            kind: [{"higher_rank": block["P"][a], "lower_rank": block["P"][b],
                    **survey.stat(aug[:, b]-aug[:, a], base)} for a, b in pp]
            for kind, pp in pairs.items()}})
        disagreements.append({"block": i, "count": sum(
            picks[a] != picks[b] for a in range(len(picks)) for b in range(a+1, len(picks)))})
    primary = survey.stat(base_sq-chosen_sq, base_sq)
    require(primary["se"] > 0, "primary design standard error is zero")
    lo, hi = primary["ci95"]
    primary["benchmark"] = BENCHMARK
    primary["class"] = "exceeded" if lo > BENCHMARK else "missed" if hi < BENCHMARK else "unresolved"
    draw_mean = np.mean(rank_draws, axis=0)
    complete = int(np.isfinite(draw_mean).sum())
    defined = all(row["spearman"] is not None for row in ranks)
    return {
        "E1": primary,
        "per_block_squared_gain": gains,
        "bounded_loss_sensitivity": survey.stat(base_bd-chosen_bd, base_bd),
        "E2": {"blocks": ranks,
               "mean_within_block_spearman": float(np.mean([row["spearman"] for row in ranks]))
               if defined else None,
               "conditional_percentile_range": np.quantile(draw_mean, [.025, .975]).tolist()
               if defined and complete == BOOT_DRAWS else None,
               "defined_replicates": complete, "total_replicates": BOOT_DRAWS},
        "E3": {"blocks": swaps},
        "observed_refit_disagreement": {"blocks": disagreements,
                                        "total": sum(row["count"] for row in disagreements),
                                        "pairs": 50},
    }


def evaluate(a):
    root = a.project_root
    for name, expected in PINNED.items():
        require(digest(root / name) == expected, f"frozen input changed: {name}")
    spec = importlib.util.spec_from_file_location("frozen_wrapper", root / "nhanes_fa_run.py")
    w = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(w)
    n2 = w.frozen()
    gate = n2.gate3a()
    fit = read(root / "nhanes_fa_stab_out/dev_fit.json")
    dev_record = read(root / "nhanes_fa_out/FAS_G_development_v04_result.json")
    blocks = read(root / "nhanes_fa_blocks_v0.4.json")["blocks"]
    choices = fit["full_fit_picks"]
    require(tuple(b["P"][j] for b, j in zip(blocks, choices)) == CHOICES,
            "development choices do not match the locked list")
    require(all(np.isfinite(s).all() and len(s) == 9 and np.argmax(s) == j
                for s, j in zip(fit["full_fit_scores"], choices)), "invalid fixed score ranks")
    scales = [np.maximum(x["target_sd"], x["floor"]) for x in dev_record["standardisation"]]
    require(len(scales) == 5 and all(len(s) == len(b["T"]) and np.isfinite(s).all() and (s > 0).all()
                                    for b, s in zip(blocks, scales)), "invalid fixed target scales")
    suffix = a.cycle if a.cycle != "synthetic" else "SYNTH"
    names = (f"FAS_{suffix}.xpt", f"DEMO_{suffix}.xpt")
    paths = [a.data_dir / n for n in names]
    if a.cycle == "H":
        require(a.authorisation is not None, "separate authorisation required before H file access")
        opening = verify_opening(a.authorisation, root)
        # One attempt per protocol checkout. Record it before inspecting or hashing either file.
        write_new(root / "NHANES_FA_ADJACENT_CYCLE_H_ATTEMPT.json",
                  {"started_utc": utc(), "authorisation_sha256": digest(a.authorisation),
                   "protocol_sha256": opening["protocol_sha256"], "runner_sha256": opening["runner_sha256"],
                   "public_commit": opening["public_commit"],
                   "timestamp_receipt_sha256": opening["timestamp_receipt_sha256"]})
        expected = {name: item["sha256"] for name, item in opening["downloads"].items()}
    else:
        expected = DEV_INPUTS if a.cycle == "G" else {}
    hashes = {name: digest(path) for name, path in zip(names, paths)}
    require(not expected or hashes == expected, "cycle input hashes do not match the recorded downloads")
    d, weight, eligible, ids = load_design(w, *paths, fit["allvars"], a.cycle, dev_record)
    survey = Survey(d, weight, eligible, w)
    fold = gate.folds_for(ids, w.SEED_DEV if a.cycle == "G" else w.SEED_EXT)
    require(fold.shape == (len(ids),) and set(fold) == set(range(n2.K_OUT)),
            "invalid outer fold partition")
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        futures = [pool.submit(block_losses, n2, gate, d, b, s, fold, i)
                   for i, (b, s) in enumerate(zip(blocks, scales))]
        fitted = [future.result() for future in futures]
    if a.cycle == "G":
        picks = [p for _, _, block_picks in fitted for p in block_picks]
        require(picks == [row["G1"] for row in dev_record["inference"]["picks"]],
                "development refit picks differ from frozen development record")
    result = summaries(fitted, survey, blocks, choices, fit["full_fit_scores"])
    if a.cycle == "G":
        for loss, key in (("sq", "E1"), ("bd", "bounded_loss_sensitivity")):
            target = dev_record["inference"][loss]
            require(np.allclose([result[key]["estimate"], result[key]["se"]],
                                [target["delta_g1_base"], target["se_g1_base"]],
                                rtol=1e-7, atol=1e-10), f"development {loss} did not reproduce")
        require(survey.counts["eligible"] == 1725 and survey.df == 17,
                "development domain counts did not reproduce")
        result["E1"].pop("class")
    elif a.cycle == "synthetic":
        result["E1"].pop("class")
    return {"status": "assessable", "cycle": a.cycle, "input_sha256": hashes,
            "runner_sha256": digest(Path(__file__)), "pinned_sha256": PINNED,
            "software": {"python": platform.python_version(), "numpy": np.__version__,
                         "pandas": pd.__version__, "scipy": scipy.__version__},
            "design": survey.counts, **result}


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--cycle", choices=("G", "H", "synthetic"), required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--workers", type=int, choices=range(1, 6), default=5)
    parser.add_argument("--authorisation", type=Path)
    a = parser.parse_args()
    if a.out.exists() or a.out.is_symlink():
        parser.error("output already exists")
    try:
        result = evaluate(a)
    except (Exception, SystemExit) as exc:
        # Diagnostics deliberately exclude raw exception text, which may contain row values.
        reason = str(exc) if isinstance(exc, Unassessable) else type(exc).__name__
        write_new(a.out, {"status": "unassessable", "cycle": a.cycle, "reason": reason})
        raise SystemExit(f"Unassessable: {reason}") from None
    write_new(a.out, result)


if __name__ == "__main__":
    main()
