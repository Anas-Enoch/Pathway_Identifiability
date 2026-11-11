"""Development-only precision diagnostic. Opens exactly two hash-pinned 2011–12 files.
Usage: python3 -I -B nhanes_value_transport_dev_gate.py PROJECT_ROOT
Writes aggregate JSON to stdout. No external evaluation entry point.
"""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
from scipy.stats import t

root = Path(sys.argv[1]).resolve()
data = root.parent / "nhanes_data"
files = {
    "fas": (data / "FAS_G.xpt", "5a8a907c9797b00d4f9d4dfbf0613e68d14d833d332f2aa3d9dfbedb3525063b"),
    "demo": (data / "DEMO_G.xpt", "eaf0525d1952626885af3e935415a1f66ad62c18698080e7354789c125af252d"),
    "blocks": (root / "nhanes_fa_blocks_v0.4.json", "e3e1d8ac636afd9690e4b5b00b884e21ccf79dad77993e4a75d8776ef3a14338"),
    "wrapper": (root / "nhanes_fa_run.py", "f85d5235cd2fd609145be9ba2b444cd2c4b394a8b3ed6c761c942e765a730297"),
    "gate": (root / "nhanes_fa/frozen/p5_gate3a.py", "4600e2f3e3bb652affadcc6f6935380d8b9a77a4237e7cb306f056f849478130"),
    "fit": (root / "nhanes_fa_stab_out/dev_fit.json", "e9d894aa2e0b549ce4c02712cbbbcb2fd679ffddbac85ce9539420125dc46669"),
}
for key, (path, expected) in files.items():
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError(f"Checksum mismatch for {key}")
def read(key):
    return json.loads(files[key][0].read_text())
def module(key):
    spec = importlib.util.spec_from_file_location(key, files[key][0])
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

w, g = module("wrapper"), module("gate")
blocks = read("blocks")["blocks"]
fit = read("fit")
columns = fit["allvars"]
d = w.load_cycle(files["fas"][0], files["demo"][0], {"O": columns, "P": [], "T": []})
eligible = d.eligible.to_numpy()
e = d.loc[eligible]
assert (len(d), len(e)) == (5643, 1725)
fold = g.folds_for(e.SEQN.astype("int64").astype(str).tolist(), w.SEED_DEV)
weight = np.where(eligible, d.WTFAS2YR.to_numpy(float), 0)
strata, psu = d.SDMVSTRA.to_numpy(), d.SDMVPSU.to_numpy()

def stat(num, den):
    y, b = np.zeros(len(d)), np.zeros(len(d))
    y[eligible], b[eligible] = num, den
    value, se, df, nh, npsu = w.ratio_linearised(y, b, weight, strata, psu)
    assert (df, nh, npsu) == (17, 14, 31)
    hw = float(t.ppf(.975, df)*se)
    return {"estimate": value, "se": se, "ci95": [value-hw, value+hw], "half_width": hw}

def matched(train, test, o, tcount):
    mu, sd = train.mean(0), train.std(0, ddof=1)
    assert np.min(sd) > 1e-12
    Z = (train-mu)/sd
    dtest = (test-mu)/sd
    S = np.cov(Z, rowvar=False, ddof=1)
    alpha = max(.05, len(S)/(len(S)+len(train)))
    S = (1-alpha)*S + alpha*np.trace(S)/len(S)*np.eye(len(S))
    I, T, P = slice(0,o), slice(o,o+tcount), slice(o+tcount,None)
    B = np.linalg.solve(S[I,I]+1e-6*np.eye(o), S[I,:])
    C = S[T,P]-S[T,I]@B[:,P]
    V = np.maximum(np.diag(S[P,P]-S[P,I]@B[:,P]),1e-12)
    score = (C*C).sum(0)/V
    R0 = np.trace(S[T,T]-S[T,I]@B[:,T])
    target_residual = dtest[:,T]-dtest[:,I]@B[:,T]
    candidate_residual = dtest[:,P]-dtest[:,I]@B[:,P]
    base = (target_residual**2).sum(1)
    augmented = ((target_residual[:,:,None]-candidate_residual[:,None,:]*(C/V)[None,:,:])**2).sum(1)
    return score, R0, base, augmented

out, bases, selected, selected_predictions = [], [], [], []
for i, block in enumerate(blocks):
    keys = block["O"]+block["T"]+block["P"]
    o, tcount = len(block["O"]), len(block["T"])
    X = e[keys].to_numpy(float)
    scores, scale, _, _ = matched(X, X[:0], o, tcount)
    np.testing.assert_allclose(scores, fit["full_fit_scores"][i], atol=1e-13, rtol=1e-12)
    order = np.argsort(-scores, kind="stable")
    assert order[0] == fit["full_fit_picks"][i]
    base, loss = np.empty(len(X)), np.empty((len(X),len(order)))
    predicted = np.empty_like(loss)
    for k in range(5):
        take = fold == k
        sc, r0, b, l = matched(X[~take],X[take],o,tcount)
        train = X[~take]
        assert np.argmax(sc) == g.g1_pick(train[:,:o], train[:,o:o+tcount], train[:,o+tcount:])
        # Recorded G1 choices agree in all folds and with the full-development choice.
        assert np.argmax(sc) == order[0]
        base[take], loss[take], predicted[take] = b, l, sc/r0
    assert np.isfinite(loss).all() and np.isfinite(base).all()
    bets = block["P"]
    gains = [{"analyte": bets[j], "score_over_schur_risk": float(scores[j]/scale),
              "gain": stat(base-loss[:,j], base),
              "fold_matched_calibration_error": stat((base-loss[:,j])-predicted[:,j]*base,base)}
             for j in order]
    leader = gains[0]["gain"]["estimate"]
    bar = .1*leader*(1-leader)
    assert bar>0
    def pair(hi,lo):
        r = stat(loss[:,lo]-loss[:,hi],base)
        return {"higher_development_rank":bets[hi], "lower_development_rank":bets[lo], **r,
                "ci_entirely_above_bar":r["ci95"][0]>bar,
                "ci_inside_symmetric_bar":r["ci95"][0]>-bar and r["ci95"][1]<bar}
    adj = [pair(hi,lo) for hi,lo in zip(order[:-1],order[1:])]
    vs = [pair(order[0],j) for j in order[1:]]
    out.append({"block":i+1, "reference_bar_plugin":bar, "gains_ranked":gains,
                "adjacent":adj, "winner_vs_others":vs})
    bases.append(base)
    selected.append(loss[:,order[0]])
    selected_predictions.append(predicted[:,order[0]] * base)
base, win = np.sum(bases,axis=0),np.sum(selected,axis=0)
print(json.dumps({
    "scope":"2011–12 exploratory same-cycle cross-fit feasibility. Intervals conditional on fits, pointwise. No protocol verdict.",
    "input_sha256": {key: item[1] for key,item in files.items()},
    "method":"Matched shrunk affine predictions, training-fold scaling and means, fixed development-ranked pairs, full design Taylor linearisation. Pooled base risk weights the blocks through the pooled denominator.",
    "design":{"eligible":1725,"strata":14,"psus":31,"df":17},
    "pooled_selected_gain":stat(base-win,base),
    "pooled_selected_score_prediction":stat(np.sum(selected_predictions,axis=0),base)["estimate"],
    "pooled_selected_calibration_error":stat(base-win-np.sum(selected_predictions,axis=0),base),
    "blocks":out},indent=2,allow_nan=False))
