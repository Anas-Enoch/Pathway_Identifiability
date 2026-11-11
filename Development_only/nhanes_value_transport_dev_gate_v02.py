"""NHANES 2011–12 only. Reproduce frozen E1, and estimate conditional E2 rank uncertainty.
Usage: python3 -I -B nhanes_value_transport_dev_gate_v02.py PROJECT_ROOT [WORKERS]
All inputs are hash-pinned. No participant output or external-cycle entry point.
"""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
from scipy.stats import nct, rankdata, t

root = Path(sys.argv[1]).resolve()
workers = int(sys.argv[2]) if len(sys.argv) > 2 else 5
if not 1 <= workers <= 5:
    raise ValueError("Use one to five workers.")
data = root.parent / "nhanes_data"
inputs = {
    "fas": (data / "FAS_G.xpt", "5a8a907c9797b00d4f9d4dfbf0613e68d14d833d332f2aa3d9dfbedb3525063b"),
    "demo": (data / "DEMO_G.xpt", "eaf0525d1952626885af3e935415a1f66ad62c18698080e7354789c125af252d"),
    "blocks": (root / "nhanes_fa_blocks_v0.4.json", "e3e1d8ac636afd9690e4b5b00b884e21ccf79dad77993e4a75d8776ef3a14338"),
    "wrapper": (root / "nhanes_fa_run.py", "f85d5235cd2fd609145be9ba2b444cd2c4b394a8b3ed6c761c942e765a730297"),
    "fit": (root / "nhanes_fa_stab_out/dev_fit.json", "e9d894aa2e0b549ce4c02712cbbbcb2fd679ffddbac85ce9539420125dc46669"),
    "record": (root / "nhanes_fa_out/FAS_G_development_v04_result.json", "f7f9f5cad1c25edcd9e395ae4198b5997880d80c58552fee4aa77b91121c2101"),
}
for key, (path, expected) in inputs.items():
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError(f"Checksum mismatch for {key}")

def read(key):
    return json.loads(inputs[key][0].read_text())

spec = importlib.util.spec_from_file_location("frozen_wrapper", inputs["wrapper"][0])
w = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w)
n2 = w.frozen()
g = n2.gate3a()
blocks = read("blocks")["blocks"]
fit = read("fit")
record = read("record")["inference"]
columns = fit["allvars"]
d = w.load_cycle(inputs["fas"][0], inputs["demo"][0], {"O": columns, "P": [], "T": []})
eligible = d.eligible.to_numpy()
assert (len(d), int(eligible.sum())) == (5643, 1725)
weight = np.where(eligible, d.WTFAS2YR.to_numpy(float), 0)
strata, psu = d.SDMVSTRA.to_numpy(), d.SDMVPSU.to_numpy()

def stat(num, den):
    y, b = np.zeros(len(d)), np.zeros(len(d))
    y[eligible], b[eligible] = num, den
    value, se, df, nh, npsu = w.ratio_linearised(y, b, weight, strata, psu)
    assert (df, nh, npsu) == (17, 14, 31)
    hw = float(t.ppf(.975, df)*se)
    return {"estimate": value, "se": se, "ci95": [value-hw, value+hw], "half_width": hw}

def block_losses(job):
    i, block = job
    O, X, Y, subjects, _ = w.arrays(d, block)
    fold = g.folds_for(subjects, w.SEED_DEV)
    losses = {k: np.full((len(Y), X.shape[1]+1), np.nan) for k in ("sq","bd")}
    picks = []
    for k in range(n2.K_OUT):
        te = fold == k
        tr = ~te
        jg = g.g1_pick(O[tr], Y[tr], X[tr])
        picks.append(jg)
        assert jg == record["picks"][i*n2.K_OUT+k]["G1"] == fit["full_fit_picks"][i]
        _, best = n2.inner_scores("rbf", O[tr], X[tr], Y[tr], n2.SEED+k)
        for j in range(X.shape[1]+1):
            args = (X[tr,j-1], X[te,j-1]) if j else ()
            pred = n2.predict_grid("rbf", O[tr], O[te], Y[tr], *args)[best[j]]
            err = (pred-Y[te])**2
            losses["sq"][te,j] = err.sum(1)
            losses["bd"][te,j] = np.minimum(err,n2.CAP).sum(1)
        print(f"block {i+1}/5, fold {k+1}/5 complete", file=sys.stderr, flush=True)
    assert all(np.isfinite(x).all() for x in losses.values())
    return losses, picks[0]

with ThreadPoolExecutor(max_workers=workers) as pool:
    fitted = list(pool.map(block_losses, enumerate(blocks)))

e1 = {}
for k in ("sq","bd"):
    base = sum(loss[k][:,0] for loss,j in fitted)
    selected = sum(loss[k][:,j+1] for loss,j in fitted)
    s = stat(base-selected,base)
    np.testing.assert_allclose([s["estimate"],s["se"]],
                               [record[k]["delta_g1_base"],record[k]["se_g1_base"]],rtol=1e-7,atol=1e-10)
    e1[k] = s

# Use all participants in each PSU together. The same replicate weights are used
# for all 45 contrasts. Draw n_h-1 PSUs per stratum with replacement and rescale.
psu_groups = []
for h in np.unique(strata):
    psu_groups.append([np.flatnonzero(((strata==h)&(psu==j))[eligible]) for j in np.unique(psu[strata==h])])
assert len(psu_groups)==14 and sum(map(len,psu_groups))==31 and min(map(len,psu_groups))>=2
NREP = 5000
rng = np.random.default_rng(2026111002)
boot_weights = np.empty((NREP,31))
offset = 0
for group in psu_groups:
    n = len(group)
    counts = rng.multinomial(n-1,np.repeat(1/n,n),size=NREP)
    boot_weights[:,offset:offset+n] = counts * n/(n-1)
    offset += n
groups = [indices for group in psu_groups for indices in group]

def sums(mat):
    weighted = mat*weight[eligible,None]
    return np.stack([weighted[index].sum(axis=0) for index in groups])

def rank_corr(x,y):
    u = rankdata(x)
    v = rankdata(y,axis=-1)
    a = u-u.mean()
    b = v-v.mean(axis=-1,keepdims=True)
    return (b@a)/np.sqrt((a*a).sum()*(b*b).sum(axis=-1))

all_pred, all_real, all_draws, by_block = [],[],[],[]
e3 = []
for i, ((loss,j), block) in enumerate(zip(fitted,blocks)):
    base = loss["sq"][:,0]
    aug = loss["sq"][:,1:]
    scores = np.asarray(fit["full_fit_scores"][i])
    order = np.argsort(-scores,kind="stable")
    # Normalization preserves within-block rank. Use matched Schur base risk
    # for the pooled rank comparison, as in the affine diagnostic.
    cov = np.asarray(fit["cov_ddof1"])
    cov = cov/np.sqrt(np.outer(cov.diagonal(),cov.diagonal()))
    alpha = max(.05,len(cov)/(len(cov)+fit["n_eligible"]))
    cov = (1-alpha)*cov+alpha*np.eye(len(cov))
    idx = {v:j for j,v in enumerate(columns)}
    oi, ti = ([idx[v] for v in block[k]] for k in ("O","T"))
    oo = cov[np.ix_(oi,oi)]+1e-6*np.eye(len(oi))
    ot = cov[np.ix_(oi,ti)]
    schur = np.trace(cov[np.ix_(ti,ti)]-ot.T@np.linalg.solve(oo,ot))
    score = scores/schur
    num, den = sums(base[:,None]-aug), sums(base[:,None])
    realized = num.sum(axis=0)/den.sum()
    replicate = (boot_weights@num)/(boot_weights@den)
    by_block.append({"block": i+1, "spearman": float(rank_corr(score,realized)),
        "leader":block["P"][j], "winner_gain":stat(base-aug[:,j],base),
        "slope_realized_on_score":float(np.polyfit(score,realized,1)[0]),
        "candidate_gains":[{"analyte":v,"score":float(s),"realized":float(r)} for v,s,r in zip(block["P"],score,realized)]})
    all_pred.extend(score); all_real.extend(realized); all_draws.append(replicate)
    bar = .10*float(realized[j])*(1-float(realized[j]))
    comparisons = {}
    for name, pairs in (("adjacent",zip(order[:-1],order[1:])),("winner_vs_others",((j,k) for k in order[1:]))):
        comparisons[name] = [{"higher_rank":block["P"][a],"lower_rank":block["P"][b],
                              **stat(aug[:,b]-aug[:,a],base)} for a,b in pairs]
    e3.append({"block":i+1,"development_reference":bar,**comparisons})

pooled_draws = rank_corr(np.array(all_pred),np.concatenate(all_draws,axis=1))
macro_draws = np.mean([rank_corr(np.asarray(all_pred)[9*i:9*(i+1)],r) for i,r in enumerate(all_draws)],axis=0)
assert np.isfinite(pooled_draws).all() and np.isfinite(macro_draws).all()
estimate = e1["sq"]["estimate"]
se = e1["sq"]["se"]
df = 17
cut = float(t.ppf(.975,df))
floor = .75*estimate
print(json.dumps({
    "scope":"2011–12 only. E1 exact re-execution of the frozen G1 and base predictions. E2/E3 are development diagnostics using all nine candidates.",
    "inputs_sha256":{k:v[1] for k,v in inputs.items()},
    "design":{"eligible":1725,"strata":14,"psus":31,"df":17},
    "E1":e1,
    "E1_proposed_minimum_retained_gain":{"retention_fraction":.75,"frozen_dev_gain":estimate,
       "floor":floor,"one_sided_alpha":.025,"assumed_df":df,
       "conditional_power_if_true_gain_equals_dev":{str(m):float(nct.sf(cut,df,(estimate-floor)/(se*m))) for m in (1,1.25,1.5,2)}},
    "E2":{"method":"5,000 stratified rescaled PSU bootstrap draws; shared participant weights across all blocks/candidates. Percentile interval conditional on predictions, score order, and folds.",
          "seed":2026111002,"pooled_spearman":float(rank_corr(np.array(all_pred),np.array(all_real))),
          "pooled_ci95":np.quantile(pooled_draws,[.025,.975]).tolist(),
          "pooled_slope_realized_on_score":float(np.polyfit(all_pred,all_real,1)[0]),
          "macro_within_block_spearman":float(np.mean([b["spearman"] for b in by_block])),
          "macro_within_block_ci95":np.quantile(macro_draws,[.025,.975]).tolist(),"blocks":by_block},
    "E3":{"status":"descriptive only","blocks":e3},
},indent=2,allow_nan=False))
