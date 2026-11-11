"""Development-only rank diagnostic for the matched affine instrument. Opens exactly two hash-pinned 2011–12 files.
Usage: python3 -I -B nhanes_value_transport_dev_gate.py PROJECT_ROOT
Writes aggregate JSON to stdout. No external evaluation entry point.
"""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
from scipy.stats import rankdata, t

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


pred, real, draws = [], [], []
for i, block in enumerate(blocks):
    cols = block["O"]+block["T"]+block["P"]
    o, tcount = len(block["O"]), len(block["T"])
    X = e[cols].to_numpy(float)
    scores, scale, _, _ = matched(X, X[:0], o, tcount)
    np.testing.assert_allclose(scores, fit["full_fit_scores"][i], atol=1e-13, rtol=1e-12)
    base, loss = np.empty(len(X)), np.empty((len(X),len(scores)))
    for k in range(5):
        te = fold==k
        sc, _, base[te], loss[te] = matched(X[~te],X[te],o,tcount)
        assert np.argmax(sc) == fit["full_fit_picks"][i]
    num, den = sums(base[:,None]-loss), sums(base[:,None])
    pred.append(scores/scale)
    real.append(num.sum(axis=0)/den.sum())
    draws.append((boot_weights@num)/(boot_weights@den))
r = rank_corr(np.concatenate(pred),np.concatenate(draws,axis=1))
m = np.mean([rank_corr(p,b) for p,b in zip(pred,draws)],axis=0)
print(json.dumps({
    "scope":"NHANES 2011–12 only. Matched affine instrument, development-only.",
    "method":"5,000 stratified rescaled PSU bootstrap draws, conditional on predictions, score order, and folds.",
    "seed":2026111002,
    "pooled_spearman":float(rank_corr(np.concatenate(pred),np.concatenate(real))),
    "pooled_ci95":np.quantile(r,[.025,.975]).tolist(),
    "pooled_slope_realized_on_score":float(np.polyfit(np.concatenate(pred),np.concatenate(real),1)[0]),
    "macro_within_block_spearman":float(np.mean([rank_corr(p,q) for p,q in zip(pred,real)])),
    "macro_within_block_ci95":np.quantile(m,[.025,.975]).tolist(),
    "within_block_spearman":[float(rank_corr(p,q)) for p,q in zip(pred,real)]},indent=2,allow_nan=False))
