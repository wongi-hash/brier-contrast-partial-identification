"""Candidate-superior mirror arms for the locked comparator study.

Only q0/q1 roles are exchanged; the outcome law and inference engine remain
unchanged. Environment variables are documented in run_all.py.
"""
from __future__ import annotations
import hashlib, importlib.util, json, os, platform, sys
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
SCEN = os.environ.get("SIMULATION_SCENARIO", "no_shift")
if SCEN not in ("no_shift", "xvary_in_class"):
    raise ValueError(SCEN)
spec = importlib.util.spec_from_file_location("comparator", HERE / "run_comparator.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
original_q0, original_q1 = m.q0f, m.q1f
def source_lp_original(x): return -0.25 + 1.40 * original_q0(x) - 0.80 * original_q1(x)
def ps_original(x, scenario): return m.np.clip(m.expit(source_lp_original(x)), 1e-5, 1-1e-5)
m.q0f, m.q1f = original_q1, original_q0
m.ps_linear_predictor, m.ps_true = source_lp_original, ps_original

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""): h.update(b)
    return h.hexdigest()

def main():
    cfg={"study":"candidate-superior mirror","scenario":SCEN,"n_source":m.N_SOURCE,
         "n_target":m.N_TARGET,"reps":m.REPS,"bootstrap":m.BOOT,"with_ci":m.WITH_CI,
         "gammas":m.GAMMAS.tolist(),"base_seed":m.BASE_SEED,
         "source_censoring_design":m.CENSOR_SOURCE_DESIGN,"censor_scale":m.CENSOR_SCALE,
         "construction":"swap q0/q1 only; preserve original pS/pT laws"}
    tag=hashlib.sha256(json.dumps(cfg,sort_keys=True).encode()).hexdigest()[:12]
    root=Path(os.environ.get("SIMULATION_OUTPUT_ROOT",str(HERE/"candidate_superior_mirror")))
    out=root/f"run_{tag}"; out.mkdir(parents=True,exist_ok=True)
    checkpoint=out/f"replicates_mirror_{SCEN}.csv"; completed=set(); blocks=[]
    if checkpoint.exists():
        x=pd.read_csv(checkpoint); completed=set(x.rep.astype(int)); blocks.append(x)
    truth=m.population_truth(SCEN)
    if truth>=0: raise AssertionError("mirror arm must be candidate-superior")
    for rep in range(m.REPS):
        if rep in completed: continue
        b=pd.DataFrame(m.one_rep(rep,SCEN,truth)); b["scenario"]="mirror_"+SCEN
        b.to_csv(checkpoint,mode="a",index=False,header=not checkpoint.exists()); blocks.append(b)
    frame=pd.concat(blocks,ignore_index=True).drop_duplicates(["scenario","rep","method","gamma"],keep="last")
    summary=m.summarize(frame); rp=out/"all_replicates.csv"; sp=out/"summary.csv"
    frame.to_csv(rp,index=False); summary.to_csv(sp,index=False)
    manifest={"status":"PASS","created_utc":datetime.now(timezone.utc).isoformat(),"config":cfg,
      "truth":truth,"code_sha256":sha(__file__),"engine_sha256":sha(HERE/"run_comparator.py"),
      "python":sys.version,"platform":platform.platform(),"outputs":{"replicates":sha(rp),"summary":sha(sp)}}
    (out/"MANIFEST.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    print("[PASS]",out)
if __name__=="__main__": main()
