"""NSCLC outcome-free evaluation runner (institutional inputs restricted).

Core computations are imported from the shared direct-versus-separate module:
conditional unpenalized Cox censoring with H=Z=(q0,q1), Breslow baseline,
G floor=.02, IPCW logistic m_S, independent calibration/target pairs bootstrap,
source nuisance refit per draw, shared draws for direct/separate endpoints, and
method-specific finite-gamma-grid sup-t critical values.
"""
from __future__ import annotations
import argparse, hashlib, importlib.util, json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np, pandas as pd

HERE=Path(__file__).resolve().parent
REPOSITORY=HERE.parents[1]
CORE=HERE.parent/"inference"/"direct_vs_separate.py"
spec=importlib.util.spec_from_file_location("direct_vs_separate",CORE)
c4=importlib.util.module_from_spec(spec)
spec.loader.exec_module(c4)
FORBIDDEN={"time","event","status","death","recur","rfstime","dtime","rtime"}
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def fixture(root):
    rng=np.random.default_rng(20260918)
    d=root/"pet_radiomic_candidate"/"seed42"
    d.mkdir(parents=True,exist_ok=True)
    for site,nt in (("Stanford",60),("VA",80)):
        n=114; z=rng.normal(size=(n,2)); q0=1/(1+np.exp(-(.4+.7*z[:,0]))); q1=np.clip(q0-.03+.02*z[:,1],.02,.98)
        tf=rng.exponential(1100,n); tc=rng.exponential(900,n); t=np.minimum(tf,tc); e=(tf<=tc).astype(int)
        pd.DataFrame({"patient_id":[f"C{i}" for i in range(n)],"time":t,"event":e,
          "q0_surv2y":q0,"q1_surv2y":q1}).to_csv(d/f"calibration_{site}.csv",index=False)
        z=rng.normal(size=(nt,2)); q0=1/(1+np.exp(-(.35+.7*z[:,0]))); q1=np.clip(q0-.035+.02*z[:,1],.02,.98)
        pd.DataFrame({"patient_id":[f"T{i}" for i in range(nt)],"q0_surv2y":q0,
          "q1_surv2y":q1}).to_csv(d/f"target_{site}.csv",index=False)

def validate(cal,target,cp,tp):
    req={"patient_id","time","event","q0_surv2y","q1_surv2y"}; treq={"patient_id","q0_surv2y","q1_surv2y"}
    if set(cal)!=req: raise ValueError(f"calibration schema mismatch {cp}: {list(cal)}")
    leaked=FORBIDDEN.intersection(target.columns)
    if leaked: raise RuntimeError(f"TARGET OUTCOME FIREWALL: forbidden columns {sorted(leaked)} in {tp}")
    if set(target)!=treq: raise ValueError(f"target schema mismatch {tp}: {list(target)}")
    if not set(cal.event.unique()).issubset({0,1}) or (cal.time<=0).any(): raise ValueError("invalid time/event")
    for frame in (cal,target):
        if not ((frame.q0_surv2y>0)&(frame.q0_surv2y<1)&(frame.q1_surv2y>0)&(frame.q1_surv2y<1)).all():
            raise ValueError("q0/q1 must be in (0,1)")

def discover(root):
    for cp in sorted(root.glob("*/seed*/calibration_*.csv")):
        model=cp.parent.parent.name; seed=int(cp.parent.name.replace("seed","")); site=cp.stem.replace("calibration_","")
        tp=cp.with_name(f"target_{site}.csv")
        if not tp.exists(): raise FileNotFoundError(tp)
        cal=pd.read_csv(cp); target=pd.read_csv(tp); validate(cal,target,cp,tp)
        yield model,site,seed,cp,tp,cal.rename(columns={"q0_surv2y":"q0","q1_surv2y":"q1"}),target.rename(columns={"q0_surv2y":"q0","q1_surv2y":"q1"})

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--input-root",type=Path)
    ap.add_argument("--synthetic",action="store_true"); ap.add_argument("--write-fixture",type=Path)
    ap.add_argument("--bootstrap",type=int,default=999)
    ap.add_argument("--out",type=Path,default=REPOSITORY/"_generated"/"nsclc")
    a=ap.parse_args()
    if a.write_fixture: fixture(a.write_fixture); print("[PASS] fixture",a.write_fixture); return
    root=(REPOSITORY/"synthetic_fixture"/"nsclc") if a.synthetic else a.input_root
    if root is None or not root.is_dir():
        raise FileNotFoundError("Restricted NSCLC inputs are not distributed. Supply --input-root conforming to code/nsclc/input_schema.md, or use --synthetic.")
    rows=[]; meta=[]; inputs=[]
    for model,site,seed,cp,tp,cal,target in discover(root):
        rr,mm=c4.analyze(f"NSCLC-{site}",model,seed,cal,target,730.0,a.bootstrap)
        rows.extend(rr); meta.append(mm); inputs.extend([{"logical_role":"calibration","sha256":sha(cp)},{"logical_role":"target_outcome_free","sha256":sha(tp)}])
    if not rows: raise RuntimeError("No complete cells found; see code/nsclc/input_schema.md")
    a.out.mkdir(parents=True,exist_ok=True); out=a.out/"gamma_sweep_direct_separate.csv"; pd.DataFrame(rows).to_csv(out,index=False)
    manifest={"status":"PASS_SYNTHETIC_CODE_PATH" if a.synthetic else "PASS_RESTRICTED_INPUT_ANALYSIS",
      "created_utc":datetime.now(timezone.utc).isoformat(),"synthetic":a.synthetic,"target_outcomes_read":False,
      "core_script_sha256":sha(CORE),"runner_sha256":sha(__file__),"bootstrap":a.bootstrap,
      "gamma_grid":c4.GAMMA_GRID.tolist(),"g_floor":c4.G_MIN,"inputs":inputs,"output_sha256":sha(out),"job_metadata":meta}
    (a.out/"MANIFEST.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    print("[PASS]",manifest["status"],"rows",len(rows),"->",out)
if __name__=="__main__": main()
