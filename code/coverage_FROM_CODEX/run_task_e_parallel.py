"""Run the eight Task E shards, then aggregate (portable public launcher)."""
from __future__ import annotations
import json, os, subprocess, sys, time
from datetime import datetime, timezone
from pathlib import Path

HERE=Path(__file__).resolve().parent
PYTHON=Path(sys.executable)
RUNNER=HERE/"RUN_TASK_E_PRIMARY_COVERAGE.py"
LOG=HERE/"task_e_logs"; STATUS=HERE/"TASK_E_STATUS.json"

def save(x):
    x["updated_utc"]=datetime.now(timezone.utc).isoformat(); STATUS.write_text(json.dumps(x,indent=2),encoding="utf-8")

def main():
    LOG.mkdir(exist_ok=True); state={"status":"RUNNING","completed":0,"total":8};save(state)
    jobs=[(nt,seed) for nt in (60,80) for seed in (42,52,62,72)]; running={}; pending=list(jobs);done=[]
    while pending or running:
        while pending and len(running)<8:
            nt,seed=pending.pop(0);tag=f"nt{nt}_seed{seed}";env=os.environ.copy();env.update({"TASK_E_NT":str(nt),"TASK_E_SEED":str(seed)})
            so=(LOG/f"{tag}.stdout.log").open("w",encoding="utf-8");se=(LOG/f"{tag}.stderr.log").open("w",encoding="utf-8")
            p=subprocess.Popen([str(PYTHON),"-u",str(RUNNER)],cwd=str(HERE),env=env,stdout=so,stderr=se,
                               creationflags=getattr(subprocess,"CREATE_NO_WINDOW",0));running[p]=(nt,seed,so,se)
        time.sleep(3)
        for p,v in list(running.items()):
            rc=p.poll()
            if rc is None:continue
            nt,seed,so,se=v;so.close();se.close();done.append({"n_target":nt,"seed":seed,"returncode":rc});del running[p]
            state["completed"]=len(done);state["failed"]=sum(x["returncode"]!=0 for x in done);state["recent"]=done[-4:];save(state)
    if any(x["returncode"]!=0 for x in done):state["status"]="FAILED_SHARDS";save(state);return 1
    env=os.environ.copy();env["TASK_E_AGG_ONLY"]="1"
    with (LOG/"aggregate.stdout.log").open("w",encoding="utf-8") as so,(LOG/"aggregate.stderr.log").open("w",encoding="utf-8") as se:
        rc=subprocess.run([str(PYTHON),"-u",str(RUNNER)],cwd=str(HERE),env=env,stdout=so,stderr=se,
                          creationflags=getattr(subprocess,"CREATE_NO_WINDOW",0)).returncode
    state["status"]="PASS" if rc==0 else "FAILED_AGGREGATION";state["aggregation_returncode"]=rc;save(state);return rc

if __name__=="__main__":raise SystemExit(main())
