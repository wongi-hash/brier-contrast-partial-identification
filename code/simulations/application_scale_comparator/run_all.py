"""Portable orchestrator for the locked nT=60/80 comparator study.

Smoke mode validates both ordinary and candidate-superior mirror generators.
Full mode regenerates all 18 locked cells and the 288-row aggregate.  Existing
checkpoints are resumed.  No model, gamma, seed, or scenario is tuned here.
"""
from __future__ import annotations
import argparse, os, subprocess, sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
COMP=("no_shift","constant_shift","xvary_in_class","sign_change","out_of_class","low_overlap","pS_misspec")
MIRROR=("no_shift","xvary_in_class")

def run(script, env):
    shown={k:env[k] for k in sorted(env) if k.startswith("SIMULATION_")}
    print("[RUN]",script.name,shown,flush=True)
    subprocess.run([sys.executable,"-u",str(script)],cwd=HERE,env=env,check=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--smoke",action="store_true",help="2 reps, 3 bootstrap draws, nT=60 only")
    ap.add_argument("--full",action="store_true",help="locked 1000 reps x 199 bootstrap, nT=60/80")
    ap.add_argument("--output",type=Path,default=HERE/"_regenerated")
    args=ap.parse_args()
    if args.smoke==args.full: ap.error("choose exactly one of --smoke or --full")
    base=os.environ.copy(); base.update(SIMULATION_SOURCE_SIZE="114",SIMULATION_WITH_CI="1",SIMULATION_SEED="42")
    if args.smoke:
        base.update(SIMULATION_REPLICATES="2",SIMULATION_BOOTSTRAP="3"); targets=(60,); comp=("no_shift",); mirror=("no_shift",)
    else:
        base.update(SIMULATION_REPLICATES="1000",SIMULATION_BOOTSTRAP="199"); targets=(60,80); comp=COMP; mirror=MIRROR
    comp_root=args.output/"comparator"; mirror_root=args.output/"candidate_superior_mirror"
    for nt in targets:
        for scenario in comp:
            env=base.copy(); env.update(SIMULATION_TARGET_SIZE=str(nt),SIMULATION_SCENARIO=scenario,SIMULATION_OUTPUT_ROOT=str(comp_root))
            run(HERE/"run_comparator.py",env)
        for scenario in mirror:
            env=base.copy(); env.update(SIMULATION_TARGET_SIZE=str(nt),SIMULATION_SCENARIO=scenario,SIMULATION_OUTPUT_ROOT=str(mirror_root))
            run(HERE/"run_candidate_superior_mirror.py",env)
    if args.full:
        env=base.copy(); env.update(SIMULATION_COMPARATOR_ROOT=str(comp_root),SIMULATION_MIRROR_ROOT=str(mirror_root),SIMULATION_AGGREGATE_ROOT=str(args.output))
        run(HERE/"aggregate_results.py",env)
        import pandas as pd
        out=pd.read_csv(args.output/"application_scale_comparator.csv")
        if len(out)!=288: raise AssertionError(f"expected 288 aggregate rows, got {len(out)}")
    print("[PASS] comparator build",args.output)
if __name__=="__main__": main()
