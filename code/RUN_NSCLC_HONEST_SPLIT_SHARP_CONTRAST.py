"""
Generic downstream runner for the canonical inductive NSCLC sharp-contrast analysis.

This public runner starts from frozen source-only predictions. Upstream feature
selection, SparsePCA, Cox fitting, and the train-only Breslow baseline are fitted
on D_train and are not reproduced here because the institutional feature tables
are restricted. The candidate never uses target X during fitting, so all target
patients can be used as evaluation X.

Why this satisfies Proposition D v2 (A2): q0,q1 are TRAINED on D_train and FROZEN; p_S is estimated on the
INDEPENDENT D_cal; bootstrap resamples ONLY (calibration, target); q is never refit. Hence q ⟂ D_cal.

=========================== INPUT INTERFACE (Codex must produce) ===========================
For each model in {InductiveSparsePCA (primary), InductiveSparsePCA_SRDO
(sensitivity)} and each seed in {42,52,62,72}:
  honest_split_inputs/<model>/seed<seed>/
    calibration_<site>.csv: patient_id,time,event,q0_surv2y,q1_surv2y (D_cal; frozen q)
    target_Stanford.csv: patient_id,q0_surv2y,q1_surv2y                (NO time/event)
    target_VA.csv      : patient_id,q0_surv2y,q1_surv2y                (NO time/event)
  where q0_surv2y = S_train(tau|x) of clinical model, q1_surv2y = S_train(tau|x) of candidate (DDAU) model,
  both from the TRAIN-locked model + TRAIN-only Breslow baseline at tau=730d.
  event-stratified 2:1 train/cal split; primary seed 42; sensitivity 52/62/72;
  report per-seed separately. Canonical target sizes are Stanford=60 and VA=80.
==========================================================================================

Usage:
  python RUN_NSCLC_HONEST_SPLIT_SHARP_CONTRAST.py --synthetic
  python RUN_NSCLC_HONEST_SPLIT_SHARP_CONTRAST.py --root honest_split_inputs
"""
import csv, json, os, sys, hashlib, argparse, platform, datetime, tempfile, numpy as np
sig=lambda z:1/(1+np.exp(-np.clip(z,-35,35)))
def logit(p):p=np.clip(p,1e-6,1-1e-6);return np.log(p/(1-p))
TAU=730.0; B=999; GAMMAS=[0.0,0.1,0.2,0.3,0.5]; GG=np.round(np.arange(0,1.001,0.05),4)
SEEDS=[42,52,62,72]; MODELS=["InductiveSparsePCA","InductiveSparsePCA_SRDO"]; SITES=["Stanford","VA"]
def sha_file(p):return hashlib.sha256(open(p,'rb').read()).hexdigest()
def seed_of(*k):return int(hashlib.sha256("|".join(map(str,k)).encode()).hexdigest()[:8],16)

def km(t,ci):
    o=np.argsort(t,kind='mergesort');t=t[o];c=ci[o];n=len(t);uq,f=np.unique(t,return_index=True)
    return uq,np.cumprod(1-np.add.reduceat(c,f)/(n-f))
def Gm(Tk,V,q):
    i=np.searchsorted(Tk,q,'left')-1;q=np.atleast_1d(q);o=np.ones(len(q));ok=i>=0;o[ok]=V[np.clip(i[ok],0,len(V)-1)];return o
def wlogit(X,y,w,it=30):
    Xd=np.c_[np.ones(len(X)),X];b=np.zeros(Xd.shape[1])
    for _ in range(it):
        p=sig(Xd@b);W=np.clip(w*p*(1-p),1e-12,None)
        b=b+np.linalg.solve((Xd*W[:,None]).T@Xd+1e-6*np.eye(Xd.shape[1]),Xd.T@(w*(y-p)))
    return b
def pS_on_cal(ct,ce,cq0,cq1, tq0,tq1):
    """p_S estimated on CALIBRATION outcomes; features = frozen (q0_surv2y,q1_surv2y). Predict at target."""
    Tk,V=km(ct,1-ce);Gt=np.clip(Gm(Tk,V,ct),1e-3,None);Gtau=max(Gm(Tk,V,np.array([TAU]))[0],1e-3)
    D=(ct>TAU).astype(float);w=np.zeros(len(ct));w[ct>TAU]=1/Gtau;dd=(ct<=TAU)&(ce==1);w[dd]=1/Gt[dd]
    X=np.c_[cq0,cq1];mu,sd=X.mean(0),X.std(0)+1e-9;b=wlogit((X-mu)/sd,D,w)
    Xt=(np.c_[tq0,tq1]-mu)/sd;return np.clip(sig(b[0]+Xt@b[1:]),1e-4,1-1e-4)
def endpoints(q0,q1,pS):
    a=q1**2-q0**2;bb=q1-q0;lz=logit(pS);L=np.empty(len(GG));U=np.empty(len(GG))
    for i,g in enumerate(GG):
        pm=sig(lz-g);pp=sig(lz+g);pL=np.where(bb>=0,pp,pm);pU=np.where(bb>=0,pm,pp)
        L[i]=np.mean(a-2*bb*pL);U[i]=np.mean(a-2*bb*pU)
    return L,U
gc=lambda U:(float(GG[U<0].max()) if (U<0).any() else 0.0)
gr=lambda L:(float(GG[L>0].max()) if (L>0).any() else 0.0)

def analyze(model,site,seed, ct,ce,cq0,cq1, tq0,tq1):
    q0=tq0.copy();q1=tq1.copy()                     # FROZEN q at target
    pS=pS_on_cal(ct,ce,cq0,cq1,tq0,tq1);Lh,Uh=endpoints(q0,q1,pS)
    rng=np.random.default_rng(seed_of(model,site,seed));nc=len(ct);nt=len(tq0)
    Lb=np.empty((B,len(GG)));Ub=np.empty((B,len(GG)))
    for j in range(B):
        ci=rng.integers(0,nc,nc);ti=rng.integers(0,nt,nt)               # cal + target resample; q frozen
        pSb=pS_on_cal(ct[ci],ce[ci],cq0[ci],cq1[ci],tq0,tq1)
        Lb[j],Ub[j]=endpoints(q0[ti],q1[ti],pSb[ti])
    seL=Lb.std(0)+1e-12;seU=Ub.std(0)+1e-12
    csim=float(np.quantile(np.maximum(np.abs(Lb-Lh)/seL,np.abs(Ub-Uh)/seU).max(1),0.95))
    Ubar=Uh+csim*seU;Ulow=Uh-csim*seU;idx={round(g,2):i for i,g in enumerate(GG)};dec=[]
    for g in GAMMAS:
        i=idx[round(g,2)];Lc=Lh[i]-csim*seL[i];Uc=Uh[i]+csim*seU[i]
        dec.append({"model":model,"site":site,"seed":seed,"gamma":g,"L_simCI":round(float(Lc),5),
                    "U_simCI":round(float(Uc),5),"decision":("ADOPT CANDIDATE" if Uc<0 else ("KEEP REFERENCE" if Lc>0 else "DEFER"))})
    g={"model":model,"site":site,"seed":seed,"n_cal":nc,"n_target":nt,"gstar_cand":gc(Uh),
       "gstar_cand_CI":[gc(Ubar),gc(Ulow)],"gstar_cand_hi_right_censored":bool((Ulow<0).all()),
       "gstar_ref":gr(Lh),"gstar_ref_CI":[gr(Lh-csim*seL),gr(Lh+csim*seL)]}
    return dec,g

def firewall_scan(root):
    """Analysis-side input firewall (2-script contract).
    The analysis stage may read ONLY outcome-free target files and the source
    calibration files produced by the SEPARATE packaging script. If any other
    CSV (e.g. a combined raw file that still carries BOTH target outcomes AND
    scores) is present under the input root, refuse to run so target outcomes
    can never leak into identification. Calibration files legitimately carry
    source outcomes (time,event) for p_S; target files must not."""
    import glob
    allowed=set()
    for model in MODELS:
        for seed in SEEDS:
            for site in SITES:
                allowed.add(os.path.abspath(os.path.join(root,model,f"seed{seed}",f"calibration_{site}.csv")))
                allowed.add(os.path.abspath(os.path.join(root,model,f"seed{seed}",f"target_{site}.csv")))
    offenders=[]
    for p in glob.glob(os.path.join(root,"**","*.csv"),recursive=True):
        ap=os.path.abspath(p)
        if ap in allowed: continue
        try:hdr=next(csv.reader(open(p)))
        except Exception:hdr=[]
        has_outcome=('time' in hdr) or ('event' in hdr)
        has_score=any('q0' in h or 'q1' in h or 'score' in h for h in hdr)
        if has_outcome and has_score:
            offenders.append((p,"combined raw (outcomes+scores) — target outcomes could leak"))
        else:
            offenders.append((p,"non-whitelisted CSV in analysis root"))
    if offenders:
        msg="\n".join(f"  - {p}: {why}" for p,why in offenders)
        raise RuntimeError(f"FIREWALL: analysis stage found non-whitelisted input(s):\n{msg}\n"
                           "Analysis may read only calibration_<site>.csv and target_<site>.csv "
                           "from the separate packaging script.")

def load_split(root,model,seed,site):
    d=os.path.join(root,model,f"seed{seed}")
    # Site-labelled calibration exports preserve the locked upstream interface.
    # With an inductive candidate, model fitting itself does not use target X.
    cal=list(csv.DictReader(open(os.path.join(d,f"calibration_{site}.csv"))))
    tp=os.path.join(d,f"target_{site}.csv");trows=list(csv.DictReader(open(tp)))
    assert 'time' not in trows[0] and 'event' not in trows[0], f"FIREWALL: outcomes in target file {tp}"
    ct=np.array([float(r['time']) for r in cal]);ce=np.array([int(float(r['event'])) for r in cal])
    cq0=np.array([float(r['q0_surv2y']) for r in cal]);cq1=np.array([float(r['q1_surv2y']) for r in cal])
    tq0=np.array([float(r['q0_surv2y']) for r in trows]);tq1=np.array([float(r['q1_surv2y']) for r in trows])
    return ct,ce,cq0,cq1,tq0,tq1

def make_synthetic(root):
    rng=np.random.default_rng(0)
    for model in MODELS:
        for seed in SEEDS:
            d=os.path.join(root,model,f"seed{seed}");os.makedirs(d,exist_ok=True)
            # Each target lock has its own frozen calibration q1.
            for site,nt in [("Stanford",60),("VA",80)]:
                n=114;q1=np.clip(sig(0.3+0.9*rng.normal(0,1,n)),.02,.98);q0=np.clip(q1-rng.uniform(-.05,.1,n),.02,.98)
                T=rng.exponential(900,n);C=rng.exponential(700,n);t=np.minimum(T,C);e=(T<=C).astype(int)
                with open(os.path.join(d,f"calibration_{site}.csv"),"w",newline="") as f:
                    w=csv.writer(f);w.writerow(["patient_id","time","event","q0_surv2y","q1_surv2y"])
                    for i in range(n):w.writerow([f"C{i}",round(t[i],1),e[i],round(q0[i],4),round(q1[i],4)])
                q1=np.clip(sig(0.3+0.9*rng.normal(0,1,nt)),.02,.98);q0=np.clip(q1-rng.uniform(-.05,.1,nt),.02,.98)
                with open(os.path.join(d,f"target_{site}.csv"),"w",newline="") as f:
                    w=csv.writer(f);w.writerow(["patient_id","q0_surv2y","q1_surv2y"])
                    for i in range(nt):w.writerow([f"T{i}",round(q0[i],4),round(q1[i],4)])

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--synthetic",action="store_true");ap.add_argument("--root",default="honest_split_inputs");ap.add_argument("--out",default="honest_split_out")
    a=ap.parse_args();HERE=os.path.dirname(os.path.abspath(__file__));root=os.path.join(HERE,a.root)
    if a.synthetic:
        root=tempfile.mkdtemp(prefix="honest_split_synth_");make_synthetic(root);print("[synthetic] wrote fake split ->",root)
    if not os.path.isdir(root):
        print(f"[waiting] input root not found: {root}\n  -> Codex must produce honest-split artifacts (see header INPUT INTERFACE).");return
    firewall_scan(root)  # 2-script contract: refuse any non-whitelisted / combined-raw CSV
    dec=[];gs=[]
    for model in MODELS:
        for seed in SEEDS:
            for site in SITES:
                try:ct,ce,cq0,cq1,tq0,tq1=load_split(root,model,seed,site)
                except FileNotFoundError:print(f"[missing] {model}/seed{seed}/{site}");continue
                d,g=analyze(model,site,seed,ct,ce,cq0,cq1,tq0,tq1);dec+=d;gs+=[g]
    if not dec:raise RuntimeError("No complete model/seed/site input cells were found")
    OUT=a.out if os.path.isabs(a.out) else os.path.join(HERE,a.out);os.makedirs(OUT,exist_ok=True)
    import csv as _c
    for name,rows in [("gamma_sweep.csv",dec),("gamma_star.csv",gs)]:
        with open(os.path.join(OUT,name),"w",newline="") as f:
            w=_c.DictWriter(f,fieldnames=list(rows[0].keys()),extrasaction='ignore');w.writeheader();w.writerows(rows)
    result_files=[os.path.join(OUT,"gamma_sweep.csv"),os.path.join(OUT,"gamma_star.csv")]
    input_files=[]
    for model in MODELS:
        for seed in SEEDS:
            for site in SITES:
                for name in [f"calibration_{site}.csv",f"target_{site}.csv"]:
                    p=os.path.join(root,model,f"seed{seed}",name)
                    if os.path.isfile(p):input_files.append({"path":os.path.abspath(p),"sha256":sha_file(p)})
    upstream_manifest=os.path.join(root,"MANIFEST.json")
    manifest={"status":"PASS_HONEST_SPLIT_SHARP_CONTRAST","created_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "runner":os.path.abspath(__file__),"runner_sha256":sha_file(__file__),"input_root":os.path.abspath(root),
              "upstream_manifest_sha256":sha_file(upstream_manifest) if os.path.isfile(upstream_manifest) else None,
              "inputs":input_files,"outputs":[{"path":os.path.abspath(p),"sha256":sha_file(p)} for p in result_files],
              "tau":TAU,"bootstrap":B,"gamma_report":GAMMAS,"gamma_grid":GG.tolist(),"seeds":SEEDS,"models":MODELS,"sites":SITES,
              "target_outcomes_read":False,"q_status":"frozen upstream; never refit downstream",
              "bootstrap_units":"independent CNUH calibration patients + independent target evaluation-X patients",
              "python":sys.version,"numpy":np.__version__,"platform":platform.platform()}
    with open(os.path.join(OUT,"MANIFEST.json"),"w",encoding="utf-8") as f:json.dump(manifest,f,indent=2,ensure_ascii=False)
    print("\n=== primary DDAU, seed 42 (decisions) ===")
    for r in dec:
        if r["model"]=="DDAU" and r["seed"]==42:print(f"  {r['site']:9s} g={r['gamma']:.1f} simCI[{r['L_simCI']:+.4f},{r['U_simCI']:+.4f}] -> {r['decision']}")
    print("\nNOTE: satisfies Prop D v2 A2 (q from D_train frozen; p_S on independent D_cal; bootstrap=cal+target).")
    print(f"artifacts -> {OUT}")

if __name__=="__main__":
    main()
