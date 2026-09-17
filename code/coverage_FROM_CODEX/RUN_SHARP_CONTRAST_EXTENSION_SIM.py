"""
Sharp DIRECT Brier-risk contrast — extension simulation (publication-grade, regime-aware).
Evaluates v6 theory (연구정리_SiM_방법론코어_v6_sharp_contrast.md) across 5 pre-registered gamma* regimes:
  candidate_superior / reference_superior / ambiguous / finite_crossing / no_crossing_within_range.
Per regime: direct vs separate width (Prop A), width formula (Prop B), gamma-nestedness (Cor B1),
robustness value gamma* (Prop C) with CI & coverage (off-by-one-safe), gamma-SIMULTANEOUS band via
source+target joint nonparametric PAIRS bootstrap (sup-t), retain/reference/defer + regret (explicit defer_cost).
Firewall: target survival outcomes used ONLY for the simulation TRUTH; never for p_S/endpoint/gamma/decision.
Config locked (hash); effective (smoke) config recorded separately.

Usage:  python RUN_SHARP_CONTRAST_EXTENSION_SIM.py                                  # full spec
        FB_REPS=8 FB_BOOT=25 python RUN_...SIM.py                                    # smoke
        FB_ONLY_SEED=42 FB_ONLY_CENS=marginal FB_ONLY_REGIME=finite_crossing python  # shard
        FB_BENCH=1 python RUN_...SIM.py                                              # per-rep timing only
"""
import numpy as np, json, os, hashlib, time, math
sig=lambda z:1/(1+np.exp(-z))
def logit(p):p=np.clip(p,1e-6,1-1e-6);return np.log(p/(1-p))
def wilson(k,n,z=1.96):
    if n==0:return(0.0,1.0)
    p=k/n;d=1+z*z/n;c=(p+z*z/(2*n))/d;h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d;return max(0,c-h),min(1,c+h)
def norm_ppf(p):
    a=[-3.969683028665376e+01,2.209460984245205e+02,-2.759285104469687e+02,1.383577518672690e+02,-3.066479806614716e+01,2.506628277459239e+00]
    b=[-5.447609879822406e+01,1.615858368580409e+02,-1.556989798598866e+02,6.680131188771972e+01,-1.328068155288572e+01]
    c=[-7.784894002430293e-03,-3.223964580411365e-01,-2.400758277161838e+00,-2.549732539343734e+00,4.374664141464968e+00,2.938163982698783e+00]
    d=[7.784695709041462e-03,3.224671290700398e-01,2.445134137142996e+00,3.754408661907416e+00];pl=0.02425
    if p<pl:q=math.sqrt(-2*math.log(p));return(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p<=1-pl:q=p-0.5;r=q*q;return(((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q/(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
    q=math.sqrt(-2*math.log(1-p));return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)

HERE=os.path.dirname(os.path.abspath(__file__))
CFG=json.load(open(os.path.join(HERE,"SHARP_CONTRAST_EXTENSION_CONFIG.json")))
SEEDS=CFG["seeds"]; REPS=int(os.environ.get("FB_REPS",CFG["reps_total"])); BOOT=int(os.environ.get("FB_BOOT",CFG["boot"]))
GAMMAS=CFG["gamma_grid"]; tau=CFG["tau"]; k=CFG["weibull_k"]; GT=CFG["gamma_true_shift"]
nS=CFG["n_source"]; nT=CFG["n_target"]; CENS=CFG["censoring_modes"]
ALPHA=CFG["decision_alpha"]; ZCI=norm_ppf(1-ALPHA/2); SIMQ=1-ALPHA; DEFER_COST=CFG["defer_cost"]; GSTEP=CFG["gamma_star_grid"]["step"]
PS_I=CFG["pS_intercept"]; PS_S=CFG["pS_slope"]
REGIMES={k2:v for k2,v in CFG["dgp_regimes"].items() if k2!="note"}
if REPS%len(SEEDS): raise ValueError("reps_total must divide by #seeds")
CFG_HASH=hashlib.sha256(json.dumps(CFG,sort_keys=True).encode()).hexdigest()[:12]
CODE_SHA=hashlib.sha256(open(os.path.abspath(__file__),'rb').read()).hexdigest()
ACT_SEEDS=[int(os.environ["FB_ONLY_SEED"])] if "FB_ONLY_SEED" in os.environ else SEEDS
ACT_CENS=[os.environ["FB_ONLY_CENS"]] if "FB_ONLY_CENS" in os.environ else CENS
ACT_REG=[os.environ["FB_ONLY_REGIME"]] if "FB_ONLY_REGIME" in os.environ else list(REGIMES)
BENCH=os.environ.get("FB_BENCH")=="1"
IS_SHARD=any(e in os.environ for e in ["FB_ONLY_SEED","FB_ONLY_CENS","FB_ONLY_REGIME"])
AGG_ONLY=os.environ.get("FB_AGG_ONLY")=="1"   # read checkpoints only; error on missing; produce canonical AGG+manifest
OUT=os.path.join(HERE,"sharp_contrast_out",f"run_{CFG_HASH}"); os.makedirs(OUT,exist_ok=True)
_manifest={"requested_reps":CFG["reps_total"],"effective_reps":REPS,"requested_boot":CFG["boot"],"effective_boot":BOOT,
           "environment_overrides":{e:os.environ.get(e) for e in["FB_REPS","FB_BOOT","FB_ONLY_SEED","FB_ONLY_CENS","FB_ONLY_REGIME","FB_BENCH","FB_AGG_ONLY"] if e in os.environ},
           "is_smoke":(REPS!=CFG["reps_total"] or BOOT!=CFG["boot"]),"is_shard":IS_SHARD,"config_sha256":CFG_HASH,"code_sha256":CODE_SHA,
           "decision_alpha":ALPHA,"z_ci":ZCI,"sim_quantile":SIMQ,"defer_cost":DEFER_COST,"regimes":list(REGIMES)}
# shards write a tagged manifest (avoid concurrent overwrite of canonical); unfiltered/agg run writes canonical
if IS_SHARD and not AGG_ONLY:
    _tag="_".join(filter(None,[os.environ.get("FB_ONLY_REGIME"),os.environ.get("FB_ONLY_CENS"),os.environ.get("FB_ONLY_SEED")]))
    json.dump(_manifest,open(os.path.join(OUT,f"CONFIG_EXECUTED.shard_{_tag}.json"),"w"),indent=2)
else:
    json.dump(_manifest,open(os.path.join(OUT,"CONFIG_EXECUTED.json"),"w"),indent=2)
XED=np.array([-np.inf,-0.8416212336,-0.2533471031,0.2533471031,0.8416212336,np.inf])
xstr=lambda x:np.clip(np.searchsorted(XED,x,"right")-1,0,len(XED)-2)
pS=lambda x:sig(PS_I+PS_S*x)
def pT(x,rng):d=np.clip(0.15+0.6*x,-GT,GT);return sig(logit(pS(x))+d)

# ---- regime state (set per regime) ----
RG={"i0":0.8,"s0":0.7,"i1":0.4,"s1":0.7,"gmax":2.0}
def q0(x):return np.clip(sig(RG["i0"]+RG["s0"]*x),1e-3,1-1e-3)
def q1(x):return np.clip(sig(RG["i1"]+RG["s1"]*x),1e-3,1-1e-3)
def a(x):return q1(x)**2-q0(x)**2
def b(x):return q1(x)-q0(x)
def make_GG():return np.round(np.arange(0.0,RG["gmax"]+1e-9,GSTEP),6)

def gen(seed,n,tgt,cens):
    r=np.random.default_rng(seed);X=r.normal(0,1,n)
    p=pT(X,r) if tgt else pS(X);lam=-np.log(np.clip(p,1e-6,1-1e-6))/tau**k
    T=(-np.log(r.random(n))/np.clip(lam,1e-6,None))**(1/k)
    C=r.exponential(2.2,n) if cens=="marginal" else r.exponential(2.2*np.exp(0.35*(xstr(X)-2)),n)
    return X,np.minimum(T,C),(T<=C).astype(int)
def km(t,ci):
    o=np.argsort(t,kind="mergesort");t=t[o];c=ci[o];n=len(t);uq,f=np.unique(t,return_index=True);return uq,np.cumprod(1-np.add.reduceat(c,f)/(n-f))
def Gm(Tk,V,q):
    i=np.searchsorted(Tk,q,"left")-1;q=np.atleast_1d(q);o=np.ones(len(q));ok=i>=0;o[ok]=V[np.clip(i[ok],0,len(V)-1)];return o
def wl(x,y,w,it=15):
    Xd=np.c_[np.ones_like(x),x];bb=np.zeros(2)
    for _ in range(it):
        p=sig(Xd@bb);W=np.clip(w*p*(1-p),1e-12,None);bb=bb+np.linalg.solve((Xd*W[:,None]).T@Xd+1e-8*np.eye(2),Xd.T@(w*(y-p)))
    return bb
def phat(Xs,Tt,dl,cens,Q=5):
    censI=1-dl
    if cens=="marginal":
        Tk,V=km(Tt,censI);gi=np.clip(Gm(Tk,V,Tt),1e-3,None);gtau=np.full(len(Tt),max(Gm(Tk,V,np.array([tau]))[0],1e-3))
    else:
        s=xstr(Xs);kms={j:km(Tt[s==j],censI[s==j]) for j in range(Q) if (s==j).sum()>2};gi=np.ones(len(Tt));gtau=np.ones(len(Tt))
        for j,(Tk,V) in kms.items():
            m=s==j;gi[m]=Gm(Tk,V,Tt[m]);gtau[m]=Gm(Tk,V,np.array([tau]))[0]
        gi=np.clip(gi,1e-3,None);gtau=np.clip(gtau,1e-3,None)
    D=(Tt>tau).astype(float);w=np.zeros(len(Tt));w[Tt>tau]=1/gtau[Tt>tau];dd=(Tt<=tau)&(dl==1);w[dd]=1/gi[dd]
    h=len(Xs)//2;b1=wl(Xs[:h],D[:h],w[:h]);b2=wl(Xs[h:],D[h:],w[h:])
    return lambda x:0.5*(sig(b1[0]+b1[1]*x)+sig(b2[0]+b2[1]*x))
def LU_over_grid(Xt,pv,grid):
    A=a(Xt);B=b(Xt);lz=logit(pv);pos=B>=0;L=np.empty(len(grid));U=np.empty(len(grid))
    for i,g in enumerate(grid):
        pm=sig(lz-g);pp=sig(lz+g);pL=np.where(pos,pp,pm);pU=np.where(pos,pm,pp);L[i]=np.mean(A-2*B*pL);U[i]=np.mean(A-2*B*pU)
    return L,U
def sep_width(Xt,pv,g):
    lz=logit(pv);pm=sig(lz-g);pp=sig(lz+g)
    def w(q):coef=1-2*q(Xt);pmax=np.where(coef>=0,pp,pm);pmin=np.where(coef>=0,pm,pp);return np.mean(coef*(pmax-pmin))
    return w(q1)+w(q0)
# gamma* boundaries (off-by-one-safe, reviewer's def)
def gstar_point(U,grid):m=grid[U<0];return float(m.max()) if m.size else 0.0          # sup{g:U<0}
def gstar_lo(U,grid):  m=grid[U<0];return float(m.max()) if m.size else 0.0           # max{g:Ubar<0}
def gstar_hi(U,grid):                                                                 # min{g:Ulow>=0}; censored if none
    m=grid[U>=0]
    return (float(m.min()),0) if m.size else (float(grid.max()),1)

def truth(GG):
    r=np.random.default_rng(1);X=r.normal(0,1,400_000);Lt,Ut=LU_over_grid(X,pS(X),GG)
    dRT=float(np.mean(a(X)-2*b(X)*pT(X,np.random.default_rng(2))))
    gt=gstar_point(Ut,GG);gt_cens=int((Ut<0).all())
    return Lt,Ut,dRT,gt,gt_cens

def one_rep(seed,cens,GG,TRUE_L,TRUE_U,dRT_true):
    Xs,Ts,ds=gen(seed,nS,False,cens);Xt,_,_=gen(seed+7,nT,True,cens)
    ph=phat(Xs,Ts,ds,cens);pv=ph(Xt);Lhat,Uhat=LU_over_grid(Xt,pv,GG)
    rng=np.random.default_rng(seed*131+(0 if cens=="marginal" else 1))
    Lb=np.empty((BOOT,len(GG)));Ub=np.empty((BOOT,len(GG)))
    for jb in range(BOOT):
        si=rng.integers(0,nS,nS);ti=rng.integers(0,nT,nT);phb=phat(Xs[si],Ts[si],ds[si],cens);pvb=phb(Xt[ti])
        Lb[jb],Ub[jb]=LU_over_grid(Xt[ti],pvb,GG)
    seL=Lb.std(0)+1e-12;seU=Ub.std(0)+1e-12
    W=np.maximum(np.abs(Lb-Lhat)/seL,np.abs(Ub-Uhat)/seU).max(1);csim=np.quantile(W,SIMQ)
    simu_cov=int(np.all((Lhat-csim*seL<=TRUE_L)&(TRUE_L<=Lhat+csim*seL)&(Uhat-csim*seU<=TRUE_U)&(TRUE_U<=Uhat+csim*seU)))
    gstar_hat=gstar_point(Uhat,GG);g_lo=gstar_lo(Uhat+csim*seU,GG);g_hi,hi_cens=gstar_hi(Uhat-csim*seU,GG)
    GS_true=gstar_point(TRUE_U,GG)
    gstar_cov=int(g_lo-1e-9<=GS_true<=g_hi+1e-9)
    GP=[g for g in GAMMAS if g<=GG.max()+1e-9];gp=0.3 if 0.3 in GP else GP[-1];i3=int(np.argmin(np.abs(GG-gp)))
    L3,U3,sL3,sU3=Lhat[i3],Uhat[i3],seL[i3],seU[i3];wdir3=U3-L3;wsep3=sep_width(Xt,pv,gp)
    point3=int((L3-ZCI*sL3<=TRUE_L[i3]<=L3+ZCI*sL3) and (U3-ZCI*sU3<=TRUE_U[i3]<=U3+ZCI*sU3))
    # Use the simultaneous critical value from the prespecified gamma-grid
    # envelope. A pointwise normal critical value can over-certify decisions.
    dec="retain" if U3+csim*sU3<0 else ("reference" if L3-csim*sL3>0 else "defer")
    reg=(max(0.0,dRT_true) if dec=="retain" else (max(0.0,-dRT_true) if dec=="reference" else DEFER_COST))
    return {"seed":seed,"cens":cens,"simu_cov":simu_cov,"point_cov_03":point3,"width_dir_03":wdir3,"width_sep_03":wsep3,
            "strict_tighter_03":int(wsep3>wdir3+1e-9),"gstar_hat":gstar_hat,"gstar_lo":g_lo,"gstar_hi":g_hi,
            "gstar_hi_censored":hi_cens,"gstar_cov":gstar_cov,"decision":dec,"regret":reg}

def valid_ckpt(p,reg,seed,cens,n):
    if not os.path.isfile(p):return None
    try:d=json.load(open(p))
    except Exception:return None
    m=d.get("meta",{});return d if (m.get("cfg")==CFG_HASH and m.get("reg")==reg and m.get("seed")==seed and m.get("cens")==cens and m.get("boot")==BOOT and m.get("n")==n) else None

def main():
    per=REPS//len(SEEDS);print(f"[cfg {CFG_HASH}] reps/seed={per} boot={BOOT} z={ZCI:.3f} simQ={SIMQ} regimes={ACT_REG}")
    if BENCH:
        for reg in ACT_REG:
            RG.update(REGIMES[reg]);GG=make_GG();TL,TU,dR,gt,gc=truth(GG)
            t=time.time();_=one_rep(999,"marginal",GG,TL,TU,dR);dt=time.time()-t
            full=dt*REPS*len(CENS);print(f"  [bench] {reg}: {dt:.2f}s/rep (boot={BOOT},|GG|={len(GG)}) -> est {full/3600:.2f}h for {REPS}x{len(CENS)}cens (this regime)")
        return
    for reg in ACT_REG:
        RG.update(REGIMES[reg]);GG=make_GG();TL,TU,dRT,GS_true,GS_cens=truth(GG)
        ORACLE="retain" if dRT<0 else ("reference" if dRT>0 else "defer")
        for cens in ACT_CENS:
            recs=[]
            for seed in ACT_SEEDS:
                ck=os.path.join(OUT,f"ckpt_{reg}_{cens}_{seed}.json");v=valid_ckpt(ck,reg,seed,cens,per)
                if v:recs+=v["recs"];print(f"  skip {reg}/{cens}/seed{seed}");continue
                if AGG_ONLY:raise FileNotFoundError(f"[AGG_ONLY] missing/invalid checkpoint {ck}; run its shard first")
                sr=[one_rep(seed*1000+i,cens,GG,TL,TU,dRT) for i in range(per)]
                json.dump({"meta":{"cfg":CFG_HASH,"reg":reg,"seed":seed,"cens":cens,"boot":BOOT,"n":per},"recs":sr},open(ck,"w"))
                recs+=sr;print(f"  done {reg}/{cens}/seed{seed} ({per})")
            if len(ACT_SEEDS)<len(SEEDS) or len(ACT_REG)<len(REGIMES):continue
            def rate(key):kk=sum(r[key] for r in recs);n=len(recs);lo,hi=wilson(kk,n,ZCI);return[round(kk/n,4),round(lo,4),round(hi,4)]
            decs=[r["decision"] for r in recs]
            agg={"regime":reg,"cens":cens,"n":len(recs),"oracle":ORACLE,"dRT_true":round(dRT,4),
                 "gamma*_true":round(GS_true,3),"gamma*_true_censored":GS_cens,
                 "simultaneous_cov":rate("simu_cov"),"pointwise_cov_g0.3":rate("point_cov_03"),
                 "strict_tighter_rate":rate("strict_tighter_03"),
                 "mean_width_direct":round(float(np.mean([r["width_dir_03"] for r in recs])),5),
                 "mean_width_separate":round(float(np.mean([r["width_sep_03"] for r in recs])),5),
                 "gamma*_hat_mean":round(float(np.mean([r["gstar_hat"] for r in recs])),3),
                 "gamma*_hat_rmse":round(float(np.sqrt(np.mean([(r["gstar_hat"]-GS_true)**2 for r in recs]))),3),
                 "gamma*_CI_cov":rate("gstar_cov"),"gamma*_hi_censored_rate":round(float(np.mean([r["gstar_hi_censored"] for r in recs])),3),
                 "decision_rates":{d:round(decs.count(d)/len(decs),3) for d in["retain","reference","defer"]},
                 "mean_regret":round(float(np.mean([r["regret"] for r in recs])),5)}
            agg["width_ratio_sep_over_dir"]=round(agg["mean_width_separate"]/max(agg["mean_width_direct"],1e-9),3)
            json.dump(agg,open(os.path.join(OUT,f"AGG_{reg}_{cens}.json"),"w"),indent=2);print(json.dumps(agg,indent=2))

if __name__=="__main__":
    t=time.time();main();print(f"[elapsed {time.time()-t:.1f}s] out={OUT}")
