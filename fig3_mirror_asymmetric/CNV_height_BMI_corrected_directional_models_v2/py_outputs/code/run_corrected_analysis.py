from __future__ import annotations
import argparse, json, math, os, shutil, textwrap, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize, minimize_scalar
from scipy.special import expit, logsumexp
from scipy.stats import spearmanr, mannwhitneyu, ttest_ind
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_CELL_VERTICAL_ALIGNMENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

SEED=20260726
rng=np.random.default_rng(SEED)

# --- Paths -------------------------------------------------------------------
# Originally hardcoded to the authoring sandbox (/mnt/data/...), which made this
# script unrunnable outside that environment. Inputs are now read from the
# package itself and outputs written to a separate directory by default.
_ap=argparse.ArgumentParser(description='Corrected directional models of reciprocal CNV effects on height and BMI.')
_ap.add_argument('--project-dir',default=str(Path(__file__).resolve().parent.parent),
                 help='Package root containing height/inputs and bmi/inputs (default: parent of code/).')
_ap.add_argument('--out-dir',default=None,
                 help='Where to write outputs (default: <project-dir>/py_outputs).')
_ap.add_argument('--in-place',action='store_true',
                 help='Write back into <project-dir>, overwriting shipped outputs, figures, report, README.md and SOURCES.md.')
_ap.add_argument('--height-src',default=None,help='Optional external directory to copy height inputs from.')
_ap.add_argument('--bmi-src',default=None,help='Optional external directory to copy BMI inputs from.')
_args=_ap.parse_args()

PROJECT=Path(_args.project_dir).resolve()
if not (PROJECT/'height'/'inputs').is_dir():
    raise SystemExit(f'--project-dir does not look like the package root (no height/inputs): {PROJECT}')
ROOT=PROJECT if _args.in_place else Path(_args.out_dir or PROJECT/'py_outputs').resolve()

# Never delete the package itself; only a dedicated output directory is cleared.
def _rmtree_resilient(path,attempts=6,pause=1.0):
    """Windows + cloud-sync clients (Dropbox, OneDrive) transiently lock files,
    which makes a plain rmtree fail with WinError 32. Retry, then give up
    gracefully rather than aborting the whole run."""
    import time
    for i in range(attempts):
        try:
            shutil.rmtree(path);return True
        except (PermissionError,OSError) as e:
            if i==attempts-1:
                print(f'[warn] could not fully clear {path} ({e}); reusing directory and overwriting in place',flush=True)
                shutil.rmtree(path,ignore_errors=True);return False
            time.sleep(pause)
if ROOT!=PROJECT and ROOT.exists(): _rmtree_resilient(ROOT)
for sub in ['height/inputs','height/outputs','height/figures','height/code','bmi/inputs','bmi/outputs','bmi/figures','bmi/code','report','code','qa']:
    (ROOT/sub).mkdir(parents=True,exist_ok=True)

# Inputs come from the package unless an external source is supplied.
H_SRC=Path(_args.height_src).resolve() if _args.height_src else PROJECT/'height'
B_SRC=Path(_args.bmi_src).resolve() if _args.bmi_src else PROJECT/'bmi'
def _stage(src_dir,names,dest_sub):
    missing=[]
    for fn in names:
        src=src_dir/'inputs'/fn
        if not src.exists(): missing.append(fn); continue
        dst=ROOT/dest_sub/fn
        if src.resolve()!=dst.resolve(): shutil.copy2(src,dst)
    if missing: raise SystemExit(f'Missing required inputs in {src_dir/"inputs"}: {", ".join(missing)}')
_stage(H_SRC,['sacks_height_all_loci_exact.csv','sacks_height_mirror_loci_exact.csv','kosmicki_height_direction_reference.csv','Sacks_supplementary_tables.xlsx'],'height/inputs')
_stage(B_SRC,['sacks_bmi_complete_reciprocal_loci_exact.csv','sacks_bmi_reciprocal_loci_audit.csv','akbari_2021_exomewide_BMI_genes.csv','akbari_2021_pLoF_direction_reference.csv','backman_2021_usage_note.csv','locus_gene_membership.csv'],'bmi/inputs')
_stage(B_SRC,['Sacks_supplementary_tables.xlsx'],'bmi/inputs')
print(f'[paths] project={PROJECT}\n[paths] output ={ROOT}',flush=True)

height=pd.read_csv(ROOT/'height/inputs/sacks_height_all_loci_exact.csv')
height_mirror=height[height.Dose_response_class.str.lower().eq('mirror')].copy().reset_index(drop=True)
bmi=pd.read_csv(ROOT/'bmi/inputs/sacks_bmi_complete_reciprocal_loci_exact.csv')

# Literature targets
burden_targets=pd.DataFrame([
    ['Mace_2017','height','deletion',-0.41,'cm_per_Mb','Published deletion-burden coefficient'],
    ['Mace_2017','BMI','deletion',0.14,'kg_m2_per_Mb','Published deletion-burden coefficient'],
    ['Auwerx_2022','height','all_CNV',np.nan,'direction_only','Greater CNV burden associated with lower physical traits including height'],
    ['Auwerx_2022','BMI/adiposity','all_CNV',np.nan,'direction_only','Greater CNV burden associated with increased adiposity'],
],columns=['source','trait','copy_state','estimate','units','interpretation'])
burden_targets.to_csv(ROOT/'code/literature_burden_targets.csv',index=False)

# ---------------- model functions ----------------
def bvlogpdf(x,y,m1,m2,v1,v2,cov):
    x=np.asarray(x);y=np.asarray(y);m1=np.asarray(m1);m2=np.asarray(m2)
    v1=np.asarray(v1);v2=np.asarray(v2);cov=np.asarray(cov)
    det=v1*v2-cov*cov
    if np.any(det<=1e-14) or np.any(~np.isfinite(det)): return None
    dx=x-m1;dy=y-m2
    quad=(v2*dx*dx+v1*dy*dy-2*cov*dx*dy)/det
    return -np.log(2*np.pi)-.5*np.log(det)-.5*quad

def mirror_components(dat,z):
    D=dat.Deletion_beta.to_numpy(float);U=dat.Duplication_beta.to_numpy(float)
    G=dat.Gene_count.to_numpy(float); seD=dat.Deletion_SE.to_numpy(float);seU=dat.Duplication_SE.to_numpy(float)
    p,m,rho,tD,tU,mus=z['p'],z['m'],z['rho'],z['tauD'],z['tauU'],z['sign_mean']
    mean_gene=p*mus*m
    var_gene=p*2*m*m-(p*mus*m)**2
    var_gene=max(var_gene,1e-12)
    v=G*var_gene
    md=G*mean_gene;mu=-rho*G*mean_gene
    vd=v+seD**2+tD**2;vu=rho*rho*v+seU**2+tU**2;cv=-rho*v
    return D,U,G,md,mu,vd,vu,cv

def fit_mirror(dat,pi_positive,nstarts=10):
    D=dat.Deletion_beta.to_numpy(float);U=dat.Duplication_beta.to_numpy(float)
    G=dat.Gene_count.to_numpy(float);seD=dat.Deletion_SE.to_numpy(float);seU=dat.Duplication_SE.to_numpy(float)
    mus=2*pi_positive-1
    def unpack(x):
        return dict(p=float(expit(x[0])),m=float(np.exp(x[1])),rho=float(np.exp(x[2])),tauD=float(np.exp(x[3])),tauU=float(np.exp(x[4])),sign_mean=mus,pi_positive=pi_positive)
    def nll(x):
        z=unpack(x);p,m,rho,tD,tU=z['p'],z['m'],z['rho'],z['tauD'],z['tauU']
        mg=p*mus*m;vg=p*2*m*m-(p*mus*m)**2
        if vg<=1e-12 or not np.isfinite(vg): return 1e100
        v=G*vg;md=G*mg;mu=-rho*G*mg;vd=v+seD**2+tD**2;vu=rho*rho*v+seU**2+tU**2;cv=-rho*v
        ll=bvlogpdf(D,U,md,mu,vd,vu,cv)
        return 1e100 if ll is None else float(-ll.sum())
    p0=.4 if pi_positive<.5 else .03
    m0=max(.02,abs(D.mean())/(max(abs(mus),.1)*G.mean()*p0))
    slope=np.polyfit(D,U,1)[0]
    rho0=max(.15,min(1.0,abs(slope)))
    x0=np.array([math.log(p0/(1-p0)),math.log(m0),math.log(rho0),math.log(.02),math.log(.05)])
    scales=np.array([1.4,1,.6,1,1]);best=None
    for i in range(nstarts):
        rg=np.random.default_rng(SEED+107*i+int(pi_positive*1000))
        r=minimize(nll,x0+rg.normal(0,scales),method='Nelder-Mead',options={'maxiter':8000,'xatol':1e-8,'fatol':1e-8})
        if best is None or r.fun<best.fun:best=r
    z=unpack(best.x);k=5
    z.update(model='Mirror alone',negative_log_likelihood=float(best.fun),logLik=float(-best.fun),n_parameters=k,AIC=float(2*k+2*best.fun),BIC=float(k*np.log(len(dat))+2*best.fun))
    return z,best.x

def frequency_effect_profile(dat,pi_positive,p_grid=None,nstarts=4):
    """Frequency-effect ridge: profile the mirror likelihood over contributor
    frequency p, optimising (m, rho, tauD, tauU) at each fixed p and recording the
    best mean absolute gene effect m.

    Reconstructed in-package. The original curve was read from
    {height,BMI}_frequency_effect_profile.csv in the upstream _rework_* packages,
    which were never shipped with this bundle, so the figure could not be
    regenerated. Recomputing from the same likelihood keeps the ridge exactly
    consistent with the mirror fit reported here.
    """
    D=dat.Deletion_beta.to_numpy(float);U=dat.Duplication_beta.to_numpy(float)
    G=dat.Gene_count.to_numpy(float);seD=dat.Deletion_SE.to_numpy(float);seU=dat.Duplication_SE.to_numpy(float)
    mus=2*pi_positive-1
    if p_grid is None: p_grid=np.linspace(.005,.95,60)
    def nll_given_p(p,x):
        # Nelder-Mead is unconstrained and can push the log-parameters far enough
        # that exp() overflows; clamp before exponentiating.
        xs=np.clip(np.asarray(x,dtype=float),-50.0,50.0)
        if not np.all(np.isfinite(xs)): return 1e100
        m,rho,tD,tU=(float(np.exp(v)) for v in xs[:4])
        mg=p*mus*m;vg=p*2*m*m-(p*mus*m)**2
        if vg<=1e-12 or not np.isfinite(vg): return 1e100
        v=G*vg;md=G*mg;mu=-rho*G*mg;vd=v+seD**2+tD**2;vu=rho*rho*v+seU**2+tU**2;cv=-rho*v
        if not (np.all(np.isfinite(vd)) and np.all(np.isfinite(vu)) and np.all(np.isfinite(cv))): return 1e100
        ll=bvlogpdf(D,U,md,mu,vd,vu,cv)
        if ll is None: return 1e100
        s=float(-ll.sum())
        return s if np.isfinite(s) else 1e100
    rows=[];x0=np.array([math.log(.05),math.log(.5),math.log(.02),math.log(.05)]);scales=np.array([1,.6,1,1])
    for i,p in enumerate(p_grid):
        best=None
        for s in range(nstarts):
            rg=np.random.default_rng(SEED+911*i+s)
            r=minimize(lambda x:nll_given_p(float(p),x),x0+rg.normal(0,scales),method='Nelder-Mead',
                       options={'maxiter':6000,'xatol':1e-8,'fatol':1e-8})
            if best is None or r.fun<best.fun: best=r
        rows.append({'mirror_gene_frequency':float(p),'mean_abs_gene_effect_SD':float(np.exp(np.clip(best.x[0],-50,50))),
                     'negative_log_likelihood':float(best.fun)})
        x0=np.clip(best.x,-50,50)
    return pd.DataFrame(rows)

def load_frequency_effect_profile(src_dir,fname,dat,pi_positive,out_path):
    """Prefer a shipped upstream profile; otherwise reconstruct and save it."""
    f=src_dir/'outputs'/fname
    if f.exists():
        print(f'[profile] using shipped {f}',flush=True)
        return pd.read_csv(f)
    print(f'[profile] {fname} not present; reconstructing ridge -> {out_path}',flush=True)
    prof=frequency_effect_profile(dat,pi_positive)
    prof.to_csv(out_path,index=False)
    return prof

def fit_directional(dat,base,trait_sign,model,max_q,nstarts=14,angle_min=5,angle_max=85):
    D,U,G,md,mu,vd,vu,cv=mirror_components(dat,base)
    if model=='mirror':
        ll=bvlogpdf(D,U,md,mu,vd,vu,cv)
        return dict(model='Mirror alone',q=0.,shift_magnitude=0.,shift_angle_deg=45.,shift_DEL=0.,shift_DUP=0.,negative_log_likelihood=float(-ll.sum()),n_parameters_added=0)
    def unpack(x):
        q=max_q*expit(x[0]);b=np.exp(x[1]);theta=np.deg2rad(angle_min+(angle_max-angle_min)*expit(x[2]))
        return float(q),float(b),float(theta)
    def nll(x):
        q,b,theta=unpack(x);bD=b*np.cos(theta);bU=b*np.sin(theta)
        if model=='continuous':
            vv=G*(2*q-q*q)
            ll=bvlogpdf(D,U,md+trait_sign*G*q*bD,mu+trait_sign*G*q*bU,
                        vd+vv*bD*bD,vu+vv*bU*bU,cv+vv*bD*bU)
            return 1e100 if ll is None else float(-ll.sum())
        p0=(1-q)**G
        ll0=np.log(np.maximum(p0,1e-300))+bvlogpdf(D,U,md,mu,vd,vu,cv)
        ll1=np.log(np.maximum(1-p0,1e-300))+bvlogpdf(D,U,md+trait_sign*bD,mu+trait_sign*bU,vd,vu,cv)
        return float(-logsumexp(np.column_stack([ll0,ll1]),axis=1).sum())
    q0=.08 if model=='continuous' else .018
    b0=.055 if model=='continuous' else .25
    # start with mostly vertical shift because the fitted data often favor that, but explore broadly
    frac=(70-angle_min)/(angle_max-angle_min)
    x0=np.array([math.log(q0/(max_q-q0)),math.log(b0),math.log(frac/(1-frac))])
    best=None
    for i in range(nstarts):
        rg=np.random.default_rng(SEED+211*i+(1 if model=='continuous' else 2))
        r=minimize(nll,x0+rg.normal(0,[1.5,1.2,1.6]),method='Nelder-Mead',options={'maxiter':8000,'xatol':1e-8,'fatol':1e-8})
        if best is None or r.fun<best.fun:best=r
    q,b,theta=unpack(best.x);bD=b*np.cos(theta);bU=b*np.sin(theta);k=3
    label='Mirror + continuous directional shift' if model=='continuous' else 'Mirror + sparse single hit'
    return dict(model=label,q=q,shift_magnitude=b,shift_angle_deg=float(np.rad2deg(theta)),shift_DEL=float(bD),shift_DUP=float(bU),negative_log_likelihood=float(best.fun),n_parameters_added=k,
                conditional_AIC=float(2*k+2*best.fun),conditional_BIC=float(k*np.log(len(dat))+2*best.fun),max_q=max_q,raw_par=best.x.tolist())

def logpred_row(row,base,trait_sign,model,pars):
    one=pd.DataFrame([row]);D,U,G,md,mu,vd,vu,cv=mirror_components(one,base)
    D=float(D[0]);U=float(U[0]);G=float(G[0]);md=float(md[0]);mu=float(mu[0]);vd=float(vd[0]);vu=float(vu[0]);cv=float(cv[0])
    if model=='mirror': return float(bvlogpdf([D],[U],[md],[mu],[vd],[vu],[cv])[0])
    q=pars['q'];bD=pars['shift_DEL'];bU=pars['shift_DUP']
    if model=='continuous':
        vv=G*(2*q-q*q)
        return float(bvlogpdf([D],[U],[md+trait_sign*G*q*bD],[mu+trait_sign*G*q*bU],[vd+vv*bD*bD],[vu+vv*bU*bU],[cv+vv*bD*bU])[0])
    p0=(1-q)**G
    l0=math.log(max(p0,1e-300))+float(bvlogpdf([D],[U],[md],[mu],[vd],[vu],[cv])[0])
    l1=math.log(max(1-p0,1e-300))+float(bvlogpdf([D],[U],[md+trait_sign*bD],[mu+trait_sign*bU],[vd],[vu],[cv])[0])
    return float(logsumexp([l0,l1]))

def fit_directional_fast(dat,base,trait_sign,model,max_q,xstart):
    D,U,G,md,mu,vd,vu,cv=mirror_components(dat,base)
    def unpack(x):
        q=max_q*expit(x[0]);b=np.exp(x[1]);theta=np.deg2rad(5+80*expit(x[2]));return float(q),float(b),float(theta)
    def nll(x):
        q,b,theta=unpack(x);bD=b*np.cos(theta);bU=b*np.sin(theta)
        if model=='continuous':
            vv=G*(2*q-q*q);ll=bvlogpdf(D,U,md+trait_sign*G*q*bD,mu+trait_sign*G*q*bU,vd+vv*bD*bD,vu+vv*bU*bU,cv+vv*bD*bU);return float(-ll.sum())
        p0=(1-q)**G;ll0=np.log(np.maximum(p0,1e-300))+bvlogpdf(D,U,md,mu,vd,vu,cv);ll1=np.log(np.maximum(1-p0,1e-300))+bvlogpdf(D,U,md+trait_sign*bD,mu+trait_sign*bU,vd,vu,cv);return float(-logsumexp(np.c_[ll0,ll1],axis=1).sum())
    r=minimize(nll,np.asarray(xstart),method='L-BFGS-B',bounds=[(-12,12),(-9,2),(-8,8)],options={'maxiter':800})
    if not r.success:r=minimize(nll,np.asarray(xstart),method='Nelder-Mead',options={'maxiter':3000})
    q,b,theta=unpack(r.x)
    return dict(q=q,shift_magnitude=b,shift_angle_deg=float(np.rad2deg(theta)),shift_DEL=float(b*np.cos(theta)),shift_DUP=float(b*np.sin(theta)),raw_par=r.x.tolist())

def conditional_loo(dat,base,trait_sign,fullfits):
    rows=[]
    for mkey,label in [('mirror','Mirror alone'),('continuous','Mirror + continuous directional shift'),('single','Mirror + sparse single hit')]:
        vals=[]
        for i in range(len(dat)):
            if mkey=='mirror': pars=None
            else:
                maxq=.8 if mkey=='continuous' else .2
                pars=fit_directional_fast(dat.drop(dat.index[i]),base,trait_sign,mkey,maxq,fullfits[mkey]['raw_par'])
            vals.append(logpred_row(dat.iloc[i],base,trait_sign,mkey,pars))
            rows.append({'model':label,'held_out_locus':dat.iloc[i].get('CNV_ID',dat.iloc[i].get('CNV_ID_raw',str(i))),'held_out_log_predictive_density':vals[-1]})
    tbl=pd.DataFrame(rows)
    sums=tbl.groupby('model',as_index=False).held_out_log_predictive_density.sum().rename(columns={'held_out_log_predictive_density':'conditional_LOO_log_predictive_density'})
    return tbl,sums

def simulate_model(dat,base,trait_sign,model,pars,nsim=8000,seed=SEED):
    rg=np.random.default_rng(seed);D,U,G,md,mu,vd,vu,cv=mirror_components(dat,base)
    # remove sampling SE from base covariance then add observation covariance through direct bivariate sampling using full covariance
    n=len(dat);stats=np.zeros((nsim,9));cloud=[]
    for s in range(nsim):
        ds=np.zeros(n);us=np.zeros(n)
        for j in range(n):
            mean=np.array([md[j],mu[j]]);cov=np.array([[vd[j],cv[j]],[cv[j],vu[j]]])
            if model=='continuous':
                q=pars['q'];bD=pars['shift_DEL'];bU=pars['shift_DUP'];vv=G[j]*(2*q-q*q)
                mean=mean+trait_sign*G[j]*q*np.array([bD,bU]);cov=cov+vv*np.outer([bD,bU],[bD,bU])
            elif model=='single':
                if rg.random()<1-(1-pars['q'])**G[j]:mean=mean+trait_sign*np.array([pars['shift_DEL'],pars['shift_DUP']])
            ds[j],us[j]=rg.multivariate_normal(mean,cov)
        stats[s]=[ds.mean(),us.mean(),np.corrcoef(ds,us)[0,1],ds.min(),ds.max(),us.min(),us.max(),np.mean(ds<0),np.mean(us<0)]
        if s<120:cloud.append(pd.DataFrame({'simulation':s,'Deletion_beta':ds,'Duplication_beta':us,'model':model}))
    return stats,pd.concat(cloud,ignore_index=True)

def model_locus_expectations(dat,base,trait_sign,model,pars):
    D,U,G,md,mu,vd,vu,cv=mirror_components(dat,base)
    out=dat.copy()
    out['mirror_expected_DEL']=md;out['mirror_expected_DUP']=mu
    if model=='mirror':ed,eu=md,mu;phit=np.zeros(len(dat))
    elif model=='continuous':
        ed=md+trait_sign*G*pars['q']*pars['shift_DEL'];eu=mu+trait_sign*G*pars['q']*pars['shift_DUP'];phit=np.repeat(np.nan,len(dat))
    else:
        phit=1-(1-pars['q'])**G;ed=md+trait_sign*phit*pars['shift_DEL'];eu=mu+trait_sign*phit*pars['shift_DUP']
    out['model_expected_DEL']=ed;out['model_expected_DUP']=eu;out['single_hit_probability']=phit
    out['residual_DEL_from_mirror']=out.Deletion_beta-md;out['residual_DUP_from_mirror']=out.Duplication_beta-mu
    out['directional_shift_score']=trait_sign*(out.residual_DEL_from_mirror+out.residual_DUP_from_mirror)/np.sqrt(2)
    return out

def burden_values(base,trait_sign,model,pars,Gvals):
    p,m,rho,mus=base['p'],base['m'],base['rho'],base['sign_mean']
    d0=p*mus*m;u0=-rho*d0
    Gvals=np.asarray(Gvals,float)
    if model=='mirror':d=np.repeat(d0,len(Gvals));u=np.repeat(u0,len(Gvals))
    elif model=='continuous':
        d=np.repeat(d0+trait_sign*pars['q']*pars['shift_DEL'],len(Gvals));u=np.repeat(u0+trait_sign*pars['q']*pars['shift_DUP'],len(Gvals))
    else:
        occ=1-(1-pars['q'])**Gvals;d=d0+trait_sign*pars['shift_DEL']*occ/Gvals;u=u0+trait_sign*pars['shift_DUP']*occ/Gvals
    return d,u

def profile_q(dat,base,trait_sign,model,maxq,fullfit,grid):
    D,U,G,md,mu,vd,vu,cv=mirror_components(dat,base)
    rows=[]
    for q in grid:
        def unpack2(x):
            b=np.exp(x[0]);theta=np.deg2rad(5+80*expit(x[1]));return b,theta
        def nll2(x):
            b,theta=unpack2(x);bD=b*np.cos(theta);bU=b*np.sin(theta)
            if model=='continuous':
                vv=G*(2*q-q*q);ll=bvlogpdf(D,U,md+trait_sign*G*q*bD,mu+trait_sign*G*q*bU,vd+vv*bD*bD,vu+vv*bU*bU,cv+vv*bD*bU);return float(-ll.sum())
            p0=(1-q)**G;ll0=np.log(np.maximum(p0,1e-300))+bvlogpdf(D,U,md,mu,vd,vu,cv);ll1=np.log(np.maximum(1-p0,1e-300))+bvlogpdf(D,U,md+trait_sign*bD,mu+trait_sign*bU,vd,vu,cv);return float(-logsumexp(np.c_[ll0,ll1],axis=1).sum())
        frac=np.clip((fullfit['shift_angle_deg']-5)/80,1e-6,1-1e-6);x0=[math.log(fullfit['shift_magnitude']),math.log(frac/(1-frac))]
        r=minimize(nll2,x0,method='Nelder-Mead',options={'maxiter':1200})
        b,th=unpack2(r.x);rows.append({'q':q,'shift_magnitude':b,'shift_angle_deg':np.rad2deg(th),'negative_log_likelihood':r.fun})
    z=pd.DataFrame(rows);z['delta_2logLik']=2*(z.negative_log_likelihood-z.negative_log_likelihood.min());return z

# ---------------- fit models ----------------
height_base,_=fit_mirror(height_mirror,.30)
bmi_base,_=fit_mirror(bmi,11/14)

Hfits={
    'mirror':fit_directional(height,height_base,-1,'mirror',0),
    'continuous':fit_directional(height,height_base,-1,'continuous',.8),
    'single':fit_directional(height,height_base,-1,'single',.2),
}
Bfits={
    'mirror':fit_directional(bmi,bmi_base,1,'mirror',0),
    'continuous':fit_directional(bmi,bmi_base,1,'continuous',.8),
    'single':fit_directional(bmi,bmi_base,1,'single',.2),
}
for fits,n in [(Hfits,len(height)),(Bfits,len(bmi))]:
    fits['mirror']['conditional_AIC']=2*fits['mirror']['negative_log_likelihood'];fits['mirror']['conditional_BIC']=2*fits['mirror']['negative_log_likelihood']

Hloo,Hloosum=conditional_loo(height,height_base,-1,Hfits)
Bloo,Bloosum=conditional_loo(bmi,bmi_base,1,Bfits)

# merge comparisons
Hcomp=pd.DataFrame([{k:v for k,v in z.items() if k!='raw_par'} for z in Hfits.values()]).merge(Hloosum,on='model',how='left')
Bcomp=pd.DataFrame([{k:v for k,v in z.items() if k!='raw_par'} for z in Bfits.values()]).merge(Bloosum,on='model',how='left')
for comp in [Hcomp,Bcomp]:
    comp['delta_conditional_AIC']=comp.conditional_AIC-comp.conditional_AIC.min()
    comp['delta_conditional_BIC']=comp.conditional_BIC-comp.conditional_BIC.min()
    comp['delta_conditional_LOO']=comp.conditional_LOO_log_predictive_density.max()-comp.conditional_LOO_log_predictive_density

# Profiles
Hqcont=profile_q(height,height_base,-1,'continuous',.8,Hfits['continuous'],np.linspace(.005,.25,45))
Hqsingle=profile_q(height,height_base,-1,'single',.2,Hfits['single'],np.linspace(.001,.08,45))
Bqcont=profile_q(bmi,bmi_base,1,'continuous',.8,Bfits['continuous'],np.linspace(.01,.8,45))
Bqsingle=profile_q(bmi,bmi_base,1,'single',.2,Bfits['single'],np.linspace(.001,.15,45))

# Locus expectations
Hpred={k:model_locus_expectations(height,height_base,-1,k,Hfits[k]) for k in Hfits}
Bpred={k:model_locus_expectations(bmi,bmi_base,1,k,Bfits[k]) for k in Bfits}

# Simulations and PPC
stat_names=['DEL_mean','DUP_mean','DEL_DUP_correlation','DEL_min','DEL_max','DUP_min','DUP_max','DEL_negative_fraction','DUP_negative_fraction']
def empirical(dat):
 D=dat.Deletion_beta.to_numpy();U=dat.Duplication_beta.to_numpy();return np.array([D.mean(),U.mean(),np.corrcoef(D,U)[0,1],D.min(),D.max(),U.min(),U.max(),np.mean(D<0),np.mean(U<0)])
def ppc_for(dat,base,sign,fits,trait):
 obs=empirical(dat);rows=[];clouds=[]
 for i,key in enumerate(['mirror','continuous','single']):
  st,cl=simulate_model(dat,base,sign,key,fits[key],1200,SEED+500*i+(0 if trait=='height' else 50));clouds.append(cl)
  for j,nm in enumerate(stat_names):rows.append({'model':fits[key]['model'],'statistic':nm,'observed':obs[j],'simulated_mean':st[:,j].mean(),'simulated_2.5%':np.quantile(st[:,j],.025),'simulated_97.5%':np.quantile(st[:,j],.975),'posterior_predictive_CDF':np.mean(st[:,j]<=obs[j])})
 return pd.DataFrame(rows),pd.concat(clouds,ignore_index=True)
Hppc,Hcloud=ppc_for(height,height_base,-1,Hfits,'height')
Bppc,Bcloud=ppc_for(bmi,bmi_base,1,Bfits,'bmi')

# Height independent evidence
hm=height[height.Dose_response_class=='mirror'];ha=height[height.Dose_response_class=='asymmetric']
size_tests=[]
for var,transform in [('CNV_length_Mb',np.log),('Gene_count',lambda x:np.log1p(x))]:
 a=transform(ha[var].to_numpy());m=transform(hm[var].to_numpy())
 size_tests.append({'variable':var,'asymmetric_mean':ha[var].mean(),'mirror_mean':hm[var].mean(),'asymmetric_median':ha[var].median(),'mirror_median':hm[var].median(),'Welch_p':ttest_ind(a,m,equal_var=False).pvalue,'Mann_Whitney_p':mannwhitneyu(ha[var],hm[var],alternative='two-sided').pvalue})
size_tests=pd.DataFrame(size_tests)
# residual score based on fitted mirror expectation and correct lower-left direction
hres=Hpred['mirror']
corr_rows=[]
for subset_name,z in [('all',hres),('mirror',hres[hres.Dose_response_class=='mirror']),('asymmetric',hres[hres.Dose_response_class=='asymmetric'])]:
 for var in ['CNV_length_Mb','Gene_count']:
  r,p=spearmanr(z.directional_shift_score,np.log1p(z[var]))
  corr_rows.append({'subset':subset_name,'variable':var,'Spearman_rho':r,'p_value':p})
Hcorr=pd.DataFrame(corr_rows)
# occupancy simulation preserving actual G
qH=Hfits['single']['q'];G=height.Gene_count.to_numpy(int);n_asym=(height.Dose_response_class=='asymmetric').sum();sim=[]
for _ in range(35000):
 hit=rng.random(len(G))<(1-(1-qH)**G)
 if hit.sum()==n_asym:
  sim.append([G[hit].mean()-G[~hit].mean(),height.CNV_length_Mb.to_numpy()[hit].mean()-height.CNV_length_Mb.to_numpy()[~hit].mean()])
sim=np.array(sim)
Hocc=pd.DataFrame([{
 'single_hit_frequency':qH,'expected_asymmetric_fraction':np.mean(1-(1-qH)**G),'observed_asymmetric_fraction':n_asym/len(G),'conditional_simulations':len(sim),
 'observed_gene_count_difference':ha.Gene_count.mean()-hm.Gene_count.mean(),'simulated_gene_count_difference_mean':sim[:,0].mean(),'simulated_gene_count_difference_2.5%':np.quantile(sim[:,0],.025),'simulated_gene_count_difference_97.5%':np.quantile(sim[:,0],.975),'p_sim_gene_diff_ge_observed':np.mean(sim[:,0]>=(ha.Gene_count.mean()-hm.Gene_count.mean())),
 'observed_length_difference_Mb':ha.CNV_length_Mb.mean()-hm.CNV_length_Mb.mean(),'simulated_length_difference_mean_Mb':sim[:,1].mean(),'simulated_length_difference_2.5%':np.quantile(sim[:,1],.025),'simulated_length_difference_97.5%':np.quantile(sim[:,1],.975),'p_sim_length_diff_ge_observed':np.mean(sim[:,1]>=(ha.CNV_length_Mb.mean()-hm.CNV_length_Mb.mean()))
}])

# Burden predictions
model_order=['mirror','continuous','single']
model_label={'mirror':'Mirror alone','continuous':'Mirror + continuous','single':'Mirror + single hit'}
def burden_table(dat,base,sign,fits,trait,target,target_units,trait_sd_native):
 rows=[];Gobs=dat.Gene_count.to_numpy(float)
 for key in model_order:
  d,u=burden_values(base,sign,key,fits[key],Gobs);dmean=float(d.mean());umean=float(u.mean())
  d_native=dmean*trait_sd_native;u_native=umean*trait_sd_native
  density=target/d_native if d_native!=0 else np.nan
  pred_u=umean/dmean*target if dmean!=0 else np.nan
  rows.append({'model':model_label[key],'DEL_effect_per_affected_gene_SD':dmean,'DUP_effect_per_affected_gene_SD':umean,'DEL_effect_per_affected_gene_native':d_native,'DUP_effect_per_affected_gene_native':u_native,'Mace_DEL_target_per_Mb':target,'Mace_units':target_units,'effective_dosage_genes_per_Mb_to_match_DEL':density,'predicted_DUP_effect_per_Mb_at_matched_density':pred_u})
 return pd.DataFrame(rows)
Hbur=burden_table(height,height_base,-1,Hfits,'height',-.41,'cm/Mb',7.99)
Bbur=burden_table(bmi,bmi_base,1,Bfits,'BMI',.14,'kg/m2/Mb',4.7)
# sensitivity by G
sens=[]
for trait,base,sgn,fits,target in [('height',height_base,-1,Hfits,-.41),('BMI',bmi_base,1,Bfits,.14)]:
 for g in range(2,101):
  for key in model_order:
   d,u=burden_values(base,sgn,key,fits[key],[g]);sens.append({'trait':trait,'Gene_count':g,'model':model_label[key],'DEL_SD_per_gene':d[0],'DUP_SD_per_gene':u[0],'predicted_DUP_per_Mb_after_DEL_calibration':target*u[0]/d[0]})
burden_sens=pd.DataFrame(sens)

# Save tables
for path,df in [
 ('height/outputs/height_mirror_parameters.csv',pd.DataFrame([height_base])),('bmi/outputs/bmi_mirror_parameters.csv',pd.DataFrame([bmi_base])),
 ('height/outputs/height_single_locus_model_comparison.csv',Hcomp),('bmi/outputs/bmi_single_locus_model_comparison.csv',Bcomp),
 ('height/outputs/height_conditional_LOO.csv',Hloo),('bmi/outputs/bmi_conditional_LOO.csv',Bloo),
 ('height/outputs/height_continuous_q_profile.csv',Hqcont),('height/outputs/height_single_hit_q_profile.csv',Hqsingle),('bmi/outputs/bmi_continuous_q_profile.csv',Bqcont),('bmi/outputs/bmi_single_hit_q_profile.csv',Bqsingle),
 ('height/outputs/height_posterior_predictive_checks.csv',Hppc),('bmi/outputs/bmi_posterior_predictive_checks.csv',Bppc),
 ('height/outputs/height_size_gene_tests.csv',size_tests),('height/outputs/height_within_class_gradient_tests.csv',Hcorr),('height/outputs/height_single_hit_occupancy_simulation.csv',Hocc),
 ('height/outputs/height_burden_predictions.csv',Hbur),('bmi/outputs/bmi_burden_predictions.csv',Bbur),('code/burden_gene_count_sensitivity.csv',burden_sens),
 ('height/outputs/height_locus_predictions_mirror.csv',Hpred['mirror']),('height/outputs/height_locus_predictions_continuous.csv',Hpred['continuous']),('height/outputs/height_locus_predictions_single_hit.csv',Hpred['single']),
 ('bmi/outputs/bmi_locus_predictions_mirror.csv',Bpred['mirror']),('bmi/outputs/bmi_locus_predictions_continuous.csv',Bpred['continuous']),('bmi/outputs/bmi_locus_predictions_single_hit.csv',Bpred['single'])]:
    df.to_csv(ROOT/path,index=False)

# Summary JSON
summary={
 'height':{'mirror':height_base,'models':{k:{kk:vv for kk,vv in v.items() if kk!='raw_par'} for k,v in Hfits.items()},'model_comparison':Hcomp.to_dict('records'),'burden':Hbur.to_dict('records'),'occupancy':Hocc.to_dict('records')[0]},
 'BMI':{'mirror':bmi_base,'models':{k:{kk:vv for kk,vv in v.items() if kk!='raw_par'} for k,v in Bfits.items()},'model_comparison':Bcomp.to_dict('records'),'burden':Bbur.to_dict('records')},
 'direction_correction':{'height':'asymmetry shifts DEL and DUP effects toward the lower-left','BMI':'asymmetry shifts DEL and DUP effects toward the upper-right'}
}
(ROOT/'code/analysis_summary.json').write_text(json.dumps(summary,indent=2,default=float))

# ---------------- figures ----------------
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'figure.dpi':160})
BLUE='#286090';RED='#b74d45';GREEN='#4f8a5b';PURPLE='#7563a8';GREY='#777777';LIGHT='#dddddd';ORANGE='#c47b28';BLACK='#222222'

def savefig(fig,pathbase):
    fig.savefig(str(pathbase)+'.png',dpi=300,bbox_inches='tight');fig.savefig(str(pathbase)+'.pdf',bbox_inches='tight');plt.close(fig)

def mirror_line(base,gmin,gmax):
    gs=np.linspace(gmin,gmax,100);mean_gene=base['p']*base['sign_mean']*base['m'];return gs*mean_gene,-base['rho']*gs*mean_gene

def panel_label(ax,s):ax.text(-.12,1.05,s,transform=ax.transAxes,fontweight='bold',fontsize=12,va='top')

# Height Figure 1
fig,axs=plt.subplots(2,2,figsize=(10,8))
ax=axs[0,0];z=height_mirror
ax.errorbar(z.Deletion_beta,z.Duplication_beta,xerr=z.Deletion_SE,yerr=z.Duplication_SE,fmt='o',ms=4,color=BLUE,ecolor='#9ab3c9',elinewidth=.7,capsize=1)
xline,yline=mirror_line(height_base,z.Gene_count.min(),z.Gene_count.max());ax.plot(xline,yline,color=BLACK,lw=2,label='Fitted mirror expectation')
ax.axhline(0,color=LIGHT,lw=.7);ax.axvline(0,color=LIGHT,lw=.7);ax.set(xlabel='Deletion effect on height (SD)',ylabel='Duplication effect on height (SD)',title='Mirror loci follow a reciprocal dose-response')
ax.legend(frameon=False,fontsize=8);panel_label(ax,'A')
ax=axs[0,1]
# profile from prior exact mirror fit file if available
prof=load_frequency_effect_profile(H_SRC,'height_frequency_effect_profile.csv',height_mirror,0.30,ROOT/'height/outputs/height_frequency_effect_profile.csv')
ax.plot(100*prof.mirror_gene_frequency,prof.mean_abs_gene_effect_SD,color=BLUE,lw=2);ax.scatter([100*height_base['p']],[height_base['m']],color=RED,zorder=3)
ax.set(xlabel='Mirror-contributor frequency (%)',ylabel='Mean absolute gene effect (SD)',title='Frequency and effect size trade off')
ax.text(.98,.95,f"Best point: {100*height_base['p']:.1f}%\nMean effect: {height_base['m']:.3f} SD",transform=ax.transAxes,ha='right',va='top',bbox=dict(boxstyle='round',fc='white',ec=LIGHT));panel_label(ax,'B')
ax=axs[1,0]
obs=empirical(height_mirror);sim,_=simulate_model(height_mirror,height_base,-1,'mirror',Hfits['mirror'],1500,SEED+99)
vals=['DEL mean','DUP mean','Correlation'];oi=[obs[0],obs[1],obs[2]];sm=[sim[:,0].mean(),sim[:,1].mean(),sim[:,2].mean()];lo=[np.quantile(sim[:,0],.025),np.quantile(sim[:,1],.025),np.quantile(sim[:,2],.025)];hi=[np.quantile(sim[:,0],.975),np.quantile(sim[:,1],.975),np.quantile(sim[:,2],.975)]
x=np.arange(3);ax.errorbar(x,sm,yerr=[np.array(sm)-np.array(lo),np.array(hi)-np.array(sm)],fmt='o',color=BLUE,capsize=4,label='Simulated 95% interval');ax.scatter(x,oi,color=RED,zorder=4,label='Observed');ax.axhline(0,color=LIGHT,lw=.7);ax.set_xticks(x,vals);ax.set_title('Mirror model reproduces central tendencies');ax.legend(frameon=False,fontsize=8);panel_label(ax,'C')
ax=axs[1,1];ax.axis('off')
text=(f"Height mirror model\n\nKosmicki sign prior: 70% of deletion/LoF effects decrease height.\n\nFitted contributor frequency: {100*height_base['p']:.1f}%\nMean absolute contributor effect: {height_base['m']:.3f} SD\nExpected contributors in a 22-gene CNV: {22*height_base['p']:.1f}\nDuplication/deletion magnitude ratio: {height_base['rho']:.2f}\n\nInterpretation: multiple modest effects partially cancel.\nDuplication effects are substantially weaker than deletion effects.")
ax.text(.03,.97,text,va='top',fontsize=10,bbox=dict(boxstyle='round,pad=.5',fc='#f7f7f7',ec=LIGHT));panel_label(ax,'D')
fig.suptitle('Height: mirror-gene model fitted to Sacks et al. mirror CNVs',fontsize=14,fontweight='bold');fig.tight_layout(rect=[0,0,1,.96]);savefig(fig,ROOT/'height/figures/Figure_1_height_mirror_model')

# Height Figure 2
fig,axs=plt.subplots(2,3,figsize=(13,8))
ax=axs[0,0]
cols=height.Dose_response_class.map({'mirror':BLUE,'asymmetric':RED})
ax.scatter(height.Deletion_beta,height.Duplication_beta,c=cols,s=35)
ax.plot(xline,yline,color=BLACK,lw=1.5);ax.add_patch(FancyArrowPatch((.05,.05),(-.22,-.22),arrowstyle='->',mutation_scale=15,color=RED,lw=2));ax.text(-.18,-.08,'Height-biased\nasymmetry',color=RED)
ax.axhline(0,color=LIGHT,lw=.7);ax.axvline(0,color=LIGHT,lw=.7);ax.set(xlabel='DEL effect (SD)',ylabel='DUP effect (SD)',title='Observed mirror and asymmetric loci');panel_label(ax,'A')
# model predicted expectations
for kk,ax,letter in [('mirror',axs[0,1],'B'),('continuous',axs[0,2],'C'),('single',axs[1,0],'D')]:
 p=Hpred[kk];ax.scatter(height.Deletion_beta,height.Duplication_beta,c=cols,s=18,alpha=.45);ax.scatter(p.model_expected_DEL,p.model_expected_DUP,color={'mirror':GREY,'continuous':GREEN,'single':PURPLE}[kk],s=25,label='Model expectation')
 for i in range(len(p)):
  ax.plot([p.model_expected_DEL.iloc[i],p.Deletion_beta.iloc[i]],[p.model_expected_DUP.iloc[i],p.Duplication_beta.iloc[i]],color=LIGHT,lw=.5)
 ax.axhline(0,color=LIGHT,lw=.7);ax.axvline(0,color=LIGHT,lw=.7);ax.set(xlabel='DEL effect (SD)',ylabel='DUP effect (SD)',title={'mirror':'Mirror alone leaves systematic lower-left residuals','continuous':'Continuous shifts grow with gene count','single':'Single hits create a second displaced component'}[kk]);panel_label(ax,letter)
ax=axs[1,1]
order=['Mirror alone','Mirror + continuous directional shift','Mirror + sparse single hit'];cc=Hcomp.set_index('model').loc[order]
x=np.arange(3);ax.bar(x,cc.delta_conditional_AIC,color=[GREY,GREEN,PURPLE],alpha=.85);ax.set_xticks(x,['Mirror','Continuous','Single hit']);ax.set_ylabel('Delta conditional AIC');ax.set_title('Single hit has the best conditional AIC');ax2=ax.twinx();ax2.plot(x,cc.conditional_LOO_log_predictive_density,'o--',color=BLACK);ax2.set_ylabel('Conditional LOO log predictive density');panel_label(ax,'E')
ax=axs[1,2]
grid=np.arange(1,70);ax.plot(grid,1-(1-Hfits['single']['q'])**grid,color=PURPLE,lw=2,label='Single-hit occupancy')
ax.scatter(height.Gene_count,(height.Dose_response_class=='asymmetric').astype(int),c=cols,s=25,alpha=.65);ax.set(xlabel='Genes in CNV',ylabel='Probability / observed class',title='Rare hits predict more asymmetry in larger CNVs');ax.text(.98,.08,f"q={100*Hfits['single']['q']:.2f}%\nExpected asymmetric={100*Hocc.expected_asymmetric_fraction.iloc[0]:.1f}%\nObserved={100*Hocc.observed_asymmetric_fraction.iloc[0]:.1f}%",transform=ax.transAxes,ha='right',bbox=dict(boxstyle='round',fc='white',ec=LIGHT));ax.legend(frameon=False,fontsize=8);panel_label(ax,'F')
fig.suptitle('Height: corrected lower-left asymmetry models fitted to Sacks et al. loci',fontsize=14,fontweight='bold');fig.tight_layout(rect=[0,0,1,.96]);savefig(fig,ROOT/'height/figures/Figure_2_height_single_locus_comparison')

# Height Figure 3 burden
fig,axs=plt.subplots(2,2,figsize=(10,8))
ax=axs[0,0];x=np.arange(3);w=.34;ax.bar(x-w/2,Hbur.DEL_effect_per_affected_gene_SD,w,label='DEL',color=RED);ax.bar(x+w/2,Hbur.DUP_effect_per_affected_gene_SD,w,label='DUP',color=BLUE);ax.axhline(0,color=BLACK,lw=.7);ax.set_xticks(x,['Mirror','Continuous','Single hit']);ax.set_ylabel('Expected height effect per affected gene (SD)');ax.set_title('Asymmetry reverses the DUP-burden direction');ax.legend(frameon=False);panel_label(ax,'A')
ax=axs[0,1];ax.bar(x,Hbur.predicted_DUP_effect_per_Mb_at_matched_density,color=[GREY,GREEN,PURPLE]);ax.axhline(0,color=BLACK,lw=.7);ax.set_xticks(x,['Mirror','Continuous','Single hit']);ax.set_ylabel('Predicted DUP height effect (cm/Mb)');ax.set_title('After calibrating each model to Macé DEL burden');ax.text(.03,.05,'Macé DEL target = -0.41 cm/Mb',transform=ax.transAxes);panel_label(ax,'B')
ax=axs[1,0];z=burden_sens[burden_sens.trait=='height'];
for model,col in [('Mirror alone',GREY),('Mirror + continuous',GREEN),('Mirror + single hit',PURPLE)]:
 zz=z[z.model==model];ax.plot(zz.Gene_count,zz.predicted_DUP_per_Mb_after_DEL_calibration,label=model,color=col,lw=2)
ax.axhline(0,color=BLACK,lw=.7);ax.set(xlabel='Genes per CNV',ylabel='Predicted DUP effect after DEL calibration (cm/Mb)',title='Burden direction is stable across CNV sizes');ax.legend(frameon=False,fontsize=8);panel_label(ax,'C')
ax=axs[1,1];ax.axis('off')
text=("Burden evidence\n\nMirror alone predicts that deletion burden lowers height, but duplication\nburden raises height. This conflicts with genome-wide burden evidence.\n\nBoth corrected asymmetry models shift duplication burden below zero\nwhile retaining a negative deletion effect.\n\nBurden alone does not distinguish continuous from single-hit asymmetry.\nSingle hit is favored when combined with the discrete locus classes,\ntheir ~30% prevalence, and their larger size and gene content.")
ax.text(.03,.97,text,va='top',fontsize=10,bbox=dict(boxstyle='round,pad=.5',fc='#f7f7f7',ec=LIGHT));panel_label(ax,'D')
fig.suptitle('Height: propagation from recurrent CNVs to genome-wide burden',fontsize=14,fontweight='bold');fig.tight_layout(rect=[0,0,1,.96]);savefig(fig,ROOT/'height/figures/Figure_3_height_burden_comparison')

# BMI Figure 1
fig,axs=plt.subplots(2,2,figsize=(10,8))
ax=axs[0,0];ax.errorbar(bmi.Deletion_beta,bmi.Duplication_beta,xerr=bmi.Deletion_SE,yerr=bmi.Duplication_SE,fmt='o',ms=4,color=BLUE,ecolor='#9ab3c9',elinewidth=.7,capsize=1)
xlineB,ylineB=mirror_line(bmi_base,bmi.Gene_count.min(),bmi.Gene_count.max());ax.plot(xlineB,ylineB,color=BLACK,lw=2);ax.add_patch(FancyArrowPatch((0,0),(.20,.20),arrowstyle='->',mutation_scale=15,color=RED,lw=2));ax.text(.08,.24,'BMI-biased\nasymmetry',color=RED);ax.axhline(0,color=LIGHT,lw=.7);ax.axvline(0,color=LIGHT,lw=.7);ax.set(xlabel='Deletion effect on BMI (SD)',ylabel='Duplication effect on BMI (SD)',title='Reciprocal effects with an upper-right shift');panel_label(ax,'A')
ax=axs[0,1];prof=load_frequency_effect_profile(B_SRC,'BMI_frequency_effect_profile.csv',bmi,11/14,ROOT/'bmi/outputs/BMI_frequency_effect_profile.csv');ax.plot(100*prof.mirror_gene_frequency,prof.mean_abs_gene_effect_SD,color=BLUE,lw=2);ax.scatter([100*bmi_base['p']],[bmi_base['m']],color=RED);ax.set(xlabel='Mirror-contributor frequency (%)',ylabel='Mean absolute gene effect (SD)',title='BMI frequency-effect ridge');ax.text(.98,.95,f"Best point: {100*bmi_base['p']:.1f}%\nMean effect: {bmi_base['m']:.3f} SD",transform=ax.transAxes,ha='right',va='top',bbox=dict(boxstyle='round',fc='white',ec=LIGHT));panel_label(ax,'B')
ax=axs[1,0]
ak=pd.read_csv(ROOT/'bmi/inputs/akbari_2021_pLoF_direction_reference.csv');counts=ak.Direction.value_counts();ax.bar(['BMI increasing','BMI decreasing'],[counts.get('increases BMI',0),counts.get('decreases BMI',0)],color=[RED,BLUE]);ax.set_ylabel('Explicitly pLoF-tested genes');ax.set_title('Akbari supplies the LoF direction prior');ax.text(.5,.9,'Posterior mean used: 11/14 = 78.6%',transform=ax.transAxes,ha='center');panel_label(ax,'C')
ax=axs[1,1];ax.axis('off')
text=(f"BMI mirror model\n\nAkbari-informed probability that LoF increases BMI: 78.6%.\n\nFitted contributor frequency: {100*bmi_base['p']:.1f}%\nMean absolute contributor effect: {bmi_base['m']:.3f} SD\nExpected contributors in the average CNV: {bmi_base['p']*bmi.Gene_count.mean():.2f}\nDuplication/deletion magnitude ratio: {bmi_base['rho']:.2f}\n\nObservation: the mirror component explains the strong negative\nDEL-DUP relationship, but predicts a more negative average DUP\neffect than the observed near-zero mean.")
ax.text(.03,.97,text,va='top',fontsize=10,bbox=dict(boxstyle='round,pad=.5',fc='#f7f7f7',ec=LIGHT));panel_label(ax,'D')
fig.suptitle('BMI: mirror-gene model fitted to exact Sacks et al. reciprocal CNVs',fontsize=14,fontweight='bold');fig.tight_layout(rect=[0,0,1,.96]);savefig(fig,ROOT/'bmi/figures/Figure_1_BMI_mirror_model')

# BMI Figure 2
fig,axs=plt.subplots(2,3,figsize=(13,8))
ax=axs[0,0];r=Bpred['mirror'];sc=ax.scatter(r.residual_DEL_from_mirror,r.residual_DUP_from_mirror,c=r.directional_shift_score,cmap='coolwarm',s=40);ax.axhline(0,color=LIGHT);ax.axvline(0,color=LIGHT);ax.add_patch(FancyArrowPatch((0,0),(.2,.2),arrowstyle='->',mutation_scale=15,color=RED,lw=2));ax.set(xlabel='DEL residual from mirror',ylabel='DUP residual from mirror',title='Both DEL and DUP buffering point upper-right');fig.colorbar(sc,ax=ax,label='Upper-right shift score');panel_label(ax,'A')
for kk,ax,letter in [('mirror',axs[0,1],'B'),('continuous',axs[0,2],'C'),('single',axs[1,0],'D')]:
 p=Bpred[kk];ax.scatter(bmi.Deletion_beta,bmi.Duplication_beta,color=BLUE,s=18,alpha=.45);ax.scatter(p.model_expected_DEL,p.model_expected_DUP,color={'mirror':GREY,'continuous':GREEN,'single':PURPLE}[kk],s=25)
 for i in range(len(p)):ax.plot([p.model_expected_DEL.iloc[i],p.Deletion_beta.iloc[i]],[p.model_expected_DUP.iloc[i],p.Duplication_beta.iloc[i]],color=LIGHT,lw=.5)
 ax.axhline(0,color=LIGHT);ax.axvline(0,color=LIGHT);ax.set(xlabel='DEL effect (SD)',ylabel='DUP effect (SD)',title={'mirror':'Mirror alone underestimates the upper-right shift','continuous':'Many tiny shifts move the full cloud','single':'Sparse hits create a displaced mixture'}[kk]);panel_label(ax,letter)
ax=axs[1,1];order=['Mirror alone','Mirror + continuous directional shift','Mirror + sparse single hit'];cc=Bcomp.set_index('model').loc[order];x=np.arange(3);ax.bar(x,cc.delta_conditional_AIC,color=[GREY,GREEN,PURPLE]);ax.set_xticks(x,['Mirror','Continuous','Single hit']);ax.set_ylabel('Delta conditional AIC');ax.set_title('Continuous has the best conditional fit');ax2=ax.twinx();ax2.plot(x,cc.conditional_LOO_log_predictive_density,'o--',color=BLACK);ax2.set_ylabel('Conditional LOO log predictive density');panel_label(ax,'E')
ax=axs[1,2];ax.plot(100*Bqcont.q,Bqcont.shift_magnitude,color=GREEN,lw=2,label='Continuous');ax.plot(100*Bqsingle.q,Bqsingle.shift_magnitude,color=PURPLE,lw=2,label='Single hit');ax.axvline(100*Bfits['continuous']['q'],color=GREEN,ls='--');ax.axvline(100*Bfits['single']['q'],color=PURPLE,ls='--');ax.set(xlabel='Asymmetry-gene frequency (%)',ylabel='Fitted shift magnitude (SD)',title='Frequency and magnitude are weakly identified');ax.legend(frameon=False);panel_label(ax,'F')
fig.suptitle('BMI: corrected upper-right asymmetry models fitted to Sacks et al. loci',fontsize=14,fontweight='bold');fig.tight_layout(rect=[0,0,1,.96]);savefig(fig,ROOT/'bmi/figures/Figure_2_BMI_single_locus_comparison')

# BMI Figure 3 burden
fig,axs=plt.subplots(2,2,figsize=(10,8))
ax=axs[0,0];x=np.arange(3);w=.34;ax.bar(x-w/2,Bbur.DEL_effect_per_affected_gene_SD,w,label='DEL',color=RED);ax.bar(x+w/2,Bbur.DUP_effect_per_affected_gene_SD,w,label='DUP',color=BLUE);ax.axhline(0,color=BLACK,lw=.7);ax.set_xticks(x,['Mirror','Continuous','Single hit']);ax.set_ylabel('Expected BMI effect per affected gene (SD)');ax.set_title('Corrected asymmetry moves both burdens toward obesity');ax.legend(frameon=False);panel_label(ax,'A')
ax=axs[0,1];ax.bar(x,Bbur.predicted_DUP_effect_per_Mb_at_matched_density,color=[GREY,GREEN,PURPLE]);ax.axhline(0,color=BLACK,lw=.7);ax.set_xticks(x,['Mirror','Continuous','Single hit']);ax.set_ylabel('Predicted DUP BMI effect (kg/m2/Mb)');ax.set_title('After calibrating each model to Macé DEL burden');ax.text(.03,.05,'Macé DEL target = +0.14 kg/m2/Mb',transform=ax.transAxes);panel_label(ax,'B')
ax=axs[1,0];z=burden_sens[burden_sens.trait=='BMI'];
for model,col in [('Mirror alone',GREY),('Mirror + continuous',GREEN),('Mirror + single hit',PURPLE)]:
 zz=z[z.model==model];ax.plot(zz.Gene_count,zz.predicted_DUP_per_Mb_after_DEL_calibration,label=model,color=col,lw=2)
ax.axhline(0,color=BLACK,lw=.7);ax.set(xlabel='Genes per CNV',ylabel='Predicted DUP effect after DEL calibration (kg/m2/Mb)',title='Upper-right shifts reverse the mirror-only DUP prediction');ax.legend(frameon=False,fontsize=8);panel_label(ax,'C')
ax=axs[1,1];ax.axis('off')
text=("Burden evidence\n\nMacé reported that deletion burden increases BMI. Mirror alone matches\nthat direction but predicts that duplication burden should lower BMI.\n\nBoth corrected upper-right models move predicted duplication burden\nto approximately zero or positive, consistent with adverse adiposity burden.\n\nBurden data do not require a discrete class. At the single-locus level,\ncontinuous asymmetry has slightly better AIC and LOO prediction, and\nthe BMI cloud does not divide into two curves.")
ax.text(.03,.97,text,va='top',fontsize=10,bbox=dict(boxstyle='round,pad=.5',fc='#f7f7f7',ec=LIGHT));panel_label(ax,'D')
fig.suptitle('BMI: propagation from recurrent CNVs to genome-wide burden',fontsize=14,fontweight='bold');fig.tight_layout(rect=[0,0,1,.96]);savefig(fig,ROOT/'bmi/figures/Figure_3_BMI_burden_comparison')

# Combined summary figure
fig,axs=plt.subplots(1,2,figsize=(12,5))
for ax,trait,comp,conclusion in [(axs[0],'Height',Hcomp,'Sparse single hit'),(axs[1],'BMI',Bcomp,'Continuous shifted mirror')]:
 cc=comp.set_index('model').loc[['Mirror alone','Mirror + continuous directional shift','Mirror + sparse single hit']];x=np.arange(3);ax.bar(x,cc.delta_conditional_AIC,color=[GREY,GREEN,PURPLE]);ax.set_xticks(x,['Mirror','Continuous','Single hit']);ax.set_ylabel('Delta conditional AIC');ax.set_title(f'{trait}: most parsimonious = {conclusion}',fontweight='bold');ax.text(.5,.94,'Lower is better',transform=ax.transAxes,ha='center')
fig.suptitle('Corrected directional model comparison across traits',fontsize=14,fontweight='bold');fig.tight_layout(rect=[0,0,1,.92]);savefig(fig,ROOT/'report/Figure_summary_model_comparison')

# ---------------- report DOCX ----------------
def shade_cell(cell,fill):
    tcPr=cell._tc.get_or_add_tcPr();shd=OxmlElement('w:shd');shd.set(qn('w:fill'),fill);tcPr.append(shd)
def set_cell_text(cell,text,bold=False,size=8):
    cell.text='';p=cell.paragraphs[0];r=p.add_run(str(text));r.bold=bold;r.font.size=Pt(size)
def add_table(doc,df,cols=None,digits=3):
    if cols is None:cols=list(df.columns)
    t=doc.add_table(rows=1,cols=len(cols));t.alignment=WD_TABLE_ALIGNMENT.CENTER;t.style='Table Grid'
    for j,c in enumerate(cols):set_cell_text(t.rows[0].cells[j],c,True,8);shade_cell(t.rows[0].cells[j],'D9EAF7')
    for _,row in df.iterrows():
        cells=t.add_row().cells
        for j,c in enumerate(cols):
            v=row[c]
            if isinstance(v,(float,np.floating)):
                txt='' if pd.isna(v) else f'{v:.{digits}f}'
            else:txt=str(v)
            set_cell_text(cells[j],txt,False,8);cells[j].vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
    doc.add_paragraph()
    return t
def add_fig(doc,path,caption,width=6.6):
    p=doc.add_paragraph();p.alignment=WD_ALIGN_PARAGRAPH.CENTER;p.add_run().add_picture(str(path),width=Inches(width))
    c=doc.add_paragraph(caption);c.alignment=WD_ALIGN_PARAGRAPH.CENTER;c.runs[0].italic=True;c.runs[0].font.size=Pt(9)
def add_callout(doc,text):
    t=doc.add_table(rows=1,cols=1);t.style='Table Grid';shade_cell(t.cell(0,0),'EDF5E9');p=t.cell(0,0).paragraphs[0];r=p.add_run(text);r.bold=True;r.font.size=Pt(10);doc.add_paragraph()

doc=Document();sec=doc.sections[0];sec.top_margin=Inches(.55);sec.bottom_margin=Inches(.55);sec.left_margin=Inches(.6);sec.right_margin=Inches(.6)
styles=doc.styles;styles['Normal'].font.name='Arial';styles['Normal'].font.size=Pt(9.5)
for s in ['Title','Heading 1','Heading 2','Heading 3']:
    styles[s].font.name='Arial'
styles['Title'].font.size=Pt(22);styles['Heading 1'].font.size=Pt(15);styles['Heading 2'].font.size=Pt(12)
p=doc.add_paragraph();p.style='Title';p.alignment=WD_ALIGN_PARAGRAPH.CENTER;p.add_run('Corrected directional models of reciprocal CNV effects on height and BMI')
p=doc.add_paragraph();p.alignment=WD_ALIGN_PARAGRAPH.CENTER;p.add_run('Mirror alone versus continuous generalized dysfunction versus sparse single-hit asymmetry').italic=True
p=doc.add_paragraph();p.alignment=WD_ALIGN_PARAGRAPH.CENTER;p.add_run('Exact Sacks et al. reciprocal-CNV effects; burden comparison with Macé et al. and Auwerx et al.')

doc.add_heading('Executive summary',1)
add_callout(doc,'Direction correction: height-biased asymmetry was constrained to move deletion and duplication effects toward the lower-left of the reciprocal-effect plot, whereas BMI-biased asymmetry was constrained to move both effects toward the upper-right.')
doc.add_paragraph('For height, mirror alone explains the reciprocal trend among mirror loci but fails for the full locus set and predicts a beneficial duplication-burden effect. Both asymmetry models correct the burden direction. The sparse single-hit model is the most parsimonious overall because it has the best conditional AIC and leave-one-locus-out prediction and, uniquely, links the same approximately 2% gene frequency to the discrete asymmetric class, its approximately 30% prevalence, and its enrichment in larger, gene-rich CNVs.')
doc.add_paragraph('For BMI, the mirror component explains the strong negative deletion-duplication relationship but predicts a more protective duplication effect than observed. Correctly oriented upper-right shifts improve the locus fit and move predicted duplication burden toward zero or obesity. The continuous model has the best conditional AIC and leave-one-locus-out prediction, while the scatter remains a single continuous cloud rather than two curves. BMI is therefore most parsimoniously described by an obesity-shifted mirror model with broadly distributed small upper-right effects, not a distinct sparse asymmetric class.')

doc.add_heading('1. Data and corrected model structure',1)
doc.add_paragraph(f'The single-locus analyses used {len(height)} independent reciprocal height loci ({len(height_mirror)} mirror and {(height.Dose_response_class=="asymmetric").sum()} asymmetric) and {len(bmi)} reciprocal BMI loci from the exact Sacks et al. supplementary statistics. The overlapping 22q11.2 A-D interval was excluded. Each locus contributed its deletion and duplication beta, standard errors, and annotated gene count.')
doc.add_paragraph('The mirror component treated a CNV effect as the additive sum of latent gene effects. Kosmicki informed the height sign prior (70% of deletion/LoF effects decrease height); Akbari informed the BMI sign prior (posterior mean 78.6% of pLoF effects increase BMI). The ordinary contributor frequency, mean absolute gene effect, and duplication-to-deletion magnitude ratio were estimated from Sacks CNVs rather than imported from the large, ascertainment-selected WES effects.')
doc.add_paragraph('After calibrating the mirror component, three models were compared on the complete locus set: mirror alone; mirror plus continuous asymmetry, in which many small directional effects add with CNV gene count; and mirror plus a sparse single hit, in which the presence of at least one rare gene creates a fixed directional displacement. In the corrected bivariate parameterization, all asymmetry shifts were restricted to the lower-left for height and the upper-right for BMI. The shift angle was estimated within that quadrant, so both deletion and duplication components remained correctly signed.')
doc.add_paragraph('Model comparisons below are conditional on the calibrated mirror component. AIC/BIC evaluate fit to the Sacks single-locus data; conditional leave-one-locus-out values test predictive stability. Burden analyses are presented separately: each fitted gene-level model was propagated to expected deletion and duplication burden, calibrated to the Macé deletion coefficient, and checked against the genome-wide direction reported by Auwerx.')

doc.add_heading('2. Height: fit to Sacks et al. single-locus results',1)
add_fig(doc,ROOT/'height/figures/Figure_1_height_mirror_model.png','Figure 1. Height mirror model. The observed mirror loci follow a negative reciprocal relationship; the model infers multiple modest contributors and weaker duplication effects.',6.8)
doc.add_paragraph(f"The fitted height mirror model estimated a contributor frequency of {100*height_base['p']:.1f}% and a mean absolute effect of {height_base['m']:.3f} SD, corresponding to approximately {22*height_base['p']:.1f} contributors in a 22-gene CNV. The duplication-to-deletion magnitude ratio was {height_base['rho']:.2f}. The exact percentage and per-gene magnitude lie on a broad frequency-effect ridge, but the robust inference is that mirror CNVs aggregate multiple modest effects and that duplication effects are substantially weaker than deletion effects.")
add_fig(doc,ROOT/'height/figures/Figure_2_height_single_locus_comparison.png','Figure 2. Height single-locus comparison. All asymmetry shifts are lower-left. The single-hit model creates a second displaced component and predicts increasing occupancy in larger CNVs.',7.1)
add_table(doc,Hcomp[['model','q','shift_magnitude','shift_DEL','shift_DUP','conditional_AIC','conditional_BIC','conditional_LOO_log_predictive_density','delta_conditional_AIC','delta_conditional_LOO']],digits=3)
doc.add_paragraph(f"The continuous model estimated q={100*Hfits['continuous']['q']:.1f}% and a mean lower-left vector magnitude of {Hfits['continuous']['shift_magnitude']:.3f} SD. The sparse model estimated q={100*Hfits['single']['q']:.2f}% and a much larger vector magnitude of {Hfits['single']['shift_magnitude']:.3f} SD. The fitted direction was predominantly downward in both models but retained the correct negative deletion component. The single-hit model improved conditional AIC by {Hcomp.loc[Hcomp.model=='Mirror + continuous directional shift','conditional_AIC'].iloc[0]-Hcomp.loc[Hcomp.model=='Mirror + sparse single hit','conditional_AIC'].iloc[0]:.2f} units relative to the continuous model and had a {Hcomp.loc[Hcomp.model=='Mirror + sparse single hit','conditional_LOO_log_predictive_density'].iloc[0]-Hcomp.loc[Hcomp.model=='Mirror + continuous directional shift','conditional_LOO_log_predictive_density'].iloc[0]:.2f}-unit advantage in conditional leave-one-locus-out prediction.")
doc.add_heading('Specific evidence favoring continuous or single-hit height asymmetry',2)
doc.add_paragraph('Evidence compatible with a continuous mechanism: a continuous additive model fits the reciprocal-effect cloud much better than mirror alone, and its small effects can produce negative deletion and duplication burden. Thus, cumulative generalized dysfunction remains biologically and statistically plausible.')
doc.add_paragraph(f"Evidence favoring a single hit: the sparse model gives the best conditional information criteria and predictive score; q={100*Hfits['single']['q']:.2f}% predicts an asymmetric-locus prevalence of {100*Hocc.expected_asymmetric_fraction.iloc[0]:.1f}%, close to the observed {100*Hocc.observed_asymmetric_fraction.iloc[0]:.1f}%; and the same occupancy mechanism explains why asymmetric CNVs contain more genes and span larger intervals. Within mirror and asymmetric classes, however, the lower-left shift score does not increase with gene count or length, arguing against a simple within-class accumulation gradient.")
add_table(doc,size_tests,digits=4);add_table(doc,Hcorr,digits=4)

doc.add_heading('3. Height: propagation to Macé/Auwerx burden results',1)
add_fig(doc,ROOT/'height/figures/Figure_3_height_burden_comparison.png','Figure 3. Height burden predictions. Mirror alone predicts positive duplication burden; both lower-left asymmetry models reverse it.',6.8)
add_table(doc,Hbur,digits=4)
doc.add_paragraph('The burden analysis is separate from the Sacks single-locus fit. Each model was calibrated to the Macé deletion coefficient of -0.41 cm per deleted Mb by allowing an effective density of dosage-responsive genes per Mb. At that same density, mirror alone predicts that duplication burden increases height, whereas both lower-left asymmetry models predict a negative duplication effect. This directional correction is consistent with the adverse genome-wide CNV-burden pattern reported by Auwerx. Burden therefore rejects mirror alone but does not by itself distinguish continuous from single-hit asymmetry.')
doc.add_paragraph('Height conclusion: considering burden and single-locus data together, the sparse single-hit model is most parsimonious. Both asymmetry models solve the burden-direction problem, but only the sparse model also naturally accounts for the distinct asymmetric group, its prevalence, and its preferential occurrence in larger and gene-rich CNVs.')

doc.add_page_break();doc.add_heading('4. BMI: fit to Sacks et al. single-locus results',1)
add_fig(doc,ROOT/'bmi/figures/Figure_1_BMI_mirror_model.png','Figure 4. BMI mirror model. The mirror component captures the negative reciprocal trend but leaves an upper-right displacement.',6.8)
doc.add_paragraph(f"The BMI mirror model estimated a contributor frequency of {100*bmi_base['p']:.1f}%, a mean absolute contributor effect of {bmi_base['m']:.3f} SD, and {bmi_base['p']*bmi.Gene_count.mean():.2f} contributors in the average CNV. Duplication effects were {bmi_base['rho']:.2f} times the deletion magnitude. The low-frequency/large-effect point estimate is weakly identified and should be interpreted as an effective architecture rather than a literal count. The important observation is that the mirror component explains the strong negative correlation but predicts a more negative average duplication effect than observed.")
add_fig(doc,ROOT/'bmi/figures/Figure_2_BMI_single_locus_comparison.png','Figure 5. BMI single-locus comparison. Corrected asymmetry moves both deletion and duplication residuals upper-right; continuous and single-hit models are compared on equal directional footing.',7.1)
add_table(doc,Bcomp[['model','q','shift_magnitude','shift_DEL','shift_DUP','conditional_AIC','conditional_BIC','conditional_LOO_log_predictive_density','delta_conditional_AIC','delta_conditional_LOO']],digits=3)
doc.add_paragraph(f"The corrected continuous model estimated a very diffuse process: q reached {100*Bfits['continuous']['q']:.0f}% under the broad 0-80% parameterization, but each upper-right shift was only {Bfits['continuous']['shift_magnitude']:.4f} SD. The more stable quantities are the per-gene shifts, {Bfits['continuous']['q']*Bfits['continuous']['shift_DEL']:.4f} SD for deletion and {Bfits['continuous']['q']*Bfits['continuous']['shift_DUP']:.4f} SD for duplication. The single-hit model estimated q={100*Bfits['single']['q']:.2f}% and a shift magnitude of {Bfits['single']['shift_magnitude']:.3f} SD. Continuous asymmetry had the best conditional AIC and exceeded the single-hit conditional leave-one-locus-out score by {Bcomp.loc[Bcomp.model=='Mirror + continuous directional shift','conditional_LOO_log_predictive_density'].iloc[0]-Bcomp.loc[Bcomp.model=='Mirror + sparse single hit','conditional_LOO_log_predictive_density'].iloc[0]:.2f} log units.")
doc.add_heading('Specific evidence favoring continuous or single-hit BMI asymmetry',2)
doc.add_paragraph('Evidence favoring a continuous mechanism: the loci form one continuous reciprocal cloud rather than two separated curves; residuals from the mirror expectation occur on both the deletion and duplication axes but share an upper-right trait direction; and the continuous model has the best conditional AIC and leave-one-locus-out prediction. The fitted frequency-effect ridge also shows that the data identify a small distributed upper-right load more clearly than a particular event frequency.')
doc.add_paragraph('Evidence compatible with single hits: a sparse upper-right mixture also improves upon mirror alone and can reproduce individual loci with a near-null deletion or duplication effect. However, there is no independent BMI class label, size enrichment, or occupancy result analogous to height, and the single-hit model does not improve predictive fit over the continuous model.')

doc.add_heading('5. BMI: propagation to Macé/Auwerx burden results',1)
add_fig(doc,ROOT/'bmi/figures/Figure_3_BMI_burden_comparison.png','Figure 6. BMI burden predictions. Mirror alone predicts protective duplication burden; both corrected upper-right models move duplication burden to near-zero or positive.',6.8)
add_table(doc,Bbur,digits=4)
doc.add_paragraph('The Macé burden result is again analyzed separately from the Sacks locus fit. Macé reported +0.14 kg/m2 of BMI per deleted Mb. Calibrating each model to this deletion coefficient, mirror alone predicts a negative duplication-burden effect. Both corrected upper-right models reverse that prediction to approximately zero or positive, which is directionally consistent with the increased-adiposity burden signal reported by Auwerx. Because both asymmetry models make this correction, the burden result does not favor a sparse BMI class.')
doc.add_paragraph('BMI conclusion: a shifted mirror model with broadly distributed small upper-right effects is most parsimonious. Sparse hits remain possible for particular loci, but the data do not require a separate asymmetric class and do not identify a stable single-hit frequency.')

doc.add_heading('6. Overall synthesis',1)
add_fig(doc,ROOT/'report/Figure_summary_model_comparison.png','Figure 7. Corrected model comparison. Height favors sparse single-hit asymmetry when all evidence is combined; BMI favors continuous upper-right asymmetry.',6.8)
doc.add_paragraph('Mirror effects are a major component of both traits: deletion and duplication effects are negatively correlated, and duplication effects are weaker than deletion effects. Directional asymmetry is nevertheless required to reconcile recurrent loci with genome-wide burden. The direction differs by trait: lower-left for height and upper-right for BMI.')
doc.add_paragraph('For height, continuous and sparse mechanisms both fit the lower-left displacement and both produce negative deletion and duplication burden. The decisive evidence is cross-scale: the sparse model also explains the discrete 30% asymmetric group and its size/gene-count enrichment without a within-class gradient. It is therefore the most parsimonious height model.')
doc.add_paragraph('For BMI, both corrected mechanisms account for DEL- and DUP-side buffering and improve the burden prediction. The absence of two curves, the slightly better conditional AIC and leave-one-locus-out performance of the continuous model, and the non-identifiability of a sparse event frequency favor a continuous obesity-shifted mirror model.')

doc.add_heading('7. Limitations',1)
doc.add_paragraph('The mirror contributor frequency and per-gene magnitude are partially confounded. The asymmetry frequency and shift magnitude are likewise weakly identified, especially for BMI. The burden translation uses an effective dosage-responsive gene density rather than a directly observed genome-wide gene count, and the BMI native-unit conversion assumes 1 BMI SD = 4.7 kg/m2. Macé provides a quantitative deletion-burden coefficient but not an equally precise duplication coefficient; Auwerx is therefore used mainly as a directional burden check. These limitations affect exact parameter values more than the qualitative comparison of model directions.')

doc.add_heading('References',1)
refs=[
 'Sacks et al. Combinatorial effects of gene dosage, polygenic background and environment on complex traits. medRxiv (2026), doi:10.64898/2026.04.30.26352063.',
 'Kosmicki et al. Rare protein-coding variation and the genetic architecture of height in >1.4 million individuals. medRxiv (2026), doi:10.64898/2026.06.22.26355163.',
 'Akbari P et al. Sequencing of 640,000 exomes identifies GPR75 variants associated with protection from obesity. Science 373:eabf8683 (2021), doi:10.1126/science.abf8683.',
 'Backman JD et al. Exome sequencing and analysis of 454,787 UK Biobank participants. Nature 599:628-634 (2021), doi:10.1038/s41586-021-04103-z.',
 'Macé A et al. CNV-association meta-analysis in 191,161 European adults reveals new loci associated with anthropometric traits. Nature Communications 8:744 (2017), doi:10.1038/s41467-017-00556-x.',
 'Auwerx C et al. The individual and global impact of copy-number variants on complex human traits. American Journal of Human Genetics 109:647-668 (2022), doi:10.1016/j.ajhg.2022.02.010.'
]
for r in refs:doc.add_paragraph(r,style=None)

docx_path=ROOT/'report/CNV_height_BMI_corrected_models_report.docx';doc.save(docx_path)

# README and sources
(ROOT/'README.md').write_text(textwrap.dedent(f'''\
# Corrected directional CNV models for height and BMI

This package reruns the full analysis after correcting the direction of asymmetric effects:

- **Height:** asymmetric/generalized-dysfunction effects move both copy states toward the **lower-left**.
- **BMI:** asymmetric/generalized-dysfunction effects move both copy states toward the **upper-right**.

For each trait, the package compares:

1. Mirror alone
2. Mirror + continuous directional asymmetry
3. Mirror + sparse single hit

Results fitted to the **Sacks et al. single-locus CNVs** are separated from predictions calibrated to the **Macé et al. deletion-burden coefficient** and checked against the global burden direction in **Auwerx et al.**

Primary report: `report/CNV_height_BMI_corrected_models_report.pdf`

Python reference implementation: `code/run_corrected_analysis.py`
R implementation: `code/run_corrected_analysis.R` and `code/corrected_model_functions.R`

Key corrected point estimates:

- Height continuous: q = {100*Hfits['continuous']['q']:.2f}%, shift = {Hfits['continuous']['shift_magnitude']:.3f} SD.
- Height single hit: q = {100*Hfits['single']['q']:.2f}%, shift = {Hfits['single']['shift_magnitude']:.3f} SD.
- BMI continuous: diffuse small upper-right load; q reaches the broad bound and should not be interpreted literally.
- BMI single hit: q = {100*Bfits['single']['q']:.2f}%, shift = {Bfits['single']['shift_magnitude']:.3f} SD.

The numerical outputs were generated with the Python reference implementation. R was unavailable in the artifact environment, so the R workflow was structurally checked but not executed.
'''))
(ROOT/'SOURCES.md').write_text('\n'.join(['# Sources','',*['- '+x for x in refs],'- Exact Sacks locus-level inputs are included under each trait input directory.']))

# Copy this running script as reference (skip when running in place, where it would copy onto itself)
_self=Path(__file__).resolve();_dst=(ROOT/'code/run_corrected_analysis.py').resolve()
if _self!=_dst: shutil.copy2(_self,_dst)
print(json.dumps({'root':str(ROOT),'height_base':height_base,'height_fits':Hfits,'bmi_base':bmi_base,'bmi_fits':Bfits,'height_comp':Hcomp.to_dict('records'),'bmi_comp':Bcomp.to_dict('records')},indent=2,default=float))
