"""Reproduce all tables and figures; no network or secret data required."""
import csv,json,math,random,time,sys,platform
from pathlib import Path
import sympy as sp
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sourcecrt.core import *
from sourcecrt.verify import verify

ROOT=Path(__file__).resolve().parents[1]
for folder in ['results','figures']: (ROOT/folder).mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Serif','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight'})


def multiplication_matrix(f,g):
    t=sp.Symbol('t');d=sp.degree(f,t)
    cols=[]
    for j in range(d):
        h=sp.Poly(sp.rem(g*t**j,f,t),t)
        cols.append([int(h.nth(i)) for i in range(d)])
    return sp.Matrix(d,d,lambda i,j:cols[j][i])


def char_mod(A,p):
    """Faddeev-LeVerrier over F_p, independently of SymPy charpoly."""
    n=A.rows;X=[[int(A[i,j])%p for j in range(n)] for i in range(n)]
    B=[[int(i==j) for j in range(n)] for i in range(n)]
    out=[]
    for k in range(1,n+1):
        C=[[sum(X[i][l]*B[l][j] for l in range(n))%p for j in range(n)] for i in range(n)]
        c=(-sum(C[i][i] for i in range(n))*pow(k,-1,p))%p
        out.append(c)
        for i in range(n):C[i][i]=(C[i][i]+c)%p
        B=C
    return out


def run():
    rng=random.Random(6312028);t=sp.Symbol('t');X=sp.Symbol('X')
    primes=list(sp.primerange(101,1000))[:48]
    rows=[];certs=[];objects=[];transcripts=[];start=time.perf_counter()
    # Fixed prereported factorial grid, not selected on outcome.
    for d in [4,8,12,16]:
      f=t**d-t-1
      assert sp.Poly(f,t).is_irreducible
      for instance in range(12):
        coeff=[rng.randint(-3,3) for _ in range(d)]
        g=sum(c*t**j for j,c in enumerate(coeff))
        A=multiplication_matrix(f,g)
        truth=[int(c) for c in A.charpoly(X).all_coeffs()[1:]]
        # Independent algebraic ground truth via resultant.
        norm=sp.Poly(sp.resultant(f,X-g,t),X)
        assert truth==[int(c) for c in norm.all_coeffs()[1:]]
        R=max(sum(abs(int(A[i,j])) for j in range(d)) for i in range(d))
        bounds=[(-math.comb(d,k)*R**k,math.comb(d,k)*R**k) for k in range(1,d+1)]
        assert all(l<=v<=u for v,(l,u) in zip(truth,bounds))
        objects.append(dict(degree=d,instance=instance,g_coefficients=coeff,multiplication_matrix=[[int(A[i,j]) for j in range(d)] for i in range(d)],norm_coefficients=truth,bounds=bounds))
        raw=[char_mod(A,p) for p in primes]
        assert all(row==[v%p for v in truth] for p,row in zip(primes,raw))
        for mode in ['dense','coherent','sparse']:
          for placement in ['round_robin','blocked']:
           for f_budget in [1,2]:
             groups=6
             owner=[i%groups if placement=='round_robin' else min(groups-1,i//8) for i in range(len(primes))]
             corrupt={(instance+h)%groups for h in range(f_budget)}
             received=[]
             # Coherent wrong input: use g+1, preserving monicity.
             wrong=[int(c) for c in (A+sp.eye(d)).charpoly(X).all_coeffs()[1:]]
             for i,(p,row) in enumerate(zip(primes,raw)):
               y=row.copy()
               if owner[i] in corrupt:
                 if mode=='dense':y=[(v+rng.randrange(1,p))%p for v in row]
                 elif mode=='coherent':y=[v%p for v in wrong]
                 else:y[-1]=(y[-1]+1)%p
               received.append(y)
             transcripts.append(dict(degree=d,instance=instance,mode=mode,placement=placement,fault_budget=f_budget,owners=owner,faulty_sources=sorted(corrupt),moduli=primes,residues=received))
             stop={'joint':None,'scalar':None,'universal':None,'clean':None}
             decoder_s=0.;verifier_s=0.;chosen=None
             for n in range(1,len(primes)+1):
               ms=primes[:n];ys=received[:n]
               Fs=grouped_supports([{i for i in range(n) if owner[i]==h} for h in range(groups)],f_budget)
               if stop['clean'] is None and math.prod(p for i,p in enumerate(ms) if owner[i] not in corrupt)>max(u-l for l,u in bounds):stop['clean']=n
               if stop['universal'] is None and separation_modulus(ms,Fs)>max(u-l for l,u in bounds):stop['universal']=n
               if stop['scalar'] is None and all(s=='unique' for s in scalar_status(ms,ys,bounds,Fs)):stop['scalar']=n
               if stop['joint'] is None:
                 ts=time.perf_counter();c=certificate(ms,ys,bounds,Fs);decoder_s+=time.perf_counter()-ts
                 if c['status']=='unique':
                   assert c['candidate']==truth
                   ts=time.perf_counter();assert verify(ms,ys,bounds,Fs,c);verifier_s=time.perf_counter()-ts
                   stop['joint']=n
                   chosen={'degree':d,'instance':instance,'mode':mode,'placement':placement,'fault_budget':f_budget,'moduli':ms,'residues':ys,'bounds':bounds,'supports':Fs,'certificate':c}
               if all(v is not None for v in stop.values()):break
             assert all(v is not None for v in stop.values()),stop
             rows.append(dict(degree=d,instance=instance,mode=mode,placement=placement,fault_budget=f_budget,
                              **stop,decoder_s=decoder_s,verifier_s=verifier_s,
                              bound_bits=max(u-l for l,u in bounds).bit_length()))
             certs.append(chosen)
    with (ROOT/'results/trials.csv').open('w') as out:
        w=csv.DictWriter(out,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (ROOT/'results/objects.json').write_text(json.dumps(objects,separators=(',',':')))
    (ROOT/'results/transcripts.json').write_text(json.dumps(transcripts,separators=(',',':')))
    (ROOT/'results/certificates.json').write_text(json.dumps(certs,separators=(',',':')))
    summary={'trials':len(rows),'objects':48,'seed':6312028,'modular_images':48*48,
             'joint_mean':np.mean([r['joint'] for r in rows]),'scalar_mean':np.mean([r['scalar'] for r in rows]),
             'universal_mean':np.mean([r['universal'] for r in rows]),
             'joint_beats_scalar':sum(r['joint']<r['scalar'] for r in rows),
             'joint_beats_universal':sum(r['joint']<r['universal'] for r in rows),
             'max_scalar_saving':max(r['scalar']-r['joint'] for r in rows),
             'seconds':time.perf_counter()-start,'python':platform.python_version(),'sympy':sp.__version__,
             'platform':platform.platform()}
    (ROOT/'results/summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
    # Figure 1: exact family of joint/scalar candidate boxes.
    fig,ax=plt.subplots(figsize=(6.3,3.4))
    ax.set_xlim(-1,12);ax.set_ylim(-1,12)
    ax.scatter([0,11,0,11],[0,0,7,7],c=['#13786b','#d16936','#d16936','#d16936'],s=70)
    ax.annotate('Jointly admissible: (0, 0)',(0,0),(1,2),arrowprops={'arrowstyle':'->'})
    ax.annotate('Scalar-only alternative',(11,0),(5,4),arrowprops={'arrowstyle':'->'})
    ax.annotate('Scalar-only alternative',(0,7),(2,9),arrowprops={'arrowstyle':'->'})
    ax.set_xlabel('First coefficient');ax.set_ylabel('Second coefficient')
    ax.set_xticks([0,7,11]);ax.set_yticks([0,7,11]);fig.savefig(ROOT/'figures/joint.pdf');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(7.1,3.2),sharey=True)
    for ax,placement in zip(axs,['round_robin','blocked']):
      for method,color,label in [('joint','#13786b','Joint exact certificate'),('scalar','#d16936','Scalar exact certificates'),('universal','#556b9c','Universal support bound')]:
        vals=[np.mean([r[method] for r in rows if r['degree']==d and r['placement']==placement]) for d in [4,8,12,16]]
        ax.plot([4,8,12,16],vals,'o-',color=color,label=label)
      ax.set_title(placement.replace('_',' ').title());ax.set_xlabel('Field degree');ax.set_xticks([4,8,12,16])
    axs[0].set_ylabel('Prime computations at stopping');axs[1].legend(fontsize=7,loc='upper left')
    fig.tight_layout();fig.savefig(ROOT/'figures/stopping.pdf');plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.3,3.2))
    for j,mode in enumerate(['dense','coherent','sparse']):
      vals=[r['scalar']-r['joint'] for r in rows if r['mode']==mode]
      ax.bar(j,np.mean(vals),color=['#13786b','#556b9c','#d16936'][j])
      ax.text(j,np.mean(vals)+.03,f'{sum(v>0 for v in vals)}/{len(vals)} improve',ha='center',fontsize=9)
    ax.set_xticks(range(3),['Dense random','Coherent wrong input','One coefficient'])
    ax.set_ylabel('Mean prime computations saved vs scalar');ax.set_ylim(0,max(1,ax.get_ylim()[1]*1.2))
    fig.tight_layout();fig.savefig(ROOT/'figures/modes.pdf');plt.close(fig)

if __name__=='__main__':run()
