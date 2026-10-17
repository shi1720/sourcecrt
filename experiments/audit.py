"""Certificate tampering, layout theorem, and out-of-model counterexample."""
import copy,itertools,json,random,math
from pathlib import Path
from sourcecrt.core import *
from sourcecrt.verify import verify

root=Path(__file__).resolve().parents[1]
certs=json.loads((root/'results/certificates.json').read_text())
checked=mutated=0
for case in certs:
    args=[case[k] for k in ['moduli','residues','bounds','supports']]
    c=case['certificate']
    assert verify(*args,c);checked+=1
    z=copy.deepcopy(c);z['candidate'][0]+=1
    assert not verify(*args,z);mutated+=1
    z=copy.deepcopy(c);z['records'].pop()
    assert not verify(*args,z);mutated+=1
    z=copy.deepcopy(c)
    r=z['records'][0]
    r['exclude_coordinate']=len(case['bounds'])
    assert not verify(*args,z);mutated+=1
rng=random.Random(6312029);layout_checks=0
for n in [5,10,20,40]:
 for h in [3,5,8]:
  for f in range(1,(h-1)//2+1):
   for rep in range(50):
    # Pairwise coprime prime-power moduli, up to very unequal bit lengths.
    ps=[2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173][:n]
    ms=[p**rng.randint(1,8) for p in ps]
    groups=greedy_layout(ms,h);Fs=grouped_supports(groups,f)
    loads=sorted(math.prod(ms[i] for i in g) for g in groups)
    k=2*f;M=math.prod(ms);S=math.prod(loads[:h-k])
    assert S==separation_modulus(ms,Fs)
    assert S**h*max(ms)**(k*(h-k)) >= M**(h-k)
    assert max(loads)<=min(loads)*max(ms)
    layout_checks+=1
# An undeclared common-mode failure can produce a valid certificate for a wrong
# answer. This tests the trust boundary, not a defect in the theorem.
ms=[5,7,11];rs=[[1],[1],[1]];bounds=[(0,1)];Fs=[(0,),(1,),(2,)]
c=certificate(ms,rs,bounds,Fs)
assert c['candidate']==[1] and verify(ms,rs,bounds,Fs,c)
# If physical truth were zero, all three sources would have failed, violating f=1.
summary=dict(certificates_checked=checked,tampered_certificates_rejected=mutated,
             exact_layout_checks=layout_checks,out_of_model_counterexample=True)
(root/'results/audit.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
