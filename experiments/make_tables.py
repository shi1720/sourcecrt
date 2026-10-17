import csv,statistics
from pathlib import Path
root=Path(__file__).resolve().parents[1]
r=list(csv.DictReader((root/'results/trials.csv').open()))
lines=[]
for mode,label in [('dense','Dense random'),('coherent','Coherent wrong input'),('sparse','Single coefficient')]:
 q=[x for x in r if x['mode']==mode]
 vals=[statistics.mean(int(x[k]) for x in q) for k in ['joint','scalar','universal']]
 wins=sum(int(x['joint'])<int(x['scalar']) for x in q)
 lines.append(label+' & 192 & '+' & '.join(f'{v:.3f}' for v in vals)+f' & {wins}'+r' \\')
(root/'paper/mode_rows.tex').write_text(r'\begin{tabular}{lrrrrr}\toprule'+'\n'+r'Mechanism & Cases & Joint & Scalar & Universal & Joint $<$ scalar\\\midrule'+'\n'+'\n'.join(lines)+'\n'+r'\bottomrule\end{tabular}'+'\n')
lines=[]
for place,label in [('round_robin','Round robin'),('blocked','Blocked')]:
 for f in ['1','2']:
  q=[x for x in r if x['placement']==place and x['fault_budget']==f]
  vals=[statistics.mean(int(x[k]) for x in q) for k in ['joint','scalar','universal']]
  lines.append(label+' & '+f+' & 144 & '+' & '.join(f'{v:.3f}' for v in vals)+r' \\')
(root/'paper/layout_rows.tex').write_text(r'\begin{tabular}{lrrrrr}\toprule'+'\n'+r'Placement & $f$ & Cases & Joint & Scalar & Universal\\\midrule'+'\n'+'\n'.join(lines)+'\n'+r'\bottomrule\end{tabular}'+'\n')
