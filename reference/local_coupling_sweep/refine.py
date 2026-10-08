from pathlib import Path
import json
from run import simulate,configs
P=Path(__file__).resolve().parent
out=[]
for index in [2,6,11]:
 for k,dt in enumerate([.05,.025]):
  out.append(simulate(configs[index],1000+index*10+k,trials=4000,dt=dt,suffix=f'_dt_{dt}'))
  (P/'refinement.json').write_text(json.dumps(out,indent=2))
