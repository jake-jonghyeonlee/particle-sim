import json,shutil
from pathlib import Path
import numpy as np,pandas as pd,zarr,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from particle_sim.config import load_manifest
from particle_sim.pipeline import run_probability
from particle_sim.io.masks import GridSpec,MaskSet
from particle_sim.schema import OpenPath
from particle_sim.failure.open_failure import simulate_open,causes_open,blocked_mask
w=Path('docs/report/work'); out=Path('docs/report/figures')
m=load_manifest(w/'case/manifest.yaml')
run_probability(w/'case/manifest.yaml',Path('configs/scenarios/strict.yaml'),w/'zero',sobol_exponent_m=6,seed=42)
run_probability(w/'case/manifest.yaml',Path('configs/scenarios/nominal.yaml'),w/'repeat',sobol_exponent_m=6,seed=42)
a=zarr.open_group(w/'run/hotspots.zarr',mode='r'); b=zarr.open_group(w/'zero/probability.zarr',mode='r'); p=zarr.open_group(w/'prob/probability.zarr',mode='r'); q=zarr.open_group(w/'repeat/probability.zarr',mode='r')
result={'zero_max_pixel_difference':float(np.abs(a['nominal_failure'][:]-b['probability'][:]).max()),'zero_union_max_pixel_difference':float(np.abs(a['union_failure'][:]-b['union_probability'][:]).max()),'repro_probability_bitwise_equal':p['probability'][:].tobytes()==q['probability'][:].tobytes(),'repro_union_bitwise_equal':p['union_probability'][:].tobytes()==q['union_probability'][:].tobytes()}
checks=[]; fig,axs=plt.subplots(2,4,figsize=(12,6))
for j,name in enumerate(['straight','detour','branch','corner']):
 wire=np.zeros((24,24),bool); wire[11:14,3:21]=1; terms=((12,3),(12,20))
 if name=='detour': wire[5:8,6:18]=1;wire[5:14,6:9]=1;wire[5:14,15:18]=1
 if name=='branch': wire[3:12,11:14]=1
 if name=='corner':
  wire[:]=0;wire[3:15,6:9]=1;wire[12:15,6:21]=1;terms=((3,7),(13,20))
 for i,(py,px) in enumerate([(0.5,0.5),(1.0,0.5)]):
  grid=GridSpec(shape_yx=wire.shape,pixel_size_y_um=py,pixel_size_x_um=px)
  masks=MaskSet(grid=grid,layers={'M2':wire},electrodes={'W':wire},input_hashes={})
  path=OpenPath(path_id='W',electrode='W',terminal_a_yx_px=terms[0],terminal_b_yx_px=terms[1])
  pred=simulate_open(m,masks,path,4.0).failure
  yy,xx=np.mgrid[:24,:24];brute=np.zeros(wire.shape,bool)
  for y,x in np.ndindex(wire.shape):brute[y,x]=causes_open(wire,blocked_mask(np.hypot((yy-y)*py,(xx-x)*px),2,0),*terms)
  diff=int(np.count_nonzero(pred!=brute));assert diff==0
  checks.append(dict(shape=name,pixel_y_um=py,pixel_x_um=px,diameter_um=4,different_pixels=diff,hotspot_pixels=int(pred.sum())))
  rgb=np.ones((24,24,3));rgb[wire]=[.28,.36,.43];rgb[pred]=[1,.39,.16]
  ax=axs[i,j];ax.imshow(rgb,origin='upper',interpolation='nearest');ax.scatter([t[1] for t in terms],[t[0] for t in terms],s=22,c='#16b9b0');ax.set_title(f'{name} | dy={py}, dx={px} um');ax.set_xticks([]);ax.set_yticks([])
fig.tight_layout();fig.savefig(out/'open_checks.png',dpi=180);plt.close(fig)
result['open_full_pixel_checks']=checks
for label,dirname in [('nominal','run'),('probability','prob'),('triple','triple')]:
 for f in (w/dirname/'report').glob('*.png'):shutil.copy2(f,out/f'{label}_{f.name}')
for dirname in ['run','prob','large-short']:
 path=w/dirname
 if path.exists():result[dirname+'_output_bytes']=sum(f.stat().st_size for f in path.rglob('*') if f.is_file() and 'report' not in f.relative_to(path).parts)
result['legacy_totals']=pd.read_csv(w/'run/legacy_vs_edt.csv').select_dtypes(include='number').sum().to_dict()
(w/'recheck.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
