from particle_sim.validation.synthetic_cases import build_case,write_case
import numpy as np
from pathlib import Path
case=build_case('pilot')
case.manifest['product_id']='SYN_BENCH_1000'
case.manifest['coordinate_system'].update(pixel_size_x_um=0.25,pixel_size_y_um=0.25,roi_width_um=250.0,roi_height_um=250.0)
case.manifest['particle_sizes_um']=list(range(1,31))+[40,50]
for mapping in (case.layer_masks,case.electrode_masks):
 for key,arr in mapping.items(): mapping[key]=np.repeat(np.repeat(arr,5,axis=0),5,axis=1)
write_case(case,Path('docs/report/work/large-case'))
