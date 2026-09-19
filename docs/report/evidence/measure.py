import json, subprocess, sys, time, resource
from pathlib import Path
name=sys.argv[1]
start=time.perf_counter()
p=subprocess.run(sys.argv[2:],capture_output=True,text=True)
elapsed=time.perf_counter()-start
usage=resource.getrusage(resource.RUSAGE_CHILDREN)
result=dict(command=sys.argv[2:],returncode=p.returncode,wall_seconds=elapsed,peak_rss_bytes=usage.ru_maxrss,stdout=p.stdout,stderr=p.stderr)
Path('docs/report/work/'+name+'.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
print(json.dumps({k:v for k,v in result.items() if k not in ('stdout','stderr')}))
sys.exit(p.returncode)
