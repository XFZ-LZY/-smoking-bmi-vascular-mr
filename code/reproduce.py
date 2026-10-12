from pathlib import Path
import os,sys,subprocess,shutil,csv,math
root=Path(__file__).resolve().parent.parent
rs=sys.argv[1] if len(sys.argv)>1 else shutil.which('Rscript')
if not rs:raise SystemExit('Supply the Rscript executable path')
src=root/'bmc_reproducibility';out=root/'reproduced_audit';out.mkdir(exist_ok=True)
subprocess.run([rs,'code/final_article_analysis.R'],cwd=src,check=True)
env=os.environ.copy();env['MR_SOURCE_ROOT']=str(src);env['MR_OUTPUT_ROOT']=str(out)
subprocess.run([sys.executable,str(root/'code/strict_mr_audit.py')],env=env,cwd=root,check=True)
# strict_mr_audit's legacy S14 is intentionally superseded by audit/S14_corrected_stroke_metadata.csv.
legacy=out/'新增核查数据/S14_样本量来源逐项核对.csv'
if legacy.exists():legacy.rename(legacy.with_name('historical_S14_do_not_use.csv'))
shutil.copy2(root/'audit/instrument_flags.csv',out/'instrument_flags.csv')
subprocess.run([rs,'code/independence_sensitivity.R','bmc_reproducibility','reproduced_audit'],cwd=root,check=True)
def read(p):return list(csv.DictReader(open(p,encoding='utf-8-sig')))
new=read(out/'independence_sensitivity_results.csv');ref=read(root/'audit/independence_sensitivity_results.csv')
assert len(new)==len(ref)==60
for a,b in zip(new,ref):
 for k in ['scheme','key','model','exposure','n','removed']:assert a[k]==b[k],(k,a,b)
 for k in ['b','se','p','fdr','conditional_F']:
  if a[k]=='NA':assert b[k]=='NA'
  else:assert abs(float(a[k])-float(b[k]))<1e-9,(k,a,b)
print('PASS: original primary estimates and 60 instrument-sensitivity rows reproduced.')
