from pathlib import Path
import csv, math, statistics, sys, json, os
from statistics import NormalDist

ROOT = Path(__file__).resolve().parent.parent
OUT = Path(os.environ.get('MR_OUTPUT_ROOT', ROOT / 'MR严格终检_20261012'))
RES = OUT / '新增核查数据'
RES.mkdir(parents=True, exist_ok=True)
SRC = Path(os.environ.get('MR_SOURCE_ROOT', ROOT / 'bmc_reproducibility'))
ND = NormalDist()
KEYS = ['cad','pad_finngen_r5','ais','laa','svs','ces']
LABELS = ['CAD','PAD','AIS','LAA','SVS','CES']
def read(p):
    with open(p,encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f))
def write(name, rows):
    with (RES/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def kept(name): return [r for r in read(SRC/'data_archive'/name) if r['mr_keep'].upper()=='TRUE']
def ivw(rows):
    x=[float(r['beta.exposure']) for r in rows]; y=[float(r['beta.outcome']) for r in rows]
    w=[1/float(r['se.outcome'])**2 for r in rows]
    den=sum(ww*xx*xx for ww,xx in zip(w,x))
    b=sum(ww*xx*yy for ww,xx,yy in zip(w,x,y))/den
    q=sum(ww*(yy-b*xx)**2 for ww,xx,yy in zip(w,x,y))
    se=math.sqrt(q/(len(rows)-1)/den)
    p=math.erfc(abs(b/se)/math.sqrt(2))
    return b,se,q,p

core={r['key']:r for r in read(SRC/'results/core_results_refitted.csv')}
diag={r['outcome_key']:r for r in read(SRC/'results/table_s3_univariable_diagnostics.csv')}
strength=[]; per_snp=[]; heterogeneity=[]; single=[]; chrom=[]; provenance=[]
reported_controls={'cad':123504,'pad_finngen_r5':206541,'ais':406111,'laa':297290,'svs':343560,'ces':355468}
cases={'cad':60801,'pad_finngen_r5':7098,'ais':34217,'laa':4373,'svs':5386,'ces':7193}
official_controls={'cad':123504,'pad_finngen_r5':206541,'ais':406111,'laa':406111,'svs':192662,'ces':406111}
official_n={'cad':184305,'pad_finngen_r5':None,'ais':440328,'laa':150765,'svs':198048,'ces':211763}
sets=[(k,kept(k+'_harmonised.csv')) for k in KEYS]
sets += [('smoking_to_bmi',kept('x_to_m_harmonised.csv')),('bmi_to_smoking',kept('reverse_bmi_to_lifetime_harmonised.csv'))]
for k,d in sets:
    fs=[(float(r['beta.exposure'])/float(r['se.exposure']))**2 for r in d]
    r2s=[f/(float(r['samplesize.exposure'])-2+f) for r,f in zip(d,fs)]
    n=statistics.mean(float(r['samplesize.exposure']) for r in d)
    r2=sum(r2s)
    fg=(n-len(d)-1)*r2/((1-r2)*len(d))
    strength.append(dict(analysis=k,nsnp=len(d),F_min=min(fs),F_mean=statistics.mean(fs),F_median=statistics.median(fs),F_max=max(fs),F_below_10=sum(f<10 for f in fs),sum_R2_bsen=r2,F_aggregate_approx=fg,N_mean=n))
    for r,f,r2 in zip(d,fs,r2s):
        per_snp.append(dict(analysis=k,SNP=r['SNP'],chromosome=r.get('chr',r.get('chr.exposure','')),beta_exposure=r['beta.exposure'],se_exposure=r['se.exposure'],N_exposure=r['samplesize.exposure'],F=f,R2_bsen=r2))
    if k not in KEYS: continue
    ref=core[k]; b,se,q,p=ivw(d)
    assert abs(b-float(ref['total_b']))<1e-12
    # Outcome heterogeneity I2; separate from exposure I2_GX used for Egger NOME.
    i2=max(0,(q-len(d)+1)/q)
    xx=[abs(float(r['beta.exposure'])) for r in d]
    wx=[1/float(r['se.exposure'])**2 for r in d]
    xm=sum(w*x for w,x in zip(wx,xx))/sum(wx)
    qx=sum(w*(x-xm)**2 for w,x in zip(wx,xx))
    i2gx=max(0,(qx-len(d)+1)/qx)
    heterogeneity.append(dict(analysis=k,nsnp=len(d),IVW_Q=q,df=len(d)-1,Q_P=diag[k]['Q_pval'],I2_outcome=i2,Q_GX=qx,I2_GX=i2gx,egger_intercept=diag[k]['egger_intercept'],egger_P=diag[k]['egger_intercept_p']))
    for j,r in enumerate(d):
        bb,ss,qq,pp=ivw(d[:j]+d[j+1:])
        single.append(dict(analysis=k,excluded_SNP=r['SNP'],n_retained=len(d)-1,b=bb,se=ss,lci=bb-1.96*ss,uci=bb+1.96*ss,p=pp,full_b=b))
    chromosomes=sorted(set(r['chr'] for r in d),key=lambda s:int(s))
    for ch in chromosomes:
        dd=[r for r in d if r['chr']!=ch]
        bb,ss,qq,pp=ivw(dd)
        chrom.append(dict(analysis=k,excluded_chromosome=ch,n_removed=len(d)-len(dd),n_retained=len(dd),b=bb,se=ss,lci=bb-1.96*ss,uci=bb+1.96*ss,p=pp,full_b=b))
    nn=sorted(set(r['samplesize.outcome'] for r in d))
    provenance.append(dict(analysis=k,accession=d[0]['id.outcome'],manuscript_cases=cases[k],manuscript_controls=reported_controls[k],manuscript_total=cases[k]+reported_controls[k],frozen_N_field=';'.join(nn),OpenGWAS_cases=cases[k],OpenGWAS_controls=official_controls[k],OpenGWAS_total_field=official_n[k],case_control_total=cases[k]+official_controls[k],status='一致' if k in ['cad','pad_finngen_r5','ais'] else ('正文对照数错误，官方与冻结N一致，已修正' if k=='svs' else '正文与官方不一致，官方N字段与病例对照合计亦不一致，待原始readme核定')))

write('S9_逐工具F.csv',per_snp);write('S10_工具强度汇总.csv',strength)
write('S11_异质性I2及I2GX.csv',heterogeneity);write('S12_逐SNP剔除完整结果.csv',single)
write('S13_逐染色体剔除.csv',chrom);write('S14_样本量来源逐项核对.csv',provenance)

# Minimum detectable OR from observed model precision, not post-hoc observed power.
# Normal approximation; 80% power; conservative six-test Bonferroni alpha 0.05/6.
z1=ND.inv_cdf(1-.05/2); z6=ND.inv_cdf(1-(.05/6)/2); zpow=ND.inv_cdf(.8)
mde=[]
for k in KEYS:
    for model in ['total','smoking','bmi']:
        r=core[k]; se=float(r[model+'_se'])
        mde.append(dict(analysis=k,model=model,se_logOR=se,OR=r[model+'_or'],CI_low=r[model+'_lci'],CI_high=r[model+'_uci'],p_raw=r[model+'_p'],p_FDR=r[model+'_fdr'],MDE_OR_80pct_alpha005=math.exp((z1+zpow)*se),MDE_OR_80pct_alpha005_div6=math.exp((z6+zpow)*se),method='fixed-SE normal approximation; two-sided; not observed-effect power'))
write('S15_最小可检测效应.csv',mde)

# Continuous traits only, independent-sample Fisher-z approximation.
steiger=[]
for k,name in [('smoking_to_bmi','x_to_m_harmonised.csv'),('bmi_to_smoking','reverse_bmi_to_lifetime_harmonised.csv')]:
    d=kept(name)
    ns=[statistics.mean(float(r['samplesize.'+s]) for r in d) for s in ['exposure','outcome']]
    rs=[]
    for s in ['exposure','outcome']:
        rs.append(sum((float(r['beta.'+s])/float(r['se.'+s]))**2/(float(r['samplesize.'+s])-2+(float(r['beta.'+s])/float(r['se.'+s]))**2) for r in d))
    z=(math.atanh(math.sqrt(rs[0]))-math.atanh(math.sqrt(rs[1])))/math.sqrt(1/(ns[0]-3)+1/(ns[1]-3))
    steiger.append(dict(direction=k,nsnp=len(d),R2_exposure=rs[0],R2_outcome=rs[1],mean_N_exposure=ns[0],mean_N_outcome=ns[1],correct_direction=rs[0]>rs[1],z=z,p=math.erfc(abs(z)/math.sqrt(2)),measurement_reliability_threshold_Rx_over_Ry=rs[0]/rs[1],assumptions='quantitative traits; independent samples; no measurement-error correction; no selection filtering'))
write('S16_连续性状Steiger近似.csv',steiger)
summary={'per_snp_F_min':min(r['F_min'] for r in strength),'weak_snps':sum(r['F_below_10'] for r in strength),'I2_GX':{r['analysis']:r['I2_GX'] for r in heterogeneity},'leave_chromosome_sign_reversals':sum(r['b']*r['full_b']<=0 for r in chrom),'steiger':steiger,'MDE_null_results':[r for r in mde if float(r['p_FDR'])>=.05],'verdict':'未通过用户严格标准；元数据和独立性验证仍有实质缺项'}
(RES/'终检数值摘要.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False,indent=2))
