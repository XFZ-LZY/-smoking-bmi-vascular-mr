from pathlib import Path
import sys, csv, math, json, os
ROOT=Path(__file__).resolve().parent.parent

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager, ticker
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.lines import Line2D

OUT=Path(os.environ.get('MR_OUTPUT_ROOT',ROOT/'MR严格终检_20261012')); FIG=OUT/'figures'; FIG.mkdir(exist_ok=True)
SRC=Path(os.environ.get('MR_SOURCE_ROOT',ROOT/'bmc_reproducibility')); RES=SRC/'results'; DATA=SRC/'data_archive'; NEW=OUT/'新增核查数据'
for p in Path('C:/Windows/Fonts').glob('arial*.ttf'): font_manager.fontManager.addfont(str(p))
BLUE='#0F4D92'; TEAL='#42949E'; RED='#B64342'; GREEN='#568A5A'; GRAY='#767676'; LIGHT='#E7E9EC'; INK='#272727'
# scientific-figure-making: project-specific implementation of the documented API.
def apply_publication_style():
    plt.rcParams.update({'font.family':['Arial','DejaVu Sans'],'font.size':10,'axes.linewidth':1.05,'axes.spines.right':False,'axes.spines.top':False,'legend.frameon':False,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white','axes.unicode_minus':True})
apply_publication_style()
def read(p):
    with open(p,encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def n(v):return float(v)
def pval(v):
    v=float(v)
    return f'{v:.2e}' if v<.001 else f'{v:.4f}'.rstrip('0').rstrip('.')
def finalize_figure(fig,stem):
    fig.canvas.draw()
    for ext in ['png','pdf','svg']:
        fig.savefig(FIG/f'{stem}.{ext}',dpi=600 if ext=='png' else 300,bbox_inches='tight',pad_inches=.10)
    plt.close(fig)
K=['cad','pad_finngen_r5','ais','laa','svs','ces']
names=['Coronary artery disease','Peripheral artery disease','Any ischaemic stroke','Large-artery stroke','Small-vessel stroke','Cardioembolic stroke']
short=['CAD','PAD','AIS','LAA','SVS','CES']; LAB=dict(zip(K,names)); SH=dict(zip(K,short))
core={r['key']:r for r in read(RES/'core_results_refitted.csv')}
methodmap={'Inverse variance weighted (multiplicative random effects)':'IVW','Weighted median':'Weighted median','Weighted mode':'Weighted mode','MR Egger':'MR-Egger'}
methods=['IVW','Weighted median','Weighted mode','MR-Egger']; cols=[BLUE,TEAL,GREEN,RED]; markers=['s','o','D','^']
allm={(r['outcome_key'],methodmap[r['method']]):r for r in read(RES/'table_s2_total_effects_all_estimators.csv')}
directions=read(RES/'table_s5_bidirectional_smoking_bmi.csv')
diag={r['outcome_key']:r for r in read(RES/'table_s3_univariable_diagnostics.csv')}
payload=[]
def error(ax,b,lo,hi,y,col=BLUE,marker='s',filled=True,size=5.6):
    ax.errorbar(b,y,xerr=[[b-lo],[hi-b]],fmt=marker,ms=size,mec=col,mew=1.1,mfc=col if filled else 'white',ecolor=col,elinewidth=1.45,capsize=2.5,zorder=5)
def forest_axis(ax,ys,log=True):
    ax.axvline(1 if log else 0,color=GRAY,lw=.9,ls=(0,(3,3)),zorder=2)
    ax.set_ylim(min(ys)-.6,max(ys)+.6);ax.set_yticks([])
    ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
    ax.tick_params(axis='x',labelsize=9,length=3)
    if log:ax.set_xscale('log');ax.xaxis.set_minor_locator(ticker.NullLocator())
def strip_axis(ax,ys):
    ax.set_ylim(min(ys)-.6,max(ys)+.6);ax.set_axis_off()
def panel_title(fig,x,y,text):fig.text(x,y,text,fontsize=12.0,weight='bold',va='top',color=INK)
def shade(ax,ys):
    for idx,y in enumerate(ys):
        if idx%2==0:ax.axhspan(y-.48,y+.48,facecolor='#F3F5F7',zorder=0)

# Figure 1: three stages, aligned columns, modest graphic hierarchy.
fig=plt.figure(figsize=(8.8,4.0));ax=fig.add_axes([0,0,1,1]);ax.set_xlim(0,8.8);ax.set_ylim(0,4);ax.axis('off')
def card(x,y,w,h,title,body,color):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.015,rounding_size=.035',facecolor='white',edgecolor=LIGHT,lw=1))
    ax.add_patch(Rectangle((x,y+h-.055),w,.055,facecolor=color,edgecolor='none'))
    ax.text(x+.15,y+h-.22,title,fontsize=11,weight='bold',color=color,va='top')
    ax.text(x+.15,y+h-.63,body,fontsize=9.3,va='top',linespacing=1.45,color=INK)
for x,label in [(.10,'1  Genetic instruments'),(3.12,'2  Estimation'),(6.15,'3  Vascular outcomes')]:ax.text(x,3.88,label,fontsize=11.3,weight='bold',va='top')
card(.1,2.07,2.35,1.28,'Lifetime smoking','UK Biobank · n = 462,690\n126 genome-wide variants',BLUE)
card(.1,.53,2.35,1.27,'Body mass index','GIANT · n = 322,154\nP < 5 × 10⁻⁸ · LD r² < 0.001',TEAL)
card(3.12,.53,2.35,2.82,'MR analysis','Total smoking associations\n\nMutually adjusted smoking\nand BMI associations\n\nBidirectional smoking–BMI',BLUE)
card(6.15,.53,2.48,2.82,'Six outcomes','Coronary artery disease\nPeripheral artery disease\nAny ischaemic stroke\nLarge-artery stroke\nSmall-vessel stroke\nCardioembolic stroke',INK)
for yy in [2.71,1.17]:ax.annotate('',xy=(3.07,2 if yy>2 else 1.5),xytext=(2.49,yy),arrowprops=dict(arrowstyle='-|>',color=GRAY,lw=1.3))
ax.annotate('',xy=(6.10,1.94),xytext=(5.53,1.94),arrowprops=dict(arrowstyle='-|>',color=GRAY,lw=1.3))
ax.text(.10,.12,'Sensitivity: median / mode / Egger / MR-PRESSO / SNP and chromosome deletion',fontsize=8.8,color=GRAY)
finalize_figure(fig,'Fig1')

# Generic full-width forest with dedicated label, estimate and P columns.
def forest_table(fig,rect,rows,title,xlabel,xlim,ticks,log=True,headers=('Outcome','OR (95% CI)','P','FDR P'),col=BLUE):
    left,bottom,width,height=rect
    showq=bool(headers[-1]); ratios=[2.60,2.50,2.65,1.10] if showq else [2.60,2.90,3.25]
    gs=fig.add_gridspec(1,len(ratios),left=left,right=left+width,bottom=bottom,top=bottom+height,width_ratios=ratios,wspace=.07)
    aa=[fig.add_subplot(gs[0,i]) for i in range(len(ratios))]; ys=list(range(len(rows)-1,-1,-1))
    for a in aa:shade(a,ys)
    for a in [a for i,a in enumerate(aa) if i!=1]:strip_axis(a,ys)
    forest_axis(aa[1],ys,log);aa[1].set_xlim(*xlim);aa[1].set_xticks(ticks,labels=[f'{t:g}' for t in ticks]);aa[1].set_xlabel(xlabel,fontsize=9,labelpad=7)
    for y,r in zip(ys,rows):
        label,b,lo,hi,p,q,color,mark=r
        aa[0].text(.02,y,label,va='center',fontsize=9.3)
        error(aa[1],b,lo,hi,y,color,mark,q is None or q<.05)
        digits=2 if log else 3
        aa[2].text(.02,y,f'{b:.{digits}f} ({lo:.{digits}f}–{hi:.{digits}f})',va='center',fontsize=9.5)
        aa[2].text(.98,y,pval(p),ha='right',va='center',fontsize=9.1,color=INK)
        if showq:aa[3].text(.97,y,pval(q) if q is not None else '—',ha='right',va='center',fontsize=9.1,color=color if q is not None and q<.05 else GRAY)
        payload.append({'figure':title,'label':label,'estimate':b,'lci':lo,'uci':hi,'p':p,'fdr_p':q})
    aa[0].text(.02,1.045,headers[0],transform=aa[0].transAxes,fontsize=9.7,weight='bold')
    aa[2].text(.02,1.045,headers[1],transform=aa[2].transAxes,fontsize=9.7,weight='bold')
    aa[2].text(.98,1.045,headers[2],ha='right',transform=aa[2].transAxes,fontsize=9.7,weight='bold')
    if showq:aa[3].text(.97,1.045,headers[3],ha='right',transform=aa[3].transAxes,fontsize=9.7,weight='bold')
    panel_title(fig,left,bottom+height+.115,title)
    return aa

rows=[]
for k in K:
    r=core[k];rows.append((LAB[k],n(r['total_or']),n(r['total_lci']),n(r['total_uci']),n(r['total_p']),n(r['total_fdr']),BLUE,'s'))
fig=plt.figure(figsize=(8.8,4.15))
forest_table(fig,(.015,.22,.965,.56),rows,'Lifetime smoking and vascular outcomes','OR per 1-SD higher smoking',(.75,3.6),[.8,1,1.5,2,3])
fig.text(.03,.02,'Random-effects IVW. Filled markers: FDR P < 0.05 within six outcomes. Lines: 95% CI.',fontsize=9,color=GRAY)
finalize_figure(fig,'Fig2')

def direction_fig(stem):
    fig=plt.figure(figsize=(8.8,6.0))
    for idx,direct in enumerate(['BMI to lifetime smoking','Lifetime smoking to BMI']):
        rr=[]
        for m,c,mark in zip(methods,cols,markers):
            r=next(r for r in directions if r['direction']==direct and methodmap[r['method']]==m)
            rr.append((m,n(r['b']),n(r['lci']),n(r['uci']),n(r['p']),None,c,mark))
        forest_table(fig,(.02,.60 if idx==0 else .13,.96,.24),rr,('a  ' if idx==0 else 'b  ')+direct,'β (SD outcome / SD exposure)',(-.025,.27) if idx==0 else (-.62,.74),[0,.1,.2] if idx==0 else [-.4,0,.4],False,headers=('Estimator','β (95% CI)','P',''))
    fig.text(.035,.015,'Each direction uses its own instrument set. Exact estimates and raw P values are shown.',fontsize=9,color=GRAY)
    finalize_figure(fig,stem)
direction_fig('Fig3');direction_fig('FigS4')

fig=plt.figure(figsize=(8.8,8.1))
rows=[]
for k in K:
    r=core[k]
    for model,c,mark,label in [('total',GRAY,'o','Total'),('smoking',BLUE,'s','BMI-adjusted')]:
        rows.append((f'{SH[k]} · {label}',n(r[model+'_or']),n(r[model+'_lci']),n(r[model+'_uci']),n(r[model+'_p']),n(r[model+'_fdr']),c,mark))
forest_table(fig,(.02,.57,.96,.31),rows,'a  Smoking: total and BMI-adjusted associations','OR per 1-SD higher smoking',(.75,3.6),[.8,1,1.5,2,3])
rows=[]
for k in K:
    r=core[k];rows.append((LAB[k],n(r['bmi_or']),n(r['bmi_lci']),n(r['bmi_uci']),n(r['bmi_p']),n(r['bmi_fdr']),TEAL,'o'))
forest_table(fig,(.02,.16,.96,.235),rows,'b  BMI associations adjusted for lifetime smoking','OR per 1-SD higher BMI',(.75,1.8),[.8,1,1.2,1.5,1.8])
fig.text(.03,.035,'Conditional F: smoking 26.6–27.6; BMI 23.5–24.2. Filled markers: within-family FDR P < 0.05.',fontsize=9,color=GRAY)
finalize_figure(fig,'Fig4')

# Supplementary Fig S1: sensitivity by outcome, manageable four rows per panel.
fig,axes=plt.subplots(2,3,figsize=(8.8,6.3));fig.subplots_adjust(left=.12,right=.985,bottom=.13,top=.88,wspace=.51,hspace=.65)
for ax,k in zip(axes.flat,K):
    for yy,m,c,mark in zip(range(3,-1,-1),methods,cols,markers):
        r=allm[k,m];b,se=n(r['b']),n(r['se']);error(ax,b,b-1.96*se,b+1.96*se,yy,c,mark,m=='IVW')
    forest_axis(ax,list(range(4)),False);ax.set_yticks(range(3,-1,-1),['IVW','Median','Mode','Egger'],fontsize=8.4)
    ax.set_title(SH[k],loc='left',weight='bold',pad=10);ax.set_xlabel('Log OR per SD',fontsize=8.5)
fig.suptitle('Total smoking associations across estimators',x=.05,ha='left',weight='bold',fontsize=12)
fig.text(.06,.025,'Lines: 95% CI. CAD MR-Egger has the opposite point-estimate direction; full results are retained.',fontsize=8.8,color=GRAY)
finalize_figure(fig,'FigS1')

# S2: three panels with text in independent columns, no text over data.
fig,axes=plt.subplots(1,3,figsize=(8.8,3.9),gridspec_kw={'width_ratios':[1.35,1.35,.55],'wspace':.42})
fig.subplots_adjust(left=.075,right=.98,bottom=.22,top=.79)
ys=list(range(5,-1,-1))
for y,k in zip(ys,K):
    r=diag[k];q,df=n(r['Q']),n(r['Q_df']);c=RED if n(r['Q_pval'])<.05 else BLUE
    axes[0].plot([df,q],[y,y],color=LIGHT,lw=2.5);axes[0].scatter(df,y,marker='D',s=25,facecolors='white',edgecolors=GRAY,zorder=4);axes[0].scatter(q,y,s=35,color=c,zorder=5)
    b,se=n(r['egger_intercept']),n(r['egger_intercept_se']);error(axes[1],b,b-1.96*se,b+1.96*se,y,RED if n(r['egger_intercept_p'])<.05 else BLUE,'o')
    axes[2].text(0,y,pval(r['Q_pval']),fontsize=8.7,va='center');axes[2].text(1,y,pval(r['egger_intercept_p']),fontsize=8.7,ha='right',va='center')
for ax in axes:ax.set_ylim(-.6,5.6)
axes[0].set_yticks(ys,short);axes[1].set_yticks([]);axes[0].set_xlabel('IVW Cochran Q');axes[1].set_xlabel('Egger intercept (95% CI)');axes[1].axvline(0,ls='--',lw=.8,color=GRAY)
axes[2].set_axis_off();axes[2].text(0,1.06,'Q P',transform=axes[2].transAxes,weight='bold',fontsize=8.8);axes[2].text(1,1.06,'Egger P',ha='right',transform=axes[2].transAxes,weight='bold',fontsize=8.8)
fig.text(.055,.94,'a  Heterogeneity',weight='bold',fontsize=11);fig.text(.44,.94,'b  Directional pleiotropy',weight='bold',fontsize=11)
fig.text(.07,.035,'Diamond: Q degrees of freedom. Red: diagnostic P < 0.05. Blue: P ≥ 0.05.',fontsize=9,color=GRAY)
finalize_figure(fig,'FigS2')

# S3: MVMR BMI sensitivity from the frozen fits.
sens=read(NEW/'S4_exact_normal_CI.csv'); sm=['IVW','Robust IVW','MVMR-Egger','Multivariable median'];sl=['IVW','Robust','Egger','Median']
fig,axes=plt.subplots(2,3,figsize=(8.8,6.3));fig.subplots_adjust(left=.12,right=.985,bottom=.13,top=.88,wspace=.49,hspace=.62)
for ax,k in zip(axes.flat,K):
    for yy,m,c,mark in zip(range(3,-1,-1),sm,[TEAL,BLUE,RED,GRAY],markers):
        r=next(r for r in sens if r['key']==k and r['method']==m);error(ax,n(r['or']),n(r['lci']),n(r['uci']),yy,c,mark,m=='IVW')
    forest_axis(ax,list(range(4)),True);ax.set_yticks(range(3,-1,-1),sl,fontsize=8.4);ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    lim=ax.get_xlim();ticks=[v for v in [.6,.8,1,1.2,1.4,1.6,1.8,2,2.2] if lim[0]<v<lim[1]];ax.set_xticks(ticks)
    ax.set_title(SH[k],loc='left',weight='bold');ax.set_xlabel('Conditional BMI OR',fontsize=8.5)
fig.suptitle('BMI associations conditional on smoking',x=.05,ha='left',weight='bold',fontsize=12)
fig.text(.06,.02,'Four MVMR estimators; 95% CI. Robust / Egger / median fits retain the archived analysis outputs.',fontsize=8.8,color=GRAY)
finalize_figure(fig,'FigS3')

def influence_grid(stem,csvname,field,title):
    rows=read(NEW/csvname); fig,axes=plt.subplots(2,3,figsize=(8.8,7.9));fig.subplots_adjust(left=.08,right=.98,bottom=.085,top=.91,wspace=.40,hspace=.38)
    for ax,k in zip(axes.flat,K):
        rr=[r for r in rows if r['analysis']==k];rr=sorted(rr,key=lambda r:n(r['b'])) if field=='excluded_SNP' else sorted(rr,key=lambda r:int(r[field]))
        yy=np.arange(len(rr));b=np.array([n(r['b']) for r in rr]);lo=np.array([n(r['lci']) for r in rr]);hi=np.array([n(r['uci']) for r in rr]);full=n(rr[0]['full_b'])
        ax.hlines(yy,lo,hi,color=BLUE,alpha=.3,lw=.5 if field=='excluded_SNP' else 1);ax.scatter(b,yy,s=5 if field=='excluded_SNP' else 10,color=BLUE,zorder=3)
        ax.axvline(full,color=RED,ls='--',lw=1);ax.axvline(0,color=GRAY,ls=':',lw=.8)
        ax.set_title(f'{SH[k]} · {len(rr)} deletions',loc='left',weight='bold',fontsize=10);ax.set_xlabel('Log OR per SD',fontsize=8.5)
        if field=='excluded_SNP':ax.set_yticks([]);ax.set_ylabel('Deleted SNP (ordered estimate)',fontsize=8)
        else:ax.set_yticks(yy,[r[field] for r in rr],fontsize=6.8);ax.set_ylabel('Deleted chromosome',fontsize=8)
    fig.suptitle(title,x=.05,ha='left',weight='bold',fontsize=12)
    fig.text(.07,.017,'Blue: re-estimated IVW with 95% CI. Red dashed: full-data IVW. Grey dotted: null.',fontsize=9,color=GRAY)
    finalize_figure(fig,stem)
influence_grid('FigS5','S12_逐SNP剔除完整结果.csv','excluded_SNP','Leave-one-SNP-out sensitivity across six outcomes')
influence_grid('FigS9','S13_逐染色体剔除.csv','excluded_chromosome','Leave-one-chromosome-out sensitivity')

ind=read(RES/'table_s8_exploratory_indirect_effects.csv');fig=plt.figure(figsize=(8.8,4.15));rr=[]
for r in ind:rr.append((LAB[r['key']],n(r['indirect_log_or']),n(r['lci']),n(r['uci']),n(r['p']),n(r['fdr']),BLUE,'s'))
forest_table(fig,(.02,.22,.96,.56),rr,'Exploratory smoking × BMI path products','Product estimate (log-odds)',(-.06,.16),[-.04,0,.04,.08,.12],False,headers=('Outcome','β (95% CI)','P','FDR P'))
fig.text(.035,.02,'A × B, zero-covariance delta approximation. These products are not confirmed mediation effects.',fontsize=8.8,color=GRAY)
finalize_figure(fig,'FigS6')

# S7 scatter: consistently orient every SNP toward greater smoking for both axes.
fig,axes=plt.subplots(2,3,figsize=(8.8,6.7));fig.subplots_adjust(left=.095,right=.99,bottom=.15,top=.85,wspace=.37,hspace=.43)
for ax,k in zip(axes.flat,K):
    d=[r for r in read(DATA/(k+'_harmonised.csv')) if r['mr_keep'].upper()=='TRUE']
    x=np.array([n(r['beta.exposure']) for r in d]);y=np.array([n(r['beta.outcome']) for r in d]);flip=np.where(x<0,-1,1);x=x*flip;y=y*flip
    ax.scatter(x,y,s=14,edgecolor=GRAY,facecolor='white',linewidth=.55,alpha=.8)
    xx=np.array([0,x.max()*1.04])
    for m,c in zip(methods,cols):ax.plot(xx,n(allm[k,m]['b'])*xx+(n(diag[k]['egger_intercept']) if m=='MR-Egger' else 0),color=c,lw=1.4)
    ax.axhline(0,color=LIGHT,lw=.8);ax.set_title(f'{SH[k]} · n = {len(d)}',loc='left',weight='bold',fontsize=10)
    ax.tick_params(labelsize=8)
fig.suptitle('SNP associations and total-effect slopes',x=.05,ha='left',weight='bold',fontsize=12)
fig.legend([Line2D([0],[0],color=c,lw=1.8) for c in cols],methods,loc='upper center',bbox_to_anchor=(.52,.94),ncol=4,fontsize=9)
fig.supxlabel('SNP–smoking association (SD)',y=.065,fontsize=10);fig.supylabel('SNP–outcome association (log-odds)',x=.01,fontsize=10)
finalize_figure(fig,'FigS7')

# S8 funnel: inverse SE of Wald ratio; not a test that pleiotropy is absent.
fig,axes=plt.subplots(2,3,figsize=(8.8,6.1));fig.subplots_adjust(left=.095,right=.99,bottom=.14,top=.89,wspace=.35,hspace=.40)
for ax,k in zip(axes.flat,K):
    d=[r for r in read(DATA/(k+'_harmonised.csv')) if r['mr_keep'].upper()=='TRUE']
    x=np.array([n(r['beta.exposure']) for r in d]);y=np.array([n(r['beta.outcome']) for r in d]);s=np.array([n(r['se.outcome']) for r in d])/np.abs(x)
    ratios=y/x;prec=1/s;ax.scatter(ratios,prec,s=14,facecolor='white',edgecolor=GRAY,lw=.6,alpha=.8)
    bb=n(core[k]['total_b']);pp=np.linspace(max(.8,prec.min()*.6),prec.max()*1.1,80)
    ax.plot(bb-1.96/pp,pp,color=LIGHT,lw=1);ax.plot(bb+1.96/pp,pp,color=LIGHT,lw=1);ax.axvline(bb,color=BLUE,lw=1.3)
    ax.set_title(SH[k],loc='left',weight='bold',fontsize=10);ax.tick_params(labelsize=8)
fig.suptitle('Funnel plots of SNP-specific Wald ratios',x=.05,ha='left',weight='bold',fontsize=12)
fig.supxlabel('SNP-specific log OR per SD smoking',y=.04,fontsize=10);fig.supylabel('Inverse SE of Wald ratio',x=.01,fontsize=10)
finalize_figure(fig,'FigS8')

expected=[f'Fig{i}' for i in range(1,5)]+[f'FigS{i}' for i in range(1,10)]
for stem in expected:
    for ext in ['png','pdf','svg']:assert (FIG/f'{stem}.{ext}').is_file()
(FIG/'figure_numeric_payload.json').write_text(json.dumps(payload,indent=2),encoding='utf-8')
(OUT/'绘图Skill与风格说明.txt').write_text('Skill: '+str(ROOT.parent/'tools/figures4papers/scientific-figure-making/SKILL.md')+'\nSource: https://github.com/ChenLiu-1996/figures4papers\n字体Arial；吸烟深蓝#0F4D92；BMI青绿#42949E；Egger暗红#B64342；总效应及参考线灰色。\n主图重组为设计图、总体森林图、双向结果图和相互调整森林图；散点独立移至S7，增加S8漏斗和S9逐染色体剔除。所有图PNG 600 dpi、PDF嵌入字体、SVG保留可编辑文字。\n所有数值来自冻结分析及新增确定性核查，没有改变主分析工具选择。\n',encoding='utf-8-sig')
print(f'Created {len(expected)} figures in PNG / PDF / SVG')
