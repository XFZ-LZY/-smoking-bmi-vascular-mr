from pathlib import Path
import sys,csv
ROOT=Path(__file__).resolve().parent.parent

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import font_manager,ticker
OUT=ROOT/'audit';FIG=OUT/'figures';FIG.mkdir(exist_ok=True)
for p in Path('C:/Windows/Fonts').glob('arial*.ttf'):font_manager.fontManager.addfont(str(p))
plt.rcParams.update({'font.family':['Arial','DejaVu Sans'],'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
rows=list(csv.DictReader(open(OUT/'independence_sensitivity_results.csv',encoding='utf-8-sig')))
keys=['cad','pad_finngen_r5','ais','laa','svs','ces'];short=['CAD','PAD','AIS','LAA','SVS','CES']
schemes=['all','primary','expanded'];marks=['s','o','^'];labels=['All instruments','Shared-pathway screen','Expanded cardiometabolic screen']
fig=plt.figure(figsize=(10.6,9.7))
gs=fig.add_gridspec(2,3,width_ratios=[.65,2,2.65],hspace=.42,wspace=.07,left=.05,right=.99,bottom=.1,top=.91)
for rowno,exposure in enumerate(['smoking','bmi']):
 color=['#0F4D92','#42949E'][rowno];la=fig.add_subplot(gs[rowno,0]);ax=fig.add_subplot(gs[rowno,1]);ta=fig.add_subplot(gs[rowno,2])
 for a in [la,ax,ta]:a.set_ylim(-.8,6.4);a.set_yticks([])
 la.axis('off');ta.axis('off');ax.spines['left'].set_visible(False);ax.set_xscale('log');ax.axvline(1,color='#767676',ls='--',lw=.9)
 ax.set_xlim(.68,4.2);ax.set_xticks([.75,1,1.5,2,3,4]);ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v,p:f'{v:g}'));ax.xaxis.set_minor_locator(ticker.NullLocator());ax.set_xlabel('Odds ratio (95% CI)')
 la.text(0,6.0,'Outcome',weight='bold');ta.text(0,6.0,'OR (95% CI)',weight='bold');ta.text(.63,6.0,'FDR P',weight='bold');ta.text(.87,6.0,'SNPs',weight='bold')
 for i,(key,label) in enumerate(zip(keys,short)):
  y=5-i
  if i%2==0:
   for a in [la,ax,ta]:a.axhspan(y-.47,y+.47,color='#F3F5F7',zorder=0)
  la.text(0,y,label,va='center',weight='bold')
  for j,s in enumerate(schemes):
   r=next(r for r in rows if r['scheme']==s and r['key']==key and r['model']=='multivariable' and r['exposure']==exposure)
   b,lo,hi,p=map(float,[r['OR'],r['OR_lo'],r['OR_hi'],r['fdr']]);yy=y+.25-j*.25
   ax.errorbar(b,yy,xerr=[[b-lo],[hi-b]],fmt=marks[j],ms=5,ecolor=color,mec=color,mfc=color if j==0 else 'white',lw=1.2,capsize=2)
   ta.text(0,yy,f'{b:.2f} ({lo:.2f}–{hi:.2f})',va='center',fontsize=8.6)
   ta.text(.63,yy,f'{p:.2e}' if p<.001 else f'{p:.3f}',va='center',fontsize=8.6)
   ta.text(.90,yy,r['n'],va='center',fontsize=8.6)
 fig.text(.05,.955-rowno*.467,('a  Smoking adjusted for BMI' if rowno==0 else 'b  BMI adjusted for smoking'),weight='bold',fontsize=12)
fig.legend([Line2D([],[],marker=m,color='#767676',mfc='#767676' if i==0 else 'white',ls='none') for i,m in enumerate(marks)],labels,loc='lower center',bbox_to_anchor=(.5,.023),ncol=3,frameon=False,fontsize=9)
for ext in ['png','pdf','svg']:fig.savefig(FIG/('FigS10_independence.'+ext),dpi=600,bbox_inches='tight',pad_inches=.1)
plt.close(fig)
