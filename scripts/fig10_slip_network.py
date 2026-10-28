"""Fig. 10 of the article (fig10_slip_network.pdf): Exact capillary illustrations: wall slip and independent parallel channels."""
import numpy as np
import matplotlib.pyplot as plt
from hbstyle import *
from verify_equations import pipe_q

fig,ax=plt.subplots(1,2,figsize=(6.6,2.65));tw=np.linspace(0,2.8,500)
bulk=np.array([pipe_q(t,1,1,.5)/np.pi for t in tw]);slip=.3*np.maximum(tw-.4,0)
ax[0].plot(tw,bulk,color='k',ls='--',label='bulk shear only')
ax[0].plot(tw,slip,color=C['red'],ls=':',label='wall sliding')
ax[0].plot(tw,bulk+slip,color=C['blue'],label='total transport')
ax[0].axvspan(.4,1,color=C['blue'],alpha=.09,lw=0,gid='rigid sliding');ax[0].axvline(1,color='grey',lw=.5,gid='bulk yield, sigma_w = sigma_y')
ax[0].text(.66,.55,'rigid\nsliding',ha='center',fontsize=8.25)
ax[0].set(xlabel=r'$\sigma_w/\sigma_y$',ylabel=r'$U/[R(\sigma_y/K)^{1/n}]$',ylim=(0,1.5),xlim=(0,4.75),xticks=[0,1,2,3,4])
for yv,lab,c in [((bulk+slip)[-1],'total transport',C['blue']),(slip[-1],'wall sliding',C['red']),(bulk[-1],'bulk shear only','k')]:
    ax[0].text(2.88,yv,lab,color=c,fontsize=9.5,va='center',ha='left')
g=np.linspace(0,5,500);radii=[1,.8,.65,.5,.35];total=np.zeros_like(g)
for j,R in enumerate(radii):
    qi=np.array([pipe_q(v*R,1,1,.5,R) for v in g])/np.pi
    total+=qi
    ax[1].plot(g,qi,color=C['grey'],lw=.6,alpha=.8,gid=f'channel R/Rmax={R:g}')
ax[1].plot(g,total,color=C['blue'],lw=1.4,label='sum of five channels')
ax[1].plot([],[],color=C['grey'],lw=.7,label='individual channels')
for R in radii: ax[1].plot([1/R],[0],ls='none',marker='|',ms=9,mew=1.4,color=C['red'],clip_on=False,gid=f'onset of channel R/Rmax={R:g}')
ax[1].set(xlabel=r'$G/G_{c,\mathrm{widest}}$',ylabel=r'$Q_{\mathrm{tot}}/[\pi R_{\max}^3(\sigma_y/K)^{1/n}]$',ylim=(0,6))
ax[1].legend(loc='upper left',fontsize=9.5);label_panels(ax,x=-.20);fig.tight_layout(w_pad=2.6);save(fig, "fig10_slip_network")
