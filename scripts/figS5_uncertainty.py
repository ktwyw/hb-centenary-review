"""Fig. S5 of the supplementary material (figS5_uncertainty.pdf): Constrained yield-stress profiles, known additive noise, fixed synthetic data."""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from hbstyle import *
rng=np.random.default_rng(20260919);g=np.logspace(-1.5,3,36);sigma=1.
y=10+2*np.sqrt(g)+rng.normal(0,sigma,len(g));tys=np.linspace(0,20,241)
gp=np.logspace(-3,3,250);fig,ax=plt.subplots(1,2,figsize=(6.6,2.75))
for mask,col,lab,ls in [(np.ones(len(g),dtype=bool),C['blue'],'broad window','-'),(g>=10,C['red'],'high-rate subset','--')]:
    vals=[];pred=[]
    for ty in tys:
        fit=least_squares(lambda z:(ty+z[0]*g[mask]**z[1]-y[mask])/sigma,[2,.5],bounds=([1e-6,.1],[100,1.5]),ftol=1e-10,xtol=1e-10,gtol=1e-10)
        vals.append(np.sum(fit.fun**2));pred.append(ty+fit.x[0]*gp**fit.x[1])
    vals=np.array(vals);pred=np.array(pred);delta=vals-vals.min();valid=delta<=3.84
    ax[0].plot(tys,delta,color=col,ls=ls,label=lab)
    ax[1].fill_between(gp,pred[valid].min(0),pred[valid].max(0),color=col,alpha=.15,gid=f'{lab}: range of curves inside the 95% profile interval')
    ax[1].semilogx(gp,pred[np.argmin(vals)],color=col,ls=ls,label=lab)
ax[0].axhline(3.84,color='grey',lw=.7,ls=':',gid='95% threshold, delta chi2 = 3.84');ax[0].axvline(10,color='k',lw=.7,ls=':',gid='true sigma_y')
ax[0].text(.6,4.3,r'$\Delta\chi^2=3.84$',fontsize=9.5)
ax[0].set(xlabel=r'fixed $\sigma_y$ (Pa)',ylabel=r'profile $\Delta\chi^2$',ylim=(0,18));
ax[1].semilogx(g,y,'o',ms=2.4,color='k',label='synthetic observations');ax[1].semilogx(gp,10+2*np.sqrt(gp),':',color='k',lw=.8,label='generating law')
ax[1].set(xlabel=r'$\dot\gamma$ (s$^{-1}$)',ylabel=r'$\sigma$ (Pa)',ylim=(0,80));ax[1].legend(loc='upper left',fontsize=9.5)
label_panels(ax,x=-.20);fig.tight_layout(w_pad=2.6);save(fig, "figS5_uncertainty")
