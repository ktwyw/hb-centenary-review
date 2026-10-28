"""Fig. 9 of the article (fig09_transport_exponents.pdf): Analytical exponent crossovers; no experimental data."""
import numpy as np
import matplotlib.pyplot as plt
from hbstyle import *
from verify_equations import pipe_q

fig,ax=plt.subplots(1,2,figsize=(6.6,2.65))
e=np.logspace(-5,4,700);n=.5
q=np.array([pipe_q(1+x,1,1,n) for x in e]);s=np.gradient(np.log(q),np.log(e))
ax[0].semilogx(e,s,color=C['blue'],label='exact pipe solution')
ax[0].axhline(3,color='k',ls=':',lw=.7,gid='near-yield limit 1+1/n');ax[0].axhline(2,color='k',ls='--',lw=.7,gid='far-from-yield limit 1/n')
ax[0].text(2e-5,3.02,r'near yield: $1+1/n=3$',fontsize=9.5,va='bottom')
ax[0].text(.003,1.98,r'far above yield: $1/n=2$',fontsize=9.5,va='top')
ax[0].set(xlabel=r'$(G-G_c)/G_c$',ylabel=r'$\mathrm{d}\log Q/\mathrm{d}\log(G-G_c)$',ylim=(1.8,3.2))
x=np.logspace(-4,5,600)
for v,c,ls in zip([.01,.1,1],[C['purple'],C['blue'],C['red']],[':','--','-']):
    ne=(.5*x**.5+v*x)/(x**.5+v*x)
    ax[1].semilogx(x,ne,color=c,ls=ls,label=rf'$\eta_s\dot\gamma_*/A={v:g}$')
ax[1].set(xlabel=r'$\dot\gamma/\dot\gamma_*$',ylabel=r'$n_{\mathrm{eff}}$',ylim=(.45,1.04));ax[1].legend(loc='upper left',fontsize=9.5)
label_panels(ax,x=-.20);fig.tight_layout(w_pad=2.6);save(fig, "fig09_transport_exponents")
