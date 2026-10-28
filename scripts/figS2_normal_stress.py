"""Fig. S2 of the supplementary material (figS2_normal_stress.pdf): Supplementary normal-stress calculation from the fits of de Cagny et al.

N_i = psi_i sigma^2 with sigma the shear stress; sigma_d = sqrt(tau_d:tau_d/2).

Author manuscript arXiv:1805.12461v2, Table 1 and SI Table A; journal identity
J. Rheol. 63, 285-290 (2019), doi:10.1122/1.5063796. No measured points plotted.
"""
from pathlib import Path
import json
import numpy as np
from hbstyle import plt,save,C,label_panels
ROOT=Path(__file__).resolve().parents[1]
values=[('Emulsion',53.56,27.,.42,3.19e-3,-2.32e-3),('Hair gel',56.5,18.,.44,3.22e-3,-1.57e-3),('Carbopol',65.51,28.13,.40,2.99e-3,-2.58e-3)]
def invariant(s,n1,n2):return np.sqrt(s*s+(n1*n1+n1*n2+n2*n2)/3)
def run():
    s=np.geomspace(20,300,220);fig,ax=plt.subplots(1,2,figsize=(6.6,3.15));meta=[];output=[s]
    for (name,ty,k,n,p1,p2),c,ls in zip(values,[C['blue'],C['red'],C['green']],['-','--','-.']):
        n1=p1*s*s;n2=p2*s*s
        ax[0].loglog(s,n1,c=c,ls=ls,lw=1.6,label=name,gid=f'{name}: N1 (thick)')
        ax[0].loglog(s,-n2,c=c,ls=ls,lw=0.7,gid=f'{name}: -N2 (thin)')
        effect=100*(invariant(s,n1,n2)/s-1)
        ax[1].semilogx(s,effect,c=c,ls=ls,label=name)
        e_y=100*(invariant(ty,p1*ty**2,p2*ty**2)/ty-1)
        ax[1].plot(ty,e_y,'o',c=c,mfc='white',ms=4,gid=f'{name}: value at sigma = sigma_y')
        meta.append({'psi2_over_psi1':p2/p1,'N1_over_sigma_at_intercept':p1*ty,'material':name,'sigma_y_Pa':ty,'K_Pa_s_n':k,'n':n,'psi1_Pa_inv':p1,'psi2_Pa_inv':p2,'normal_N1_at_HB_intercept_Pa':p1*ty**2,'normal_N2_at_HB_intercept_Pa':p2*ty**2,'equivalent_stress_increment_at_intercept_percent':e_y})
        output.extend([n1,n2,effect])
    ax[0].set(xlabel=r'Shear stress $\sigma$ (Pa)',ylabel=r'$N_1$ (thick), $-N_2$ (thin) (Pa)');ax[0].legend(fontsize=9,loc='lower right')
    ax[0].text(.03,.94,'Isotropic HB: $N_1=N_2=0$\nUC Saramito: $N_2=0$',transform=ax[0].transAxes,va='top',fontsize=9)
    ax[1].set(xlabel=r'Shear stress $\sigma$ (Pa)',ylabel=r'$100(\sigma_d/\sigma-1)$ (%)')
    ax[1].text(.04,.94,r'Circles: $\sigma=\sigma_y$',transform=ax[1].transAxes,va='top',fontsize=9)
    label_panels(ax,x=-.17)
    from matplotlib.ticker import NullFormatter
    for a in ax:
        a.xaxis.set_minor_formatter(NullFormatter()); a.set_xticks([20,50,100,200]); a.set_xticklabels(['20','50','100','200'])
    fig.tight_layout();save(fig, "figS2_normal_stress")
    (ROOT/'data/normal_stress_reconstruction.json').write_text(json.dumps({'source':'de Cagny et al., J. Rheol. 63, 285 (2019), Table 1 and SI Table A','doi':'10.1122/1.5063796','meaning':'Consequences of the published fits; stress interval 20-300 Pa.','materials':meta},indent=2))
    np.savetxt(ROOT/'data/normal_stress_curves.csv',np.column_stack(output),delimiter=',',header='shear_Pa,'+','.join(f'{x}_{v}' for x in ['emulsion','hairgel','carbopol'] for v in ['N1_Pa','N2_Pa','eq_increment_percent']),comments='')
    print(json.dumps(meta,indent=2))
if __name__=='__main__':run()
