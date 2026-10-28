"""Deterministic equation checks for the review; not a CFD validation suite."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

def pipe_q(tw, ty, K, n, R=1.):
    if tw <= ty: return 0.
    d=tw-ty
    return np.pi*R**3/(K**(1/n)*tw**3)*(d**(3+1/n)/(3+1/n)+2*ty*d**(2+1/n)/(2+1/n)+ty**2*d**(1+1/n)/(1+1/n))

def pipe_phi(tw,ty,K,n,R=1.):
    if tw<=ty: return 0.
    p=ty/tw
    return np.pi*R**3*n*(tw/K)**(1/n)*(1-p)**(1+1/n)*((1-p)**2/(3*n+1)+2*p*(1-p)/(2*n+1)+p*p/(n+1))

def thin_q(G,h,ty,K,n):
    """Depth-integrated thin-layer flux as printed in Eq. (thin_flux)."""
    a=max(h-ty/G,0.)
    return (G/K)**(1/n)*(a**(2+1/n)/(2+1/n)+(h-a)*a**(1+1/n)/(1+1/n))

def couette_omega(T,ty,K,n,Ri,Ro):
    """Inner angular velocity from Eq. (couette), stationary outer cylinder, torque per unit length T."""
    if T<=2*np.pi*Ri**2*ty: return 0.
    ry=np.sqrt(T/(2*np.pi*ty)) if ty>0 else np.inf
    return quad(lambda r:((T/(2*np.pi*r*r)-ty)/K)**(1/n)/r,Ri,min(ry,Ro),epsabs=1e-14)[0]

def check_close(a,b,tol=1e-9):
    assert np.isclose(a,b,rtol=tol,atol=1e-12), (a,b)

def run():
    checks=[]
    for n in [.35,.5,.7,1.,1.4]:
        for phi in [0.,.1,.5,.9,.99]:
            tw=2.;ty=phi*tw;K=.8;R=.7
            q=pipe_q(tw,ty,K,n,R)
            check_close(q,pipe_phi(tw,ty,K,n,R))
            rp=phi*R
            qint=np.pi*quad(lambda r:((tw*r/R-ty)/K)**(1/n)*r*r,rp,R,epsabs=1e-12)[0]
            check_close(q,qint,2e-8)
            G=2*tw/R
            diss=2*np.pi*quad(lambda r:(ty*((tw*r/R-ty)/K)**(1/n)+K*((tw*r/R-ty)/K)**((n+1)/n))*r,rp,R,epsabs=1e-12)[0]
            check_close(G*q,diss,2e-8)
    checks+=['25 pipe-formula equivalence cases','25 direct pipe quadratures','25 pressure-work/dissipation balances']
    G,R,mu=3.,.4,2.
    check_close(pipe_q(G*R/2,0,mu,1,R),np.pi*R**4*G/(8*mu))
    checks.append('Newtonian Poiseuille limit')
    for phi in [.1,.4,.8]:
        tw=2.; check_close(pipe_q(tw,phi*tw,mu,1,R),np.pi*R**3*tw/(4*mu)*(1-4*phi/3+phi**4/3))
    checks.append('Bingham Buckingham-Reiner limit')
    for n in [.4,.7,1.]:
        e1,e2=1e-5,2e-5
        slope=np.log(pipe_q(1+e2,1,1,n)/pipe_q(1+e1,1,1,n))/np.log(e2/e1)
        check_close(slope,1+1/n,2e-5)
    checks.append('Near-onset exponents')
    for n in [.4,.7,1.]:
        ty,K,Ri,Ro,delta=2.,.7,1.,2.,1e-5
        ti=ty+delta;to=ti*(Ri/Ro)**2
        omega=quad(lambda t:(t-ty)**(1/n)/t,max(ty,to),ti,epsabs=1e-30)[0]/(2*K**(1/n))
        asym=n*delta**(1+1/n)/(2*(n+1)*K**(1/n)*ty)
        check_close(omega/asym,1.,1e-5)
    checks.append('Partly yielded Couette near-onset asymptotics')
    radius=.8;yield_stress=2.
    torque=2*np.pi*yield_stress*radius**3*quad(lambda t:np.sin(t)**2,0,np.pi)[0]
    check_close(torque,np.pi**2*yield_stress*radius**3)
    checks.append('Ideal spherical cavern torque prefactor')
    rho,U,D=1200.,1.3,.1;n=1.;K=.5
    re=rho*U**(2-n)*D**n/(K*8**(n-1)*((3*n+1)/(4*n))**n)
    check_close(re,rho*U*D/K); checks.append('Newtonian Metzner-Reed limit')
    # Thin-layer flux, Eq. (thin_flux): Newtonian limit of the printed formula, then the yield
    # case against direct double quadrature of the basal-layer velocity profile.
    h=.3;check_close(thin_q(G,h,0.,mu,1.),G*h**3/(3*mu)); checks.append('Newtonian film limit of Eq. thin_flux')
    nfilm=0
    for n in [.4,.7,1.,1.3]:
        for ty in [0.,.2,.6]:
            G,h,K=2.,.5,.9
            a=max(h-ty/G,0.)
            rate=lambda z:((G*(h-z)-ty)/K)**(1/n) if z<a else 0.
            u=lambda z:quad(rate,0,min(z,a),epsabs=1e-13)[0]
            qnum=quad(u,0,a,epsabs=1e-13)[0]+(h-a)*u(a)
            check_close(thin_q(G,h,ty,K,n),qnum,1e-7); nfilm+=1
    checks.append(f'{nfilm} thin-layer flux quadratures with yield (Eq. thin_flux)')
    # Plane channel flux (Sec. VII A) against quadrature of the printed velocity profile.
    nch=0
    for n in [.4,.7,1.,1.4]:
        for phi in [0.,.2,.6,.9]:
            tw,K,H=2.,.8,.6;ty=phi*tw
            uy=lambda y:n/(n+1)*H/tw*(1/K)**(1/n)*((tw-ty)**((n+1)/n)-(max(tw*abs(y)/H-ty,0.))**((n+1)/n))
            qnum=2*quad(uy,0,H,epsabs=1e-13)[0]
            qform=2*H**2*n*(tw/K)**(1/n)*(1-phi)**((n+1)/n)*((1-phi)/(2*n+1)+phi/(n+1))
            check_close(qform,qnum,1e-8); nch+=1
    checks.append(f'{nch} plane-channel flux quadratures (Sec. VII A)')
    # Couette, Eq. (couette): full relation against closed forms (Newtonian, power law, Bingham),
    # for fully and partly yielded gaps.
    Ri,Ro=1.,1.4
    for T in [10.,30.]:
        check_close(couette_omega(T,0.,mu,1.,Ri,Ro),T/(4*np.pi*mu)*(1/Ri**2-1/Ro**2))
        for n in [.5,.8]:
            K=.7
            check_close(couette_omega(T,0.,K,n,Ri,Ro),n/2*(T/(2*np.pi*K))**(1/n)*(Ri**(-2/n)-Ro**(-2/n)),1e-8)
        for ty in [.5,2.,4.]:
            if T<=2*np.pi*Ri**2*ty: continue
            rm=min(np.sqrt(T/(2*np.pi*ty)),Ro)
            bing=T/(4*np.pi*mu)*(1/Ri**2-1/rm**2)-ty/mu*np.log(rm/Ri)
            check_close(couette_omega(T,ty,mu,1.,Ri,Ro),bing,1e-8)
    checks.append('Full Couette relation: Newtonian, power-law and Bingham (full and partial yielding)')
    # Tensor algebra of the 3D HB law for uniaxial and planar extension, Eqs. (uniaxial_hb), (planar_hb).
    for n in [.4,1.,1.5]:
        for ty in [0.,3.]:
            K,ed=1.7,.8
            for D,sigma in [(np.diag([-ed/2,-ed/2,ed]),lambda t:t[2,2]-t[0,0]),(np.diag([ed,-ed,0.]),lambda t:t[0,0]-t[1,1])]:
                gd=np.sqrt(2*np.sum(D*D)); tau=(ty/gd+K*gd**(n-1))*2*D
                expected=(np.sqrt(3)*ty+3**((n+1)/2)*K*ed**n) if D[2,2]!=0 else (2*ty+2**(n+1)*K*ed**n)
                check_close(sigma(tau),expected)
    ed=.8;De=np.diag([-ed/2,-ed/2,ed]);gd=np.sqrt(2*np.sum(De*De));tau=2*mu*De
    check_close(gd,np.sqrt(3)*ed);check_close(tau[2,2]-tau[0,0],3*mu*ed)
    checks.append('Uniaxial and planar extension stresses with yield (Eqs. uniaxial_hb, planar_hb) and Newtonian limit')
    # Normal-stress invariant, Eq. (normal_invariant), from a traceless deviator (Supplementary Note S1).
    for s,N1,N2 in [(50.,3.,-1.5),(60.,10.,-7.),(1.,0.,0.)]:
        S=np.array([[(2*N1+N2)/3,s,0.],[s,(-N1+N2)/3,0.],[0.,0.,(-N1-2*N2)/3]])
        check_close(np.trace(S),0.); check_close(.5*np.sum(S*S),s*s+(N1*N1+N1*N2+N2*N2)/3)
    checks.append('Normal-stress invariant identity (Eq. normal_invariant)')
    # Spreading (slump) coefficient C=225/(128 pi^2): volume of h^2=2 ty (R-r)/(rho g) by quadrature.
    rho_c,g_c,R_c,ty_c=2300.,9.81,.35,20.
    V=quad(lambda r:2*np.pi*r*np.sqrt(2*ty_c*(R_c-r)/(rho_c*g_c)),0,R_c,epsabs=1e-14)[0]
    check_close(ty_c,225/(128*np.pi**2)*rho_c*g_c*V**2/R_c**5,1e-8)
    checks.append('Spreading coefficient C=225/(128 pi^2) by volume quadrature')
    # Sphere yield-number conventions: Y_g=3 ty/(2 drho g R) and Y_R=ty/(drho g R), so Y_g=3 Y_R/2.
    ty_s,drho,R_s=3.5,1500.,.0025
    check_close(3*ty_s/(2*drho*9.81*R_s),1.5*ty_s/(drho*9.81*R_s))
    checks.append('Sphere yield-number conventions Y_g=3Y_R/2')
    # Plug-radius diagnostic, Eq. (rate_plug): radius at which the HB-Papanastasiou rate equals eps.
    for m in [10,1000]:
        eps=.003; ty,K,n,tw=1.,1.,.5,2.
        rad=(ty*(-np.expm1(-m*eps))+K*eps**n)/tw
        recovered=brentq(lambda g:ty*(-np.expm1(-m*g))+K*g**n-tw*rad,0,1,xtol=1e-14)
        check_close(recovered,eps);checks.append(f'Rate-threshold inversion m={m}')
    Q,L,D,ty,K,n,rho=.01,500,.1,5,.5,.7,1200
    tw=brentq(lambda s:pipe_q(s,ty,K,n,D/2)-Q,ty+1e-10,1000,xtol=1e-13)
    U=Q/(np.pi*(D/2)**2)
    eps=1e-5
    slope_Q=(np.log(pipe_q(tw*(1+eps),ty,K,n,D/2))-np.log(pipe_q(tw*(1-eps),ty,K,n,D/2)))/(np.log(1+eps)-np.log(1-eps))
    np_local=1/slope_Q; Kp=tw/(8*U/D)**np_local; Re=8*rho*U**2/tw
    examples={'pipe':{'wall_stress_Pa':tw,'pressure_drop_kPa':4*tw*L/D/1000,'mean_speed_m_s':U,'plug_radius_mm':D/2*ty/tw*1000,'local_n':np_local,'local_K':Kp,'Re_MR':Re},
       'spread':{'yield_stress_R035_Pa':225*2300*9.81*.0055**2/(128*np.pi**2*.35**5),'yield_stress_R030_Pa':225*2300*9.81*.0055**2/(128*np.pi**2*.30**5)},
       'sphere':{'Yg_critical':.143,'YR_critical':.143/1.5,'yield_stress_Pa':(.143/1.5)*1500*9.81*.0025}}
    check_close(4*tw*L/D, (16/Re)*2*rho*U*U*L/D)
    checks.append('Fanning pressure-loss identity and worked examples')
    return {'status':'PASS','checks':checks,'worked_examples':examples,'note':'Analytical and synthetic regression tests only; no experimental or multidimensional CFD validation.'}

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',default=str(Path(__file__).resolve().parents[1]/'results'/'equation_verification.json'),
                    help='where to write the JSON record (default: results/equation_verification.json)')
    ap.add_argument('--no-write',action='store_true',help='run the checks and print the result without writing any file')
    args=ap.parse_args()
    result=run()
    if not args.no_write: Path(args.out).parent.mkdir(exist_ok=True); Path(args.out).write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
