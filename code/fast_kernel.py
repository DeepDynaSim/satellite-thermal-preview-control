"""Numba acceleration of the explicitly specified sampled simulation; no fast-math."""
import numpy as np
from numba import njit

@njit(cache=True)
def thermal(T,tau,direction,eta,flux,p):
    ca,G,alpha,S,cosbeta,es,w,IR=p
    z=np.tanh(tau/w);chi=(1+direction*z)/2;dc=direction*(1-z*z)/(2*w)
    solar=(alpha-eta)*S*flux*cosbeta*chi
    td=np.empty(2)
    td[0]=(solar+IR-G*(T[0]-T[1])-es*(T[0]**4-81.))/ca
    td[1]=(IR+G*(T[0]-T[1])-es*(T[1]**4-81.))/ca
    tdd=np.empty(2)
    tdd[0]=(-G-4*es*T[0]**3)/ca*td[0]+G/ca*td[1]+(alpha-eta)*S*flux*cosbeta*dc/ca
    tdd[1]=G/ca*td[0]+(-G-4*es*T[1]**3)/ca*td[1]
    return td,tdd,dc

@njit(cache=True)
def heat_next(T,tau,direction,eta,flux,dt,p):
    k1=thermal(T,tau,direction,eta,flux,p)[0]
    k2=thermal(T+dt*k1/2,tau+dt/2,direction,eta,flux,p)[0]
    k3=thermal(T+dt*k2/2,tau+dt/2,direction,eta,flux,p)[0]
    k4=thermal(T+dt*k3,tau+dt,direction,eta,flux,p)[0]
    return T+dt*(k1+2*k2+2*k3+k4)/6

@njit(cache=True)
def projection(candidate,a,b,baseline,limits,rdiag):
    u=np.minimum(np.maximum(candidate,-limits),limits)
    if a@u<=b+1e-10:return u,0.,0.
    base=np.minimum(np.maximum(baseline,-limits),limits)
    if a@base<=b+1e-10:
        alpha=min(1.,max(0.,(b-a@base)/max(a@(u-base),1e-30)))
        return base+alpha*(u-base),1.,0.
    minv=-np.sum(np.abs(a)*limits)
    if minv>b+1e-10:return -limits*np.sign(a),1.,minv-b
    lo=0.;hi=1.
    for _ in range(60):
        test=np.minimum(np.maximum(candidate-hi*a/rdiag,-limits),limits)
        if a@test<=b:break
        hi*=2
    for _ in range(32):
        mid=(lo+hi)/2;test=np.minimum(np.maximum(candidate-mid*a/rdiag,-limits),limits)
        if a@test>b:lo=mid
        else:hi=mid
    return np.minimum(np.maximum(candidate-hi*a/rdiag,-limits),limits),1.,0.

@njit(cache=True)
def run_core(kind,steps,dt,te,direction,eta,flux,ptrue,pnom,T,y,
    Ad,Bd,Ed,rqplant,tip,selector,rq,F1,F2,A,B,K,P,Q,W,R,sc,limits,portgain,
    w0,b0,w1,b1,w2,b2,yscale,oracle,learn_bias):
    rec=np.empty((steps+1,18));cost=0.;active=0.;infeasible=0
    nd=len(rq);d=len(rqplant);rd=np.diag(R).copy()
    for j in range(steps+1):
        t=j*dt;tau=t-te
        td,tdd,dc=thermal(T,tau,direction,eta,flux,pnom)
        delta=T[0]-T[1];dd=td[0]-td[1];dda=tdd[0]-tdd[1]
        rr=np.concatenate((rq*delta,rq*dd));yc=selector@y;x=yc-rr
        f=F1*dd+F2*dda;ulqr=-K@x
        if kind==0:
            u=np.zeros(4);u[0]=-150*y[0]-550*y[d]
        elif kind==1:u=-K@yc
        elif kind==2:u=ulqr.copy()
        elif kind==3:
            index=min(int(t/.1),len(oracle)-2);frac=(t-index*.1)/.1
            u=ulqr+(1-frac)*oracle[index]+frac*oracle[index+1]
        else:
            features=np.array([(T[0]+T[1]-560)/110,delta/40,np.tanh(tau/18),direction,(eta-.30)/.08,(flux-1.)/.2])
            pred=(np.tanh(np.tanh(features@w0+b0)@w1+b1)@w2+b2)*yscale
            activity=np.sum(td*td)+(.3*dc)**2;gate=activity/(activity+1e-8)
            u=ulqr+gate*pred
            if learn_bias!=0:
                for k in range(1,4):u[k]+=learn_bias*np.sin(.7*t)
        auxiliary=np.zeros(4) if kind==0 else -portgain@y
        u=u+auxiliary
        fe=f+B@auxiliary
        a=2*B.T@P@x
        b=-.25*(x@W@x)+2*abs(x@P@fe)-2*x@P@(A@x+f)
        act=0.;slack=0.
        if kind==5:u,act,slack=projection(u,a,b,ulqr+auxiliary,limits,rd)
        else:u=np.minimum(np.maximum(u,-limits),limits)
        error=y[:d]-rqplant*delta
        vibration=tip@error[1:]
        static_tip=tip@rqplant[1:]*delta;actual_tip=tip@y[1:d]
        vdot=2*x@P@(A@x+B@u+f)
        margin=vdot+.25*(x@W@x)-2*abs(x@P@fe)
        running=x@Q@x+u@R@u
        rec[j,0]=t;rec[j,1]=T[0];rec[j,2]=T[1];rec[j,3]=delta
        rec[j,4]=static_tip;rec[j,5]=actual_tip;rec[j,6]=vibration;rec[j,7]=y[0];rec[j,8]=y[d]
        rec[j,9:13]=u;rec[j,13]=running;rec[j,14]=act;rec[j,15]=slack;rec[j,16]=margin;rec[j,17]=np.sqrt(np.sum((x/sc)**2))
        active+=act
        if slack>1e-8:infeasible+=1
        if j<steps:
            cost+=dt*running
            Tn=heat_next(T,tau,direction,eta,flux,dt,ptrue)
            y=Ad@y+Bd@u+Ed*((delta+Tn[0]-Tn[1])/2)
            T=Tn
    return rec,cost,active,infeasible
