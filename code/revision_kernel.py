"""Streaming, deterministic benchmark kernel. No fast-math; native-step metrics."""
import numpy as np
from numba import njit
from fast_kernel import thermal, heat_next, projection

@njit(cache=True)
def features(T,tau,direction,eta,flux):
    return np.array([(T[0]+T[1]-560)/110,(T[0]-T[1])/40,np.tanh(tau/18),direction,(eta-.3)/.08,(flux-1)/.2])

@njit(cache=True)
def forecast_one(T,tau,direction,eta,flux,p,kernels,h):
    state=T.copy();out=np.zeros(4)
    for j in range(len(kernels)):
        mid=heat_next(state,tau+j*h,direction,eta,flux,h/2,p)
        td,tdd,_=thermal(mid,tau+(j+.5)*h,direction,eta,flux,p)
        out+=kernels[j,:,0]*(td[0]-td[1])+kernels[j,:,1]*(tdd[0]-tdd[1])
        state=heat_next(state,tau+j*h,direction,eta,flux,h,p)
    return out

@njit(cache=True)
def forecast_many(T,tau,direction,eta,flux,p,kernels,h):
    out=np.zeros((len(T),4))
    for i in range(len(T)):
        out[i]=forecast_one(T[i],tau[i],direction[i],eta[i],flux[i],p,kernels,h)
    return out

@njit(cache=True)
def polynomial(xx,coef,powers):
    v=np.ones(len(powers))
    for i in range(len(powers)):
        for j in range(6):v[i]*=xx[j]**powers[i,j]
    return v@coef

@njit(cache=True,nogil=True)
def run_stream(kind,steps,dt,te,direction,eta,flux,ptrue,pnom,T,y,
    Ad,Bd,Ed,rqplant,tip,selector,rq,F1,F2,A,B,K,P,Q,W,R,sc,limits,portgain,
    w0,b0,w1,b1,w2,b2,yscale,activation,oracle,learn_bias,gate_eps,gate_on,filter_on,
    coef,powers,stride,core_delay,residual_delay,driver_alpha,noise,observer,
    obsAd,obsBd,obsEd,obsL,obsC,sensor,normals):
    # kind 0 hub PD, 1 fixed target LQR, 2 thermal LQR, 3 precomputed preview,
    # 4 neural, 6 cubic polynomial, 7 collocated passive damping with hub PD.
    nd=len(rq);d=len(rqplant);rd=np.diag(R).copy()
    rec=np.empty((steps//stride+1,18));nr=0
    cost=0.;sumtip=0.;sumangle=0.;sumu=np.zeros(4);peaku=np.zeros(4)
    sat=np.zeros(4);clip=np.zeros(4);runs=np.zeros(4);maxrun=np.zeros(4)
    active=0.;infeasible=0.;sumcor=0.;maxcor=0.;maxmargin=-1e300;trueviol=0.
    zhat=selector@y;uactual=np.zeros(4)
    corebuf=np.zeros((max(1,core_delay+1),4));resbuf=np.zeros((max(1,residual_delay+1),4))
    trueerrmax=0.;failure=0
    for j in range(steps+1):
        t=j*dt;tau=t-te
        # Fixed deterministic draws indexed at 1 ms. Biases persist within a run.
        ix=min(int(t/.001),len(normals)-1)
        Tm=T+noise[:2]+noise[2:4]*normals[ix,:2]
        td,tdd,dc=thermal(Tm,tau,direction,eta,flux,pnom)
        delta=T[0]-T[1];dd=td[0]-td[1];dda=tdd[0]-tdd[1]
        rr=np.concatenate((rq*(Tm[0]-Tm[1]),rq*dd));yc=selector@y
        if observer:
            measurement=sensor@y+noise[4:9]+noise[9:14]*normals[ix,2:7]
            # Innovation correction only at independent sensor update boundaries.
            if j==0 or int((t-dt)/.001)!=int(t/.001):zhat+=obsL@(measurement-obsC@zhat)
            state=zhat
        else:state=yc
        x=state-rr;xt=yc-rr
        f=F1*dd+F2*dda;ulqr=-K@x
        if kind==0 or kind==7:
            core=np.zeros(4);core[0]=-150*state[0]-550*state[nd]
        elif kind==1:core=-K@state
        elif kind==2:core=ulqr.copy()
        elif kind==3:
            index=min(int(t/.1),len(oracle)-2);frac=(t-index*.1)/.1
            core=ulqr+(1-frac)*oracle[index]+frac*oracle[index+1]
        else:
            xx=features(Tm,tau,direction,eta,flux)
            if kind==6:pred=polynomial(xx,coef,powers)
            elif kind==8:
                distances=np.sum((w0-xx)**2,axis=1)
                basis=np.concatenate((np.ones(1),xx,np.exp(-b0[0]*distances)))
                pred=basis@w1
            else:
                h1=xx@w0+b0;h1=np.tanh(h1) if activation==0 else np.maximum(h1,0.)
                h2=h1@w1+b1;h2=np.tanh(h2) if activation==0 else np.maximum(h2,0.)
                pred=(h2@w2+b2)*yscale
            actv=np.sum(td*td)+(.3*dc)**2;gate=actv/(actv+gate_eps) if gate_on else 1.
            core=ulqr+gate*pred
            if learn_bias!=0:
                for k in range(1,4):core[k]+=learn_bias*np.sin(.7*t)
        corebuf[j%len(corebuf)]=core
        auxiliary=np.zeros(4) if kind==0 else -portgain@y
        if noise.size>=20 and kind!=0:
            auxiliary[1:]+=noise[14:17]*normals[ix,7:10]+noise[17:20]
        resbuf[j%len(resbuf)]=auxiliary
        delayedcore=corebuf[(j-core_delay)%len(corebuf)]
        delayedres=resbuf[(j-residual_delay)%len(resbuf)]
        candidate=delayedcore+delayedres
        fe=f+B@delayedres;a=2*B.T@P@x
        bound=-.25*(x@W@x)+2*abs(x@P@fe)-2*x@P@(A@x+f)
        act=0.;slack=0.
        if filter_on:u,act,slack=projection(candidate,a,bound,ulqr+delayedres,limits,rd)
        else:u=np.minimum(np.maximum(candidate,-limits),limits)
        correction=np.sqrt(np.sum(((u-np.clip(candidate,-limits,limits))/limits)**2))
        active+=act;infeasible+=slack>1e-8;sumcor+=correction**2;maxcor=max(maxcor,correction)
        # First-order driver state; ideal driver is alpha=1.
        uactual+=driver_alpha*(u-uactual)
        u=uactual.copy()
        vibration=tip@(y[:d]-rqplant*delta)[1:]
        truef=F1*dd+F2*dda
        margin=2*xt@P@(A@xt+B@u+truef)+.25*(xt@W@xt)-2*abs(xt@P@(truef+B@delayedres))
        running=xt@Q@xt+u@R@u
        sumtip+=vibration*vibration;sumangle+=y[0]*y[0];sumu+=u*u;peaku=np.maximum(peaku,np.abs(u))
        hit=np.abs(u)>=limits*(1-1e-9);sat+=hit;clip+=np.abs(candidate)>limits+1e-10
        runs=(runs+dt)*hit;maxrun=np.maximum(maxrun,runs)
        maxmargin=max(maxmargin,margin);trueviol+=margin>1e-8
        trueerrmax=max(trueerrmax,np.max(np.abs((state-yc)/sc)))
        if j%stride==0:
            rec[nr,0]=t;rec[nr,1:3]=T;rec[nr,3]=delta
            rec[nr,4]=tip@rqplant[1:]*delta;rec[nr,5]=tip@y[1:d];rec[nr,6]=vibration
            rec[nr,7]=y[0];rec[nr,8]=y[d];rec[nr,9:13]=u
            rec[nr,13]=running;rec[nr,14]=act;rec[nr,15]=slack;rec[nr,16]=margin;rec[nr,17]=np.sqrt(np.sum((xt/sc)**2));nr+=1
        if not np.isfinite(running) or running>1e12 or abs(vibration)>.5:
            failure=1;break
        if j<steps:
            cost+=dt*running;Tn=heat_next(T,tau,direction,eta,flux,dt,ptrue)
            y=Ad@y+Bd@u+Ed*((delta+Tn[0]-Tn[1])/2)
            if observer:zhat=obsAd@zhat+obsBd@u+obsEd*((delta+Tn[0]-Tn[1])/2)
            T=Tn
    count=j+1
    totals=np.array([cost,np.sqrt(sumtip/count),np.sqrt(sumangle/count)*206264.806,
        active/count,infeasible,maxmargin,trueviol/count,sumcor**.5/count**.5,maxcor,trueerrmax,float(failure),j*dt])
    channel=np.vstack((peaku,np.sqrt(sumu/count),sat/count,clip/count,maxrun))
    return rec[:nr],totals,channel
