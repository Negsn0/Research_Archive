#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int Z,N,A,n1,n2;
double ll,eigen,hbarm,hh,rmin,rmax,imax,im,hbarc,mcc,hbaro,alpha,r0,R,a,u0,uls,V0,aa;
double fac(int n);
double Gamma(double a);
double hyper(int n,double gamma,double x);
double WaveFunc(int n,double x);
double daikei(double (*func)(),double a,double b,int n);
double naiseki(double x);
double sympson( double (*func)(),double a,double b,int n);
double pot(double r);
double harm(double r);

int main(void){
    int i,n,count;
    double x,s,ss;
    ll=2;
    
    /* 
        n1,n2はfor文で行列要素にするためあとでforforに回す。 
        l,mについて言及しない理由は単にある状態(|\qsi>)の波動関数を与えているため、
        lは決まっており、球対称ポテンシャルであるから等方性を持つため。
    */
    n1=2,n2=2;
    
    rmax=15.,rmin=0.001,imax=400.,hh=(rmax-rmin)/imax;
    im=100;
    
    hbaro=41./pow(A,1/3),mcc=939,hbarc=197.33,hbarm=2*mcc/(hbarc*hbarc);
    alpha=sqrt(mcc*hbaro)/hbarc;
    
    /* WSポテンシャルの定数 */
    A=100;
    V0=51.;
    r0=1.27,aa=0.67,R=r0*pow(A,1./3.);

    for(int k=1;k<10;k++){
        double aaaaaaa=k+0.01;
        s=daikei(naiseki,0,10,10000);
        printf("%15.8e\n",s);
    }
    return 0;
}

double hyper(int n,double gamma,double x){
    int s,nn;
    double y0,y1,y2;
    y0=1.;
    y1=1.-x/gamma;
    nn=-n;
    if(n==0){
        return y0;
    }
    for(s=2;s<=nn;s++){
        y2=((2.*(s-1)+gamma-x)*y1-(s-1)*y0)/(s-1+gamma);
        y0=y1;
        y1=y2;
    }
    return y1;
}

double fac(int n){
    int k;
    double z1,z2;
    z1=1;
    if(n==0)
        return z1;
    for(k=1;k<=n;k++){
        z2=z1*k;
        z1=z2;
    }
    return z1;
}

double Gamma(double a){
    int k;
    double k1,aa,b;
    b=modf(a,&aa);
    if(b<0.5){
    return fac(aa-1);
    }else{
        return fac(2*aa)/(pow(2,2*aa)*fac(aa))*pow(M_PI,0.5);
    }
}

double WaveFunc(int n,double x){
    double Nnl,coef;
    Nnl=sqrt((2*Gamma(ll+1.5+n)/(fac(n)*pow(Gamma(ll+1.5),2))));
    coef=sqrt(alpha)*Nnl;
    return coef*pow(x,ll+1)*exp(-0.5*x*x)*hyper(-n,ll+1.5,x*x);
}

double daikei(double (*func)(),double a,double b,int n){
    int i;
    double dx,s;
    dx=(b-a)/n;
    s=0.5*func(a);
    for(i=1;i<n;i++)
        s+=func(a+dx*i);
    s+=0.5*func(b);
    return s*dx;
}

double sympson( double (*func)(),double a,double b,int n){
    int i;
    double dx,s;
    dx=(b-a)/n;
    s=func(a);
    for(i=1;i<n;i++)
        if(i%2==0){
            s+=4.*func(a+dx*i);
        }else{
            s+=2.*func(a+dx*i);
        }
    s+=func(b);
    return dx*s/3.;
}

double WS(double r){
    return -1./(1.+exp((r-R)/aa));
}

double harm(double r){
    return -0.5*mcc*hbarm*hbarm/(hbarc*hbarc)*r*r;
}

double pot(double x){
    return WS(x)+harm(x);
}


/*
    行列要素を計算するときのやつ
    演算子の引数はalphaで割る
*/
double naiseki(double x){
    return WaveFunc(n1,x)*pot(x/alpha)*WaveFunc(n2,x);
}