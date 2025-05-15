#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ALPHA_INV 137.035999  // ファイン構造定数の逆数
#define NDIM 120 //n=<10 l=0~5 (6) spin 2
int N_n=10,N_l=6,N_s=2;
int A,Z,N; // 核子数 
int T,l; //陽子、中性子の切り替え 1中性子、-1陽子
//量子数
int global_n1,global_n2,global_l1,global_l2,global_s1,global_s2;
//積分範囲
double r_max;
// 物理定数 
double hbarm,hbarc,mcc,hbaro,alpha; 
// WSポテンシャル定数 
double V0,r0,aa,R; 
// ポテンシャル定数
double u0,uls; 
double hyper(int n,double gamma,double x);
double fac(int n);
double Gamma(double a);
double WaveFunc(double x,int n,int l);
double daikei(double (*func)(double),double a,double b,int n);
double naiseki(double x);
double WSnaiseki(double r);
double LSnaiseki(double r);
double Coulnaiseki(double r);
double H_element(int n1,int l1, int s1,int n2,int l2,int s2);
double Ur(double r);
double WS(double r);
double dWS(double r);
double LS(int l1,int l2,int s1,int s2);
double U_coul(double r);
void yacobi(int n, double aa[NDIM][NDIM],double tt[NDIM][NDIM],double eps);

int main(void){
    //核種
    Z=50;N=50;A=Z+N;
    T=0;
    //物理定数
    mcc=939;hbarc=197.33;hbarm=2*mcc/(hbarc*hbarc);alpha=sqrt(mcc * hbaro) / hbarc;
    r_max=50.0;
    //核種ごとに異なる定数
    u0=-51+33*(N-Z)/A;uls=22-14*(N-Z)/A;
    //WS定数
    r0 = 1.27;
    aa = 0.67;
    R = r0 * pow(A, 1. / 3.);

    double H[NDIM][NDIM],tt[NDIM][NDIM];

    for(int i=0;i<NDIM;i++){
        for(int j=0;j<NDIM;j++){
            H[i][j]=0.0;
        }
    }

    for(int n1=0;n1<N_n;n1++){
        for(int l1=0;l1<N_l;l1++){
            for(int s1=0;s1<N_s;s1++){
                int i =((n1*N_l)+l1)*N_s+s1;
                for(int n2=0;n2<N_n;n2++){
                    for(int l2=0;l2<N_l;l2++){
                        for(int s2=0;s2<N_s;s2++){
                            int j =((n2*N_l)+l2)*N_s+s2;
                            H[i][j]=H_element(n1,l1,s1,n2,l2,s2);
                            if(i!=j){
                                H[j][i]=H[i][j];
                            }
                        }
                    }
                }
            }
        }
    }
    yacobi(10,H,tt,1.e-8);

    double eigenvalues[NDIM];
    for(int i=0; i<NDIM; i++){
        eigenvalues[i] = H[i][i];
    }

    // 簡単なバブルソート（効率は良くないが例示目的）
    for(int i=0; i<NDIM-1; i++){
        for(int j=0; j<NDIM-i-1; j++){
            if(eigenvalues[j] > eigenvalues[j+1]){
                double temp = eigenvalues[j];
                eigenvalues[j] = eigenvalues[j+1];
                eigenvalues[j+1] = temp;
            }
        }
    }

    // ソート後のエネルギー固有値の出力
    printf("\nSorted Energy Eigenvalues:\n");
    for(int i=0; i<NDIM; i++){
        printf("E_sorted[%d] = %lf\n", i, eigenvalues[i]);
    }
    return 0;

}

//ポテンシャル系
double WS(double r){
    return -1./(1+exp((r-R)/aa));
}

double dWS(double r){
    return (exp((r-R)/aa))/(aa*pow((1+exp((r-R)/aa)),2));
}
//s1やs2は0,1で与えられる。
double LS(int l1,int l2,int s1,int s2){
    if(s1==s2&&l1==l2){
        int l=l1;
        double s=-0.5+s1;
        return 0.5*((l+s)*(l+s+1)-l*(l+1)-s*(s+1));
    }
    return 0.0;
}

double U_coul(double r){
    if (r<=R){
        return 0.5*(Z-1)*(3*R*R-r*r)/(pow(R,3.)*ALPHA_INV);
    }else{
        return (Z-1)/(r*ALPHA_INV);
    }
}

//全ポテンシャル
double Ur(double r){
    double pot=u0*WS(r)+uls*r0*r0*dWS(r)*LS(global_l1,global_l2,global_s1,global_s2)/r-0.5*mcc*pow(hbaro/hbarc,2)*r*r;
    if(T==-1){
        pot += U_coul(r);
    }
    return pot;
}


//一体ポテンシャルによるエネルギー固有値
double E0(int n,int l){
    return hbaro*(2*n+l+1.5);
}

double WSnaiseki(double r){
    return WaveFunc(r,global_n1,global_l1)*u0*WS(r)*WaveFunc(r,global_n2,global_l2);
}

double LSnaiseki(double r){
    return WaveFunc(r,global_n1,global_l1)*(dWS(r)/r)*WaveFunc(r,global_n2,global_l2);
}

double Coulnaiseki(double r){
    return WaveFunc(r,global_n1,global_l1)*U_coul(r)*WaveFunc(r,global_n2,global_l2);
}

//見直してほしい
/*
double H_element(int n1,int l1,int s1,int n2,int l2,int s2){
    global_n1=n1,global_l1=l1,global_s1=s1;
    global_n2=n2,global_l2=l2,global_s2=s2;
    double val=0;
    if (n1==n2&&l1==l2){
        return E0(n1,l1)+daikei(naiseki,0.0,r_max,1000);
    }else{
        return daikei(naiseki,0.0,r_max,1000);
    }
}

*/
double H_element(int n1,int l1,int s1,int n2,int l2,int s2){
    global_n1=n1,global_l1=l1,global_s1=s1;
    global_n2=n2,global_l2=l2,global_s2=s2;
    double val=0.0;
    //HO対角成分
    if (n1 == n2&&l1 == l2&&s1 == s2){
        val += E0(n1,l1);
    }
    //WS(球対称)pt
    if (l1 == l2&&s1 == s2){
        val+=daikei(WSnaiseki,0.0,50,2000);
    }
    //LSポテンシャル
    if (l1 == l2&&s1 == s2){
        double LSfac=uls*r0*r0*LS(global_l1,global_l2,global_s1,global_s2);
        val+=LSfac*daikei(LSnaiseki,0.0,50,2000);
    }
    //Coulomb
    if (T==-1&&l1 == l2&&s1 == s2){
        val+=daikei(Coulnaiseki,0.0,50,2000);
    }

    return val;
}

//行列要素周辺
double daikei(double (*func)(double),double a,double b,int n){
    int i;
    double dx,s;
    dx=(b-a)/n;
    s=0.5*func(a);
    for(i=1;i<n;i++)
        s+=func(a+dx*i);
    s+=0.5*func(b);
    return s*dx;
}

double naiseki(double x){
    return WaveFunc(x,global_n1,global_l1)*Ur(x)*WaveFunc(x,global_n2,global_l2);
}


//波動関数構成要素
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

double WaveFunc(double x,int n,int l){
    double Nnl,coef;
    Nnl=sqrt((2*Gamma(l+1.5+n)/(fac(n)*pow(Gamma(l+1.5),2))));
    coef=sqrt(alpha)*Nnl;
    return coef*pow(alpha*x,l+1)*exp(-0.5*alpha*alpha*x*x)*hyper(-n,l+1.5,alpha*alpha*x*x);
}

void yacobi(int n, double aa[NDIM][NDIM],double tt[NDIM][NDIM],double eps){
    int p,q,i,j,k;
    /*ttの初期化*/
    for(i=0;i<NDIM;i++){
        for(j=0;j<NDIM;j++){
            if (i==j)
            {
                tt[i][j]=1.0;
            }else
            {
                tt[i][j]=0.0;
            }
            
            
        }
    }
    int count = 0;
    double aaa[NDIM][NDIM],ttt[NDIM][NDIM];
    double maxvalue = 1.0;
    double phi,K,c,s;
    while (count < 10000 && maxvalue > eps){
        /*aaの非対角成分で絶対値が最も大きいaa[i][j]を検索し、maxvalueに代入*/
        maxvalue=0;
        for(i=0;i<NDIM;i++){
            for(j=0;j<NDIM;j++){
                if(i!=j&&maxvalue<fabs(aa[i][j])){
                    maxvalue=fabs(aa[i][j]);
                    p=i;
                    q=j;
                }
            }
        }
        /*phiを求める*/
        K=2*aa[p][q]/(aa[p][p]-aa[q][q]);
        phi=0.5*atan(K);
        c=cos(phi);
        s=sin(phi);
        /*aa,ttを対角化するのに一旦aaa,tttに計算結果を保持させる*/
        for(i=0;i<NDIM;i++){
            for(j=0;j<NDIM;j++){
                aaa[i][j]=aa[i][j];
                ttt[i][j]=tt[i][j];
                
            }
        }
        for(i=0;i<NDIM;i++){
            aaa[p][i]=aa[p][i]*c+aa[q][i]*s;
            aaa[i][p]=aa[p][i]*c+aa[q][i]*s;
            aaa[q][i]=-aa[p][i]*s+aa[q][i]*c;
            aaa[i][q]=-aa[p][i]*s+aa[q][i]*c;
            ttt[i][p]=tt[i][p]*c-tt[i][q]*s;
            ttt[i][q]=tt[i][p]*s+tt[i][q]*c;
        }
        aaa[p][p]=aa[p][p]*c*c+aa[q][q]*s*s+aa[p][q]*sin(2*phi);
        aaa[q][q]=aa[p][p]*s*s+aa[q][q]*c*c-aa[p][q]*sin(2*phi);
        aaa[p][q]=-0.5*(aa[p][p]-aa[q][q])*sin(2*phi)+aa[p][q]*cos(2*phi);
        aaa[q][p]=-0.5*(aa[p][p]-aa[q][q])*sin(2*phi)+aa[p][q]*cos(2*phi);

        /*aaa,tttの中身をaa,ttに移し替える*/
        for(i=0;i<NDIM;i++){
            for(j=0;j<NDIM;j++){
                aa[i][j]=aaa[i][j];
                tt[i][j]=ttt[i][j];
            }
        }
        count++;
    }
}
