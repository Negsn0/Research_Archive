#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ALPHA_INV 137.035999  // ファイン構造定数の逆数

int A,Z,N; // 核子数 
int T,l; //陽子、中性子の切り替え 1中性子、-1陽子
//量子数
int global_n1,global_n2,global_l1,global_l2;
double s;
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
double H_element(int n1,int l1,int n2,int l2);
double Ux(double x);
double WS(double x);
double dWS(double x);
double LS(int l,double s);
double U_coul(double x);


int main(void) {
    // 計算対象の核種の質量数を設定
    Z=50;
    N=50;
    A = Z + N;
    // 初期化
    FILE* output = fopen("sn_levels_single2.txt", "w");
    if (output == NULL) {
        printf("ファイルを開けませんでした。\n");
        return 1;
    }

    hbaro = 41. / pow(A, 1. / 3.); 
    mcc = 939;
    hbarc = 197.33;
    hbarm = 2 * mcc / (hbarc * hbarc);
    alpha = sqrt(mcc * hbaro) / hbarc;
    r_max = 50.0;

    // ウッズ・サクソンポテンシャル定数
    V0=50.;
    u0 = -51 + 33 * (N - Z) / A;
    uls = 22 - 14 * (N - Z) / A;
    r0 = 1.27;
    aa = 0.67;
    R = r0 * pow(A, 1. / 3.);
    // ヘッダーを出力
    fprintf(output, "A, Z, N, Particle, (n, l, s), Energy (MeV)\n");
    printf("Calculating for Nucleus: A=%d, Z=%d, N=%d\n", A, Z, N);

    // 陽子の準位計算
    T = -1;  // 陽子
    for (int n = 0; n <= 2; n++) {  // 主量子数 n の範囲
        for (int l = 0; l <= n + 1; l++) {  // 軌道角運動量 l の範囲
            for (int s = -1; s <= 1; s += 2) {  // スピン s = ±1/2
                double energy = H_element(n, l, n, l);
                fprintf(output, "%d, %d, %d, Proton, (%d, %d, %+.1f), %.6f\n",
                        A, Z, N, n, l, s * 0.5, energy);
            }
        }
    }

    // 中性子の準位計算
    T = 1;  // 中性子
    for (int n = 0; n <= 3; n++) {  // 主量子数 n の範囲
        for (int l = 0; l <= n + 2; l++) {  // 軌道角運動量 l の範囲
            for (int s = -1; s <= 1; s += 2) {  // スピン s = ±1/2
                double energy = H_element(n, l, n, l);
                fprintf(output, "%d, %d, %d, Neutron, (%d, %d, %+.1f), %.6f\n",
                        A, Z, N, n, l, s * 0.5, energy);
            }
        }
    }

    fclose(output);
    printf("計算が完了しました。結果は 'sn_levels_single.txt' に保存されました。\n");
    return 0;
}


// WSポテンシャル周辺とその微分
double WS(double x){
    return -V0/(1+exp((x-R)/aa));
}

double dWS(double x){
    return -V0*(exp((x-R)/aa))/(aa*pow((1+exp((x-R)/aa)),2));
}

// 軌道角運動量l,スピン角運動量s,全角運動量j
double LS(int ll,double ss){
    double j=ll+0.5*ss;
    return 0.5*(j*(j+1)-ll*(ll+1)-0.75);
}

//Coul.ポテンシャル
double U_coul(double x){
    if (x>R){
        return hbarc*(Z-1)/(x*ALPHA_INV);
    }else{
        return hbarc*(Z-1)*0.5*(3.-(pow(x/R,2)))/(ALPHA_INV*R);
    }
}
// ポテンシャル 
double Ux(double x){
    double pot= u0*WS(x)+uls*r0*r0*dWS(x)*LS(global_l1,s) -0.5*mcc*pow(hbaro/hbarc,2)*x*x;
    if (T==-1){
        pot += U_coul(x);
    }
    return pot;
}

//一体ポテンシャルによるエネルギー固有値
double E0(int n,int l){
    return hbaro*(2*n+l+1.5);
}

/* 波動関数と行列要素を求めるためのもの */
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
    return WaveFunc(x,global_n1,global_l1)*Ux(x)*WaveFunc(x,global_n2,global_l2);
}

double H_element(int n1,int l1,int n2,int l2){
    global_n1=n1;
    global_l1=l1;
    global_n2=n2;
    global_l2=l2;
    if (n1==n2&&l1==l2){
        return E(n1,l1)+daikei(naiseki,0,r_max,10000);
    }else{
        return daikei(naiseki,0,r_max,10000);
    }
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
