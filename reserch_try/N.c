#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int Ndeg(int Nmax);
int main(void){
    printf("%d\n",Ndeg(8));
    return 0;
}

/* ブロック行列の次元数を求める */
int Ndim(int N){
    int l,i,count;
    count=0;
    if(N%2==0){
        for(l=0;l<=N;l+=2){
            ++count;
        }
    }else{
        for(l=1;l<=N;l+=2){
            ++count;
        }
    }
    return count;
}


