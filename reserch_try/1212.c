#include <stdio.h>
#include <stdlib.h>
int l;
int main(void){
    int i;
    for(i=0;i<5;i++){
        l=i;
        printf("%lf\n",l+1.5);
    }
    return 0;
}