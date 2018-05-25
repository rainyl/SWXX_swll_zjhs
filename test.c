#include <stdio.h>
#include <string.h>
#include <stdlib.h>
int main(void){
    int d=0;
    char c='0';
    FILE *fp;
    if ((fp=fopen("test.dat","wb+"))==NULL)
    {
        printf("open error!\n");
        exit(1);
    }
    for (size_t i = 0; i < 5; i++)
    {
        scanf("%d%c", &d, &c);
    fprintf(fp, "%d\t%c\n", d, c);
    }
    
    fclose(fp);
    return 0;
}