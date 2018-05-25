#include <stdio.h>
#include <string.h>
#include <stdlib.h>
int main(void)
{
	FILE *fp_data=NULL;
	double Z[13]={0},Q[13]={0};
	if((fp_data=fopen("./data.txt","r"))==NULL)
	{
		printf("打开文件“data.txt”失败！\n");
		exit(0);
	}
	for(int i=0;i<13;i++)
	{
		fscanf(fp_data,"%lf,%lf",&Z[i],&Q[i]);
		printf("Z==%lf,Q==%lf\n",Z[i],Q[i]);
	}
	rewind(fp_data);
	//fseek(fp_data,0,SEEK_END);
	char c='0';
	int N=0;
	while(((c=fgetc(fp_data))!=EOF))
	{
		printf("enter read \n");
		if(c=='\n')
		{
			printf("%d\n",c);
			N++;
		}
	}
	printf("N==%d\n",N);
	fclose(fp_data);
	/**
	double Z[13] = {2.71, 2.64, 2.88, 2.94, 3,    3.21, 3.1,  3.34, 3.47, 3.61, 3.52, 3.76, 3.87};
	double Q[13] = {4.47, 2.19, 11,   14.3, 18.5, 32.1, 24.2, 43.5, 62.4, 86.7, 71.4, 119,   148};
	
	for(int i=0;!feof(fp_data);i++)	
		fscanf(fp_data,"%lf\t%lf",&Z[i],&Q[i]);
**/
		return 0;
}
