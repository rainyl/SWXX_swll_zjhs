#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
void input(double *,double *, int);
void sort(double * ,double *, int);//排序函数
void calc_xy(double *, double *, int);//
double calc_alpha(double *, double **, int, int);//
int main(void)
{
	/**变量的定义**/
	int m_jie = 1;//拟合阶数a
	int m=1;
	int n = 0, c = 0;
	double **P;
	double *Z, *Q,  *X, *Y,  *alpha ;
	FILE *fp_data=NULL;
	if((fp_data=fopen("./data.txt","r"))==NULL)
	{
		printf("打开文件“data.txt”失败！\n");
	}
	/**获取行数**/
	while(((c=fgetc(fp_data))!=EOF))
	{
		if(c=='\n')
			n++;
	}
	printf("n==%d\n",n);
	Z=(double *)malloc(sizeof(double)*n);
	Q=(double *)malloc(sizeof(double)*n);
	rewind(fp_data);
	for(int i=0;i<n;i++)
	{
		fscanf(fp_data,"%lf\t%lf",&Z[i],&Q[i]);
		printf("Z==%lf\tQ==%lf\n",Z[i],Q[i]);
	}
	/*输入拟合的阶数*/
	printf("请输入拟合的阶数:\n");
	scanf("%d",&m_jie);
	m_jie=m_jie+1;
	P=(double **)malloc(sizeof(double *)*m_jie);
	for(int i=0;i<m_jie;i++)
		P[i]=(double *)malloc(sizeof(double)*n);
	/**初始化二维数组P[][]**/
	for (int i = 0; i < m_jie; i++)
	{
		for (int j = 0; j < n; j++)
		{
			P[i][j] = 1;
			printf("p[%d][%d]==%lf\n",i,j,P[i][j]);
		}
	}
	/**排序**/
	sort(Z,Q,n);
	/**将Z，Q的值赋给X，Y**/
	X=(double *)malloc(sizeof(double)*n);
	Y=(double *)malloc(sizeof(double)*n);
	for (int i = 0; i < n; i++)
	{
		X[i] = Z[i];
		Y[i] = Q[i];
	}
	free(Z);
	free(Q);
	/**计算X，Y**/
	calc_xy(X, Y, n);
	alpha=(double *)malloc(sizeof(double)*m_jie);
	
		/**计算alpha**/
		alpha[m-1] = calc_alpha(X, P, m-1, n);
		printf("<debug--main>  p_alpha[%d]=%lf\n", m-1, alpha[m-1]);
		
	//system("pause");
	return 0;
}
void sort(double *p_z, double *p_q, int N)
{
	printf("<debug--sort> enter sort\n");
	double tmp = 0;
	for (int i = 0; i < N - 1 ; i++)
	{
		for (int j = 0; j < N-1-i; j++)
		{
			if (p_z[j]>p_z[j+1])
			{
				tmp = p_z[j];
				p_z[j] = p_z[j+1];
				p_z[j+1] = tmp;
				tmp = p_q[j];
				p_q[j] = p_q[j + 1];
				p_q[j + 1] = tmp;
			}
		}
	}
	printf("<debug--sort> leave  sort\n");
}
/**计算X，Y**/
void calc_xy(double *p_x, double *p_y, int N)
{
	double Z_0 = 2.3;
	//计算X=Ln[Z-Z_0]，Y=Ln[Q]
	for (int i = 0; i < N; i++)
	{
		p_x[i] = log(p_x[i] - Z_0);
		//printf("<debug>X[%d]==%lf\n", i, p_x[i]);
		p_y[i] = log(p_y[i]);
		//printf("<debug>Y[%d]==%lf\n", i, p_y[i]);
	}

}

double calc_alpha(double *p_x, double **p_p, int flag_alpha, int N)
{
	/**
	 * 	此函数计算alpha并返回alpha
	 * 	sum_up表示计算式分子，p_p_alpha[N]用来在不影响原始p字符串的情况下拷贝用于计算的一行P
	 * 	flag_alpha用来控制上述的一行p_p_alpha[N],再每次调用此函数后增加一
	**/
	double alpha = 0, sum_up = 0, sum_down = 0, *p_p_alpha;
	//static int flag_alpha = 0;     //改进： 传递参数flag_alpha,由上层函数控制
	p_p_alpha=(double *)malloc(sizeof(double)*N);
	//for (int j = 0; j < N; j++)
//.cache	{
//		p_p_alpha[j] = p_p[flag_alpha][j];
//		//printf("<debug--alpha>  p_p[%d][%d]==%lf\n", flag_alpha, j, p_p[flag_alpha][j]);
//	}
	/**计算α[i]**/
	for (int j = 0; j < N; j++)
	{
		sum_up = sum_up + p_x[j] * p_p[flag_alpha][j] * p_p[flag_alpha][j];
		sum_down = sum_down + p_p[flag_alpha][j] * p_p[flag_alpha][j];
	}
	alpha = sum_up / sum_down;
	//printf("<debug--alpha> alpha=%lf\n", alpha);
	//flag_alpha++;
	return alpha;
}
