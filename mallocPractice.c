#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<malloc.h>
void calc_P(double *,double *, double **, double *,int, int);
double calc_alpha(double *, double **, int, int);
int main(void)
{
	int n=13, m_jie=0;
	double **P=NULL,*alpha=NULL;
	scanf("%d",&m_jie);
	m_jie=m_jie+1;
	P=(double **)malloc(sizeof(double *)*m_jie);
	alpha=(double *)malloc(sizeof(double)*m_jie);
	if(P==NULL)
		exit(1);
	
	for(int i=0;i<m_jie;i++)
	{
		P[i]=(double *)malloc(sizeof(double)*n);
		if(P[i]==NULL)
			exit(1);
	}

	for(int i=0;i<m_jie;i++)
		for(int j=0;j<n;j++)
		{
			P[i][j]=1;
			printf("P[%d][%d]==%lf\n",i,j,P[i][j]);
		}

	for(int m=1;;m++)
	{/**计算alpha**/
		alpha[m-1] = calc_alpha(X, P, m-1, n);
	//	printf("<debug--main>  p_alpha[%d]=%lf\n", m-1, alpha[m-1]);
	}
	for(int i=0;i<m_jie;i++)
	{
		printf("<debug>--%d\n",i);
		free(P[i]);
	}
		free(P);
	P=NULL;
	return 0;	
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
	for (int j = 0; j < N; j++)
	{
		p_p_alpha[j] = p_p[flag_alpha][j];
		//printf("<debug--alpha>  p_p[%d][%d]==%lf\n", flag_alpha, j, p_p[flag_alpha][j]);
	}
	/**计算α[i]**/
	for (int j = 0; j < N; j++)
	{
		sum_up = sum_up + p_x[j] * p_p_alpha[j] * p_p_alpha[j];
		sum_down = sum_down + p_p_alpha[j] * p_p_alpha[j];
	}
	alpha = sum_up / sum_down;
	//printf("<debug--alpha> alpha=%lf\n", alpha);
	//flag_alpha++;
	return alpha;
}
