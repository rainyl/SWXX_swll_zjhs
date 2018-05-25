#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
//#define N 13
void input(double *,double *, int);
void sort(double * ,double *, double *, int);//排序函数
/**
 * **********计算函数*******************
 * calc_xy    --->   计算X，Y
 * calc_P     --->   计算P
 * calc_apb   --->   计算alpha,P,beta
 * calc_alpha --->   计算alpha
 * calc_beta  --->   计算beta
 * calc_Y_c   --->   计算最终的水位Y_c
 * calc_A     --->   计算系数A[i]
 * calc_Q_c   --->   计算最终的流量Qc
 * *************************************
**/
void calc_xy(double *, double *, int, double);//
double calc_P(double *,double *, double **, double *,int, int);//
double calc_alpha(double *, double **, int, int);//
double calc_beta(double **, int,  int );//
double calc_Y_c(double **, double *,int,int);//
double calc_A(double **, double *, int, int);//
double calc_Q_c(double *, int );
double calc_Q_unknow(double, int, double *, double *, double *);
/**
 * *******检验函数*********
 * FH    -->   符号检验
 * SX    -->   适线检验
 * PLSZ  -->   偏离数值检验
 * T     -->   t(学生氏)检验
 * ***********************
**/
int chk_FH(double *, double *, int);
int chk_SX();
int chk_PLSZ();
int chk_T();//对于本题目，由于没有历年数据，因此无法应用t（学生式）检验
int main(void)
{ 
	/**变量的定义**/
	int res_FH = 0, res_SX = 0, res_PLSZ = 0;
	int m_jie = 1;//拟合阶数
	int isContinue=0;//判断是否继续拟合
	int n = 0, c = 0;
	char fileName[10];
	double **P;
	double *Z=NULL, *Q=NULL, *area=NULL,  *X=NULL, *Y=NULL, *Y_c=NULL, *Qc=NULL, *alpha=NULL, *beta=NULL, *A=NULL ;
	double Q_unknow=0;
	FILE *fp_data=NULL;
	printf("\n本程序使用前请将数据以文本文件的形式存储，水位与流量之间以分页符（tab）隔开，一行一组数据\n\n");
	printf("*********请输入待处理数据的文件名：");
	scanf("%s",fileName);
	if((fp_data=fopen(fileName,"r"))==NULL)
	{
		printf("打开文件“%s”失败！\n",fileName);
		exit(0);
	}else
	{
		printf("\n\n×××××××××文件打开成功，即将开始计算×××××××××\n\n");
		/**获取行数**/
		while(((c=fgetc(fp_data))!=EOF))
		{
			if(c=='\n')
				n++;
		}
		Z=(double *)malloc(sizeof(double)*n);
		Q=(double *)malloc(sizeof(double)*n);
		area=(double *)malloc(sizeof(double)*n);
		rewind(fp_data);
		/**读入数据**/
		for(int i=0;i<n;i++)
		{
			fscanf(fp_data,"%lf\t%lf\t%lf",&Z[i],&Q[i],&area[i]);
		}
	}
	/**输入函数**/
	//input(Z,Q);
	/**排序**/
	printf("×××××××××开始对水位Z与对应的流量Q进行排序×××××××××\n");
	sort(Z,Q,area,n);
	//** debug,验证排序是否成功
	printf("×××××××××排序完成，计算平均流速结果如下×××××××××\n");
	printf("---N---|----Z----|----Q----|----A----|----U----\n");
	for (int i = 0; i < n; i++)
	{
		printf("%7d|", i);
		printf("%9.2lf|", Z[i]);
		printf("%9.2lf|", Q[i]);
		printf("%9.2lf|", area[i]);
		printf("%9.2lf\n", Q[i]/area[i]);
	}
//	**/
	/**将Z，Q的值赋给X，Y**/
	X=(double *)malloc(sizeof(double)*n);
	Y=(double *)malloc(sizeof(double)*n);
	for (int i = 0; i < n; i++)
	{
		X[i] = Z[i];
		Y[i] = Q[i];
	}	
	double Z_0 = 0;
	printf("*********请输入Z_0:");
	scanf("%lf",&Z_0);

	/**计算X，Y**/
	calc_xy(X, Y, n, Z_0);
	/**
	for (int i = 0; i < n; i++)
	{
		printf("X[%d] = %lf\t",i,X[i]);
		printf("Y[%d] = %lf\n",i,Y[i]);
	}
	**/
	printf("默认最大拟合阶数为5阶，是否修改？（修改请输入1，否则请输入0）  ");
	scanf("%d",&m_jie);
	if(m_jie==0)
		m_jie=5;
	else
	{
		printf("请输入想要拟合的阶数:");
		scanf("%d",&m_jie);
	}
	//m_jie++;
	P=(double **)malloc(sizeof(double *)*(m_jie+1));
	for(int i=0;i<m_jie+1;i++)
		P[i]=(double *)malloc(sizeof(double)*n);
	/**初始化二维数组P[][]**/
	for (int i = 0; i < m_jie+1; i++)
	{ 
		for (int j = 0; j < n; j++)
		{ 
			P[i][j] = 1;
		}
	}


	alpha=(double *)malloc(sizeof(double)*(m_jie+1));
	beta=(double *)malloc(sizeof(double)*m_jie);
	for(int i=0;i<m_jie;i++)
	{
		beta[i]=0;
	}
	A=(double *)malloc(sizeof(double)*(m_jie+1));
	Y_c=(double *)malloc(sizeof(double)*n);
	Qc=(double *)malloc(sizeof(double)*n);
//	m_jie--;
	for(int m=1;;m++)
	{	     
		/**计算alpha**/
		alpha[m-1] = calc_alpha(X, P, m-1, n);
		/**计算P[][]**/
		for(int j=0;j<n;j++)
		{
			P[m][j] = calc_P(X, alpha, P, beta,m, j);
	//		printf("<<debug--main>>P[%d][%d]==%lf\n",m,j,P[m][j]);
		}
		/**计算beta**/
		beta[m] = calc_beta(P, m,  n);
	//	printf("<debug--main> beta=%lf\n", beta[m]);
		/**计算系数A[]**/
		if(m==1)
		{ 
			A[m-1]=calc_A(P, Y, m-1, n);
			A[m]=calc_A(P,Y,m, n);
		}else
		{
			A[m]=calc_A(P, Y, m, n);
		}
		printf("×××××××××【%d】阶拟合曲线  Y=",m);
		for(int j=0;j<=m;j++)
		{ 
			printf("%.2lf *  ", A[j]);
			printf("P[%d] + ", j);
		}
		printf("\n\n");
		/**计算拟合Yc**/
 		for (int i = 0; i < n; i++)
		{ 
			Y_c[i] = calc_Y_c(P, A,i,m);
		//	printf("<debug--main>  Y_c[%d]=%lf\n", i, Y_c[i]);
		}	
		/*计算拟合流量Qc*/
 		for (int i = 0; i < n; i++)
		{  
			Qc[i] = calc_Q_c(Y_c,i);
			//printf("Qc[%d]=%lf\n", i, Qc[i]);
		}

		printf("×××××××××实测水位、实测流量、拟合流量分别如下×××××××××\n");
		printf("---N---|----Z----|----Q----|----Qc----\n");
		for(int i=0;i<n;i++)
		{
			printf("%7d|", i);
			printf("%9.2lf|", Z[i]);
			printf("%9.2lf|", Q[i]);
			printf("%9.2lf\n", Qc[i]);
		}
		/** 
		 * chk_FH()为符号检验 
		 * chk_SX()为适线检验
		 * chk_PLSZ()为偏离数值检验
		**/
		res_FH = chk_FH(Q, Qc,n);
		if (res_FH==0)
 		{
			printf("=========符号检验失败！=========\n\n");
			//continue;
		}
		else
 		{
			printf("=========符号检验通过!=========\n\n");
		}
		res_SX = chk_SX(Q, Qc, n);
		if (res_SX==0)
 		{
			printf("=========适线检验失败！=========\n\n");
		}
		else
 		{
			printf("=========适线检验通过!=========\n\n");
		}
	
		res_PLSZ = chk_PLSZ(Q, Qc, n);
		if (res_PLSZ==0)
		{
			printf("=========偏离数值检验失败！=========\n\n");
		}
		else
		{
			printf("=========偏离数值检验通过!=========\n\n");
		}
		printf("<debug--main>m==%d\tm_jie==%d\tsContinue==%d\t\n",m,m_jie,isContinue);
		if(m<m_jie)
		{
			int flag=0;
			for(int i=0;;i++)
			{
				printf("[%d]阶拟合完成且检验结果如上\n\t[0]退出\n\t[1]继续\n\t[2]已知水位计算拟合流量\n请输入你的选择：",m);
				scanf("%d",&isContinue);
				if(isContinue==0||isContinue==1||isContinue==2)
				{
					break;
				}else
				{
					printf("输入错误！请重新输入！\n");
				}
			}
			if(isContinue==1)
			{
				printf("=======继续集合【%d】阶曲线========\n\n",m+1);
			}else if(isContinue==0)
			{
				printf("【%d】阶拟合结束！程序即将退出！\n",m);
				break;
			}else
			{
				printf("=======已知水位，开始计算拟合流量========\n\n");
				for(int i=0;;i++)
				{

					Q_unknow = calc_Q_unknow(Z_0, m, alpha, beta, A);
					printf("拟合流量 Q_unknow == %.2lf\n计算完成\n\t[0]退出计算并退出程序\n\t[1]退出计算并继续拟合\n\t[2]继续计算\n请选择：",Q_unknow);
					for(int j=0;;j++)
					{
						scanf("%d",&isContinue);
						if(isContinue==0||isContinue==1||isContinue==2)
							break;
						else
							printf("输入错误，请重新输入：");
					}
					if(isContinue==0)
					{	flag=0;
						break;
					}else if(isContinue==1)
					{
						flag=1;
						break;
					}else
						continue;
				}
			if(isContinue==0&&flag==0)
				break;
			}
		}else
		{
			int flag=0;
			for(int i=0;;i++)
			{
				printf("[%d]阶拟合完成且检验结果如上\n\t[0]退出\n\t[1]已知水位计算拟合流量\n请输入你的选择：",m);
				scanf("%d",&isContinue);
				if(isContinue==0)
				{	
					printf("程序结束，即将退出\n");
					break;
				}
				else if(isContinue==1)
				{
					printf("=======已知水位，开始计算拟合流量========\n\n");
					for(int i=0;;i++)
					{
	
						Q_unknow = calc_Q_unknow(Z_0, m, alpha, beta, A);
						printf("拟合流量 Q_unknow == %.2lf\n\n计算完成\n\t[0]退出计算并退出程序\n\t[1]继续计算\n请选择：",Q_unknow);
						for(int j=0;;j++)
						{
							scanf("%d",&isContinue);
							if(isContinue==0||isContinue==1)
								break;
							else
								printf("输入错误，请重新输入：");
						}
						if(isContinue==0)
						{	
							flag=0;
							break;
						}
					}
				}else
					printf("输入错误，请重新输入\n");
				if(isContinue==0&&flag==0)
					break;
			}
			break;
		}

	}
	free(Z);
	Z=NULL;
	free(Q);
	Q=NULL;
	free(Y_c);
	Y_c=NULL;
	free(Qc);
	Qc=NULL;
	//system("pause");
	return 0;
}

/**输入函数**/
void input(double *p_z, double *p_q, int N)
{



	printf("begain input!\n开始水位输入！\n");
	for (int i = 0; i < N;i++)
	{
		printf("请输入第[%d]个水位Z[%d]:", i, i);
		scanf("%lf", p_z++);
		//printf("<debug>x=%d\n", p_x);
	}
	for (int i = 0; i < N;i++)
	{
		printf("请输入第[%d]个流量Q[%d]:", i+1, i+1);
		scanf("%lf", p_q++);
	}
}
/**排序函数，采用冒泡排序法**/
void sort(double *p_z, double *p_q, double *area, int N)
{
	double tmp = 0;
	//printf("\n开始排序\n");
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
				tmp=area[j];
				area[j]=area[j+1];
				area[j+1]=tmp;
			}
		}
	}
	//printf("排序结束\n");
}

/**计算X，Y**/
void calc_xy(double *p_x, double *p_y, int N,double Z_0)
{
	//计算X=Ln[Z-Z_0]，Y=Ln[Q]
	for (int i = 0; i < N; i++)
	{
		p_x[i] = log(p_x[i] - Z_0);
		//printf("<debug>X[%d]==%lf\n", i, p_x[i]);
		p_y[i] = log(p_y[i]);
		//printf("<debug>Y[%d]==%lf\n", i, p_y[i]);
	}

}

/**
 * 计算Pi
 * 需要Xi，alpha，beta以及P[i-1]
**/
double calc_P(double *p_x,double *p_alpha,double **p_p,double *p_beta,int m, int j)
{
	/**
	 * flag_p用来控制p_p（包括alpha数组）
	 * 由于计算P[i]时使用的alpha为alpha[i-1],使用的beta为beta[i-2],因此令flag_p初始值为2
	 * 一次计算一组P
	**/
	//printf("=========开始计算P[%d][%d]=========\n\n",m,j);
	int k=0;
	double res_P=0;
	if(m==1)
		k=m-1;
	else
		k=m-2;
	/**debug
	printf("P[%d][%d]==%lf\n",m,j,p_p[m][j]);
	printf("P[%d][%d]==%lf\n",m-1,j,p_p[m-1][j]);
	printf("alpha[%d]==%lf\n",m-1,p_alpha[m-1]);
	printf("beta[%d]==%lf\n",m-1,p_beta[m-1]);
	printf("P[%d][%d]==%lf\n",k,j,p_p[k][j]);
	**/
	res_P = p_p[m-1][j] * (p_x[j] - p_alpha[m-1]) - p_beta[m-1] * p_p[k][j];
	return res_P;
//	printf("=========P[%d][1~%d]计算结束！=========\n\n\n",flag_p-1,N);
}

/**
 * 计算alpha
 *
**/
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

/**
 * 计算beta   需要P
 * beta=sum(Pi^2)/sum(P[i-1]^2)
**/
double calc_beta(double **p_p, int m, int N)
{
	/**此处用静态变量flag_beta来控制用于计算的两行，每次使用此函数fla_beta增加一
	 * 注：此设计不利于函数的方便使用
	 * 后面改进
	**/
	//static int flag_beta = 2;
	double beta = 0, sum_up = 0, sum_down = 0;
	//p_p_beta=(double *)malloc(sizeof(double)*2);
	//for(int i=0;i<N;i++)
	//	p_p_beta[i]=(double *)malloc(sizeof(double)*N);
	/**
	for (int i = 0; i < m+1; i++)
	{
		for (int j = 0; j < N; j++)
		{
			printf("<debug--beta>  beat_pp[%d][%d]==%lf\n", i, j, p_p[i][j]);
		}
	}
	**/
	/*由于计算beta只需要最新的两行P，因此这里cp最新的两行P
	for (int j = 0; j < N; j++)
	{
		p_p_beta[0][j] = p_p[flag_beta-2][j];
		p_p_beta[1][j] = p_p[flag_beta-1][j];
	}
	*/
	/**正经计算beta的分子与分母**/
	for (int j = 0; j < N; j++)
	{
		sum_up = sum_up + p_p[m][j] * p_p[m][j];
		sum_down = sum_down + p_p[m-1][j] * p_p[m-1][j];
	}
	beta = sum_up / sum_down;
	return beta;
}

/**
 * 计算系数A[0],A[1],A[2]...
 * 需要P以及Y
**/
double calc_A(double **p_p, double *p_y, int flag_a, int N)
{
	double A=0, sum_up = 0, sum_down = 0, *p_p_A;
	p_p_A=(double *)malloc(sizeof(double)*N);
	for (int i = 0; i < N; i++)
	{
		p_p_A[i] = p_p[flag_a][i];
	}
	for (int i = 0; i < N;i++)
	{
	sum_up += p_p_A[i] * p_y[i];
	sum_down += p_p_A[i] * p_p_A[i];
	}
	A = sum_up / sum_down;
	free(p_p_A);
	return A;
}

/** 
 * 计算Yc
 * 用到的数据为Pi，A      公式为Y=A[0]P[0][0]+A[1]P[1][0]+A[2]P[2][0]+...
 * 其中j用来指定  A[j]   P[i][j]
**/
double calc_Y_c(double **p_p, double *p_a, int j,int m)
{
	//printf("<debug--Y_c>  enter Y_c\n");
	double sum_Y_c = 0;
		sum_Y_c = 0;
		for (int i=0 ; i < m + 1; i++)
		{
			//printf("<debug--y_c>  A[%d]=%lf\n", i,p_a[i]);
			//printf("<debug--y_c>  P[%d][%d]=%lf\n",i, j, p_p[i][j]);
			sum_Y_c += p_a[i] * p_p[i][j];
			//printf("<debug--y_c>  sum_Y_c=%lf\n", sum_Y_c);
		}
		return sum_Y_c;
	//printf("<<debug>>Y_c==%lf\n\n",Y_c);
}

/**
 * 计算Qc
 * 用到的数据为Yc，此处传递的i是用来制定计算Q[i]
**/
double calc_Q_c(double *p_Y_c, int i)
{
	double Q_c = 0;
	Q_c = exp(p_Y_c[i]);
	//printf("<debug--Qc>  Q_c=%lf\n", Q_c);
	return Q_c;
}

/**
 * 符号检验
 * 需要Q,Qc
**/
int chk_FH(double *Q, double *Q_c, int N)
{
	int result = 0, flag = 0;
	int accdny = 1;
	float k_0 = 0, k_1 = 0, k_min = 0, XZXSP[3] = {0.05, 0.10, 0.25}, XZXSP_tmp = 0;
	double u = 0, u_std[3] = {1.96, 1.64, 1.15};
	printf("=========开始符号检验=========\n");
	for (int i = 0; i < N; i++)
	{
		if (Q[i] > Q_c[i])
	  	{
			k_0++;
	  	} else if (Q[i] < Q_c[i])
	  	{
      		k_1++;
      	} else
	  	{
      		k_0 += 0.5;
     		k_1 += 0.5;
      	}
  	}
	if (k_0<=k_1)
	{
		k_min = k_0;
	}else
	{
		k_min = k_1;
	}
	printf("---------k_0=%f\tk_1=%f---------\n", k_0,k_1);
	printf("---------默认显著性水平α=0.05，是否接受？(接受请输入1，否则为0):");
	scanf("%d", &accdny);
	if (accdny == 0 )
	{
		printf("---------请输入显著性水平α(0.05or0.10or0.25):---------\n");
		scanf("%f", &XZXSP_tmp);
		for (int i = 0; i < 3; i++)
		{
			if (XZXSP_tmp==XZXSP[i])
			{
				flag = i;
			}
		}
	}
	//printf("<debug--fh>  N==%lf\tk_min==%lf\t\n", (double)N, (double)k_min);
	u = (0.5 * (double)N - (double)k_min - 0.5) / (0.5 * sqrt((double)N));
	printf("---------u==%lf\tu_std==%f---------\n", u, u_std[flag]);
	if (u < u_std[flag])
	{
		result = 1;
	}else
	{
		result = 0;
	}
	return result;
}
/**
 * 适线检验
 * 需要Qc
**/
int chk_SX(double *Q, double *Qc, int N) {
	double *r, u = 0, u_std[2] = {1.64, 1.28};
	float k_0 = 0, k_1 = 0, k_min=0, XZXSP[3] = {0.05, 0.10, 0.25}, XZXSP_tmp = 0;
	int result = 0, flag=0;
	int  accdny = 1;
	r=(double *)malloc(sizeof(double)*(N));
	printf("========进入适线检验=========\n");
	for (int i = 0; i < N; i++)
	{
		r[i] = Q[i] - Qc[i];
//		printf("<debug--sx>Q[%d]==%lf\tQc[%d]==%lf\tr[%d]==%lf\n",i,Q[i],i,Qc[i],i,r[i]);
		if (i>=1&&r[i]*r[i-1]<0)
		{
			k_1 += 1;
		}else if(i>=1&&r[i]*r[i-1]>0)
		{
			k_0 += 1;
		}
	}
	printf("<debug--SX> k_0==%f\tk_1==%f\n",k_0,k_1);
	if (k_0<k_1)
	{
		k_min = k_0;
	}else
	{
		k_min = k_1;
	}
	printf("---------默认显著性水平α=0.05，是否接受？(接受请输入1，否则为0):");
	scanf("%d", &accdny);
	if (accdny == 0)
	{
		printf("---------请输入显著性水平α(0.05or0.10or0.25):");
		scanf("%f", &XZXSP_tmp);
		for (int i = 0; i < 3; i++)
		{
			if (XZXSP_tmp==XZXSP[i])
				{
					flag = i;
				}
		}
	}
	//printf("<debug--sx>  N==%lf\tk_min==%lf\t\n", (double)N, (double)k_min);
	u = (0.5 * ((double)(N - 1)) - (double)k_min - 0.5) / (0.5 * sqrt((double)(N - 1)));
	printf("---------u==%lf\tu_std==%f---------\n", u, u_std[flag]);
	if (u < u_std[flag])
	{
		result = 1;
	}else
	{
		result = 0;
	}
	free(r);
	r=NULL;
	return result;
}

/**
 * ===此处进行偏离数值检验===
 * t=P_ave/S_p_ave   --->   统计量
 * P_ave=sum_P[i]/N   --->   测点与关系曲线的偏离值的平均值
 * P[i]=(Q[i]-Qc[i])/Qc[i]   --->   测点与曲线的偏离值
 * S_p_ave=S/sqrt(N)=sqrt(sum((Pi-P_ave)^2)/(n*(n-1)))   --->   P_ave的标准差
 * 其中，sum((P[i]-P_ave)^2)有一个简便算法
 * 即   sum((P[i]-P_ave)^2)=sum(P[i]^2)-N*P_ave^2   注：参考 《一种求sum((P[i]-P_ave)^2)的简便算法》，王君，苑洪洁，于汪洋，《黑龙江水专学报》第30卷第4期
 *
**/
int chk_PLSZ(double *Q,double *Qc, int N)
{
	/**此处P[N]为保存P[i]的数组，sum_p为求平均值时的和，sum_P_2为P[i]的平方**/
	double *P, sum_P = 0, P_ave = 0, S_p_ave = 0, sum_P_2 = 0, t = 0;
	/**此处t_std[][]为显著性水平以及自由度的临界值t[1-α/2]表**/
	double t_std[4][11] = { {4.30, 3.18, 2.78, 2.57, 2.45, 2.31, 2.23, 2.13, 2.09, 2.04, 2.00},
				{2.92, 2.35, 2.13, 2.02, 1.94, 1.86, 1.81, 1.75, 1.73, 1.70, 1.67},
				{1.89, 1.89, 1.53, 1.48, 1.44, 1.40, 1.37, 1.34, 1.33, 1.31, 1.30},
				{1.39, 1.39, 1.19, 1.10, 1.13, 1.11, 1.09, 1.07, 1.06, 1.06, 1.05}};
	/**K[11]为上述表的K值，XZXSP[4]为上述表的α值，此处采用‘《基于最小二乘法的绳套型水位流量关系最优定线研究》(董晓华)’一文中提供的相关表**/
	float k[11] = {2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 60}, XZXSP[4] = {0.05, 0.10, 0.20, 0.30}, XZXSP_tmp = 0;
	int accdny = 1;//用于判断用户接受or拒绝默认值
	/*flag_alpha为定位上述表行的标志，lag_k为定位上述表列的标志，result为最终判断结果，成功返回1，否则返回0，k_calc为用于计算的K，其值为N-1*/
	int flag_alpha=0, flag_k=0, result = 0, k_calc=0;
	printf("=========进入偏离数值检验！=========\n");
	printf("=========默认显著性水平α=0.05，是否接受？(接受请输入1，否则为0):");
	scanf("%d", &accdny);
	//printf("<debug--plsz>  accdny==%d\n", accdny);
	if (accdny == 0)
	{
			printf("---------请输入显著性水平α(0.05 or 0.10 or 0.20 or 0.30):\n");
			scanf("%f", &XZXSP_tmp);
			for (int i = 0; i < 4; i++)
			{
				if (XZXSP_tmp==XZXSP[i])
				{
					flag_alpha = i;
				}
				else
				{
					printf("---------输入显著性水平错误! 请重新输入！\n");
					return 0;
				}
			}

	}
	/*计算P[i]*/
	P=(double *)malloc(sizeof(double)*N);
	for (int i = 0; i < N; i++)
	{
		P[i] = (Q[i] - Qc[i]) / Qc[i];
		//printf("<debug--plsz>  P[%d]==%lf\t", i, P[i]);
		sum_P += P[i];
		sum_P_2 += P[i] * P[i];
		//printf("<debug--plsz>  sum_P_2==%lf\n", sum_P_2);
	}
	P_ave = sum_P / N;
	//printf("<debug--plsz>  P_ave=%lf\n", P_ave);
	S_p_ave = sqrt((sum_P_2 - N * P_ave * P_ave)/(N*(N-1)));
	t = P_ave / S_p_ave;
	//printf("<debug--plsz>  S_P_ave=%lf\n", S_p_ave);
	/*取绝对值*/
	if (t<0)
	{
		t = -t;
	}
	//printf("<debug--plsz>  t=%lf\n", t);
	k_calc = N - 1;
	/*计算并比较t与t_std*/
	for (int i = 0; i < 11; i++)
	{
		//由于上表前五列为连续K，因此不用参与else里的计算，此处用k_calc判断
		if (k_calc <= 6 && k_calc==k[i])
		{
			flag_k = i;
			//判断结果
			if (t<t_std[flag_alpha][flag_k])
			{
				result = 1;
			}
		//若所得k_calc在表中没有直接对应的数值，则采用按权重分配相邻两数据之间的差值的办法
		}else if (k_calc > 6 && k_calc >= k[i] && k_calc <= k[i+1])
		{
			int n = 0;
			n = k[i + 1] - k[i];
			float t_std_calc = 0;
			t_std_calc = t_std[flag_alpha][i] + ((k_calc-k[i])/(n)) * (t_std[flag_alpha][i + 1] - t_std[flag_alpha][i]);
			printf("<debug--plsz>  t==%lf\tt_std-calc=%lf\n", t,  t_std_calc);
			//判断结果
			if (t<t_std_calc)
			{
				result = 1;
			}
		}
	}
	return result;
}

double calc_Q_unknow(double Z_0, int m, double *alpha, double *beta, double *A)
{
	double Q_unknow=0, Y_unknow=0;
	double **P_unknow=NULL;
	double X[1]={0}, Z=0;
	printf("进入Q_unknow\n");
	printf("请输入水位Z=");
	scanf("%lf",&Z);
	X[0]=log(Z - Z_0);
//	printf("<debug__Q_unknow> X[0]=%lf\n",X[0]);
	P_unknow=(double **)malloc(sizeof(double *) * (m+1));
	for(int i=0;i<m+1;i++)
		P_unknow[i]=(double *)malloc(sizeof(double )*1);
	/** 初始化p_unknow **/
	for(int i=0;i<m+1;i++)
		P_unknow[i][0]=1;
//	printf("<debug__Q_unknow> P_unknow[1][0]=%lf\n",P_unknow[1][0]);
//		printf("<debug__Q_unknow> P_unknow[%d][%d]=%lf\n",m,0,P_unknow[m][0]);
	/** 计算P——unknow**/	
	for(int i=1;i<m+1;i++)
		P_unknow[i][0] = calc_P(X,alpha,P_unknow,beta,m,0);
//		printf("<debug__Q_unknow> P_unknow==%lf\n",P_unknow[m][0]);
	for(int i=0;i<m+1;i++)
		Y_unknow=Y_unknow+A[i]*P_unknow[i][0];
//	printf("<debug__Q_unknow> Q_unknow==%lf\n",Y_unknow);
	/** 计算拟合流量Q **/
	Q_unknow=exp(Y_unknow);
//	printf("<debug__Q_unknow> Q_unknow==%lf\n",Q_unknow);
	return Q_unknow;
}
