//------------------------------------------------------------
// 標準的な差分進化のアルゴリズムを実装したプログラム
//------------------------------------------------------------
//------------------------------------------------------------
//「error C2065: 'M_PI' : 定義されていない識別子です」への対応
//　VisualStudioの場合は以下のコメントを外してコンパイル

//bestを降順に並び替えて100p%の確率(p=0.05より上位5位まで)でランダムにbestを選ぶ．
//やらなければならないことはbestの並び替え，bestが上位100p%になるようにする．

// 以下はメルセンヌツイスターで使用する
#define NP	100	//最大個体数
#define D		30				//最大次元数（問題の次元数）
#define FUNC_NO		1					//最適化問題の種類
#define RANGE			5.12				//最適化問題の定義域
#define MAXGRNRATION	5000		//最大繰り返し回数
#define DE_ALGORITHM_NO	1	//DEのアルゴリズム
#define EXTIME			5//試行回数
#define Terminate		1.0e-11			//終了条件
#define Fl 0.1
#define Fu 0.9
#define Tf 0.1
#define Tcr 0.1
#define P_BEST 0.10		//標準は0.05(5%)

double sum = 0;		//合計
double ave[4][6];		//平均
double stv[4][6];	//標準偏差
double aaa[EXTIME][NP];
double nFitness_bot[NP];
//------------------------------------------------------------

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//#include <iostream>
#include <vector>

// 以下はメルセンヌツイスターで使用する
#include "random.h"
#include "func.h"
//------------------------------------------------------------
// メルセンヌツイスターで使用する関数の概説
// genrand_int32() //符号なし32ビット長整数
// genrand_int31() //符号なし31ビット長整数
// genrand_real1() //一様実乱数[0,1] (32ビット精度)
// genrand_real2() //一様実乱数[0,1) (32ビット精度)
// genrand_real3() //一様実乱数(0,1) (32ビット精度)
// genrand_res53() //一様実乱数[0,1) (53ビット精度)

#define _CRT_SECURE_NO_DEPRECATE
#pragma warning(disable:4996)
//------------------------------------------------------------
// 設定パラメータ 実験時にこの部分を書き換える
//------------------------------------------------------------


//------------------------------------------------------------
//ファイル出力1　最良値の推移
//------------------------------------------------------------
void Output_To_File1(void)
{
	int i, j;					//繰返し用変数.
	FILE *fp;					//ファイルポインタ
	char filename[80];			//ファイル名
	time_t timer;				//時間計測用
	struct tm *t_st;			//時間計測用
	time(&timer);				//時間の取得
	t_st = localtime(&timer);	//時間の変換
	sprintf_s(filename, "DE_gBestHistory%04d%02d%02d%02d%02d_DE_NO%d_FUNC_NO%d_NP%d_D%d_C%lf.csv",
		t_st->tm_year + 1900, t_st->tm_mon + 1,
		t_st->tm_mday, t_st->tm_hour, t_st->tm_min, DeAlgorithmNo, Func_No, Np, d, C);
	fp = fopen(filename, "a");
	for (i = 0; i < MaxGrnration; i++) {
		for (j = 0; j < EXTIME; j++) {
			fprintf(fp, "%20.70lf,", gBestHistory[j][i]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);

}
//------------------------------------------------------------
//集団の多様性の評価
//------------------------------------------------------------
void Calc_Diversity()
{
	int i, j;					//繰返し用変数.

	double div = 0;	//分散

	//合計
	for (i = 0; i < EXTIME; i++) {
		for (j = 0; j < Np; j++) {
			sum += aaa[i][j];
		}
	}
	printf("sum%lf,", sum);
	if (sum != 0) {
		//平均
		ave[Func_No - 1][DeAlgorithmNo - 1] = sum / (Np*EXTIME);

		//分散
		for (i = 0; i < EXTIME; i++) {
			for (j = 0, div = 0.0; j < Np; j++) {
				div += pow((ave[Func_No - 1][DeAlgorithmNo - 1] - aaa[i][j]), 2.0);
			}
		}
		stv[Func_No - 1][DeAlgorithmNo - 1] += sqrt(div / (Np*EXTIME));
	}
	else {
		stv[Func_No - 1][DeAlgorithmNo - 1] = 0;
	}
	printf("stv%f,ave%f\n", stv[Func_No - 1][DeAlgorithmNo - 1], ave[Func_No - 1][DeAlgorithmNo - 1]);
	sum = 0;


	for (i = 0; i < EXTIME; i++) {
		for (j = 0; j < Np; j++) {
			aaa[i][j] = 0;
		}
	}


}

bool hoge = 0;
void Calc_Diversity2()
{
	int i, j;					//繰返し用変数.

	double div = 0;	//分散

	//合計
	for (i = 0; i < EXTIME; i++) {
		sum += gBestHistory[i][MaxGrnration - 1];
	}
	printf("sum%lf,", sum);
	if (sum != 0) {
		//平均
		ave[Func_No - 1][DeAlgorithmNo - 1] = sum / (EXTIME);

		//分散
		for (i = 0; i < EXTIME; i++) {
			div += pow((ave[Func_No - 1][DeAlgorithmNo - 1] - gBestHistory[i][MaxGrnration - 1]), 2.0);
		}
		stv[Func_No - 1][DeAlgorithmNo - 1] += sqrt(div / (EXTIME));
	}
	else {
		stv[Func_No - 1][DeAlgorithmNo - 1] = 0;
	}
	printf("");
	for (i = 0; i < MaxGrnration;i++) {
//		QuickSort(gBestHistoryMed,0,EXTIME-1,i);
		bubbleSort_cnter(gBestHistoryMed,i);

	}
	printf("");
	for (i = 0; i < MaxGrnration; i++) {
		gBestHistoryMedian[Func_No - 1][DeAlgorithmNo - 1][i] = gBestHistoryMed[i][(int)((EXTIME/2)+0.5)-1];
		printf("%.11lf\n", gBestHistoryMed[i][((int)EXTIME / 2) - 1]);

		for (j = 0; j < EXTIME; j++) {
			gBestHistoryMed[i][j] = 0;
		}
	}

	printf("stv%f,ave%f\n", stv[Func_No - 1][DeAlgorithmNo - 1], ave[Func_No - 1][DeAlgorithmNo - 1]);
	sum = 0;
	hoge = 0;

}
//------------------------------------------------------------
//ファイル出力2　分散値の推移
//------------------------------------------------------------
void Output_To_File2(void)
{
	int i, j;					//繰返し用変数.
	FILE *fp;					//ファイルポインタ
	char filename[50];			//ファイル名
	time_t timer;				//時間計測用
	struct tm *t_st;			//時間計測用
	time(&timer);				//時間の取得
	t_st = localtime(&timer);	//時間の変換
	sprintf_s(filename, "DE_pDiversity%04d%02d%02d%02d%02d_DE_NO%d_FUNC_NO%d.txt",
		t_st->tm_year + 1900, t_st->tm_mon + 1,
		t_st->tm_mday, t_st->tm_hour, t_st->tm_min, DeAlgorithmNo, Func_No);
	fp = fopen(filename, "a");
	for (i = 0; i < MaxGrnration; i++) {
		for (j = 0; j < EXTIME; j++) {
			fprintf(fp, "%20.35lf\t", pDiversity[j][i]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
	fclose(fp);
}
//------------------------------------------------------------
//ファイル出力3　最適解発見世代，収束率
//------------------------------------------------------------
void Output_To_File3(void)
{
	register int i;
	FILE *fp;
	char filename[50];
	time_t timer;
	struct tm *t_st;
	time(&timer);
	t_st = localtime(&timer);
	sprintf_s(filename, "DE_sRate_gTable%04d%02d%02d%02d%02d_DE_NO%d_FUNC_NO%d_c=%lf.txt",
		t_st->tm_year + 1900, t_st->tm_mon + 1,
		t_st->tm_mday, t_st->tm_hour, t_st->tm_min, DeAlgorithmNo, Func_No, C);
	fp = fopen(filename, "a");
	for (i = 0; i < EXTIME; i++) {
		fprintf(fp, "%6d %d\n", gTable[i], sRate[i]);
	}
	fclose(fp);
}

void Output_To_File4(void)
{
	int i, j;					//繰返し用変数.
	double gBestHistorySum = 0;
	FILE *fp;					//ファイルポインタ
	char filename[100];			//ファイル名
	time_t timer;				//時間計測用
	struct tm *t_st;			//時間計測用
	time(&timer);				//時間の取得
	t_st = localtime(&timer);	//時間の変換
	sprintf_s(filename, "DE_Uf_rand_DE_NO%d_FUNC_NO%d_NP%d_D%d_ave.txt",
		DeAlgorithmNo, Func_No, NP, D);
	fp = fopen(filename, "a");

	for (i = 0; i < MaxGrnration; i++) {
		for (j = 0; j < EXTIME; j++) {
			fprintf(fp, "%.5lf\t", Uf_rand_History[j][i]);
			printf("%lf\n", Uf_rand_History[j][i]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}

void Output_To_File5(void)
{
	int i, j;					//繰返し用変数.
	double gBestHistorySum = 0;
	FILE *fp;					//ファイルポインタ
	char filename[100];			//ファイル名
	time_t timer;				//時間計測用
	struct tm *t_st;			//時間計測用
	time(&timer);				//時間の取得
	t_st = localtime(&timer);	//時間の変換
	sprintf_s(filename, "gBestHistoryAll_Func_No%d_NP%d_D%d_MaxG%d_Extime%d_ave.csv", Func_No, NP, D, MaxGrnration, EXTIME);
	fp = fopen(filename, "a");

	for (i = 0; i < MaxGrnration; i++) {
		for (j = 0; j < DeAlgorithmNo; j++) {
			fprintf(fp, "%.11lf,", gBestHistoryAll[Func_No - 1][j][i] / EXTIME);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}

void Output_To_File6(void)
{
	int i, j;					//繰返し用変数.
	double gBestHistorySum = 0;
	FILE *fp;					//ファイルポインタ
	char filename[100];			//ファイル名
	time_t timer;				//時間計測用
	struct tm *t_st;			//時間計測用
	time(&timer);				//時間の取得
	t_st = localtime(&timer);	//時間の変換
	sprintf_s(filename, "gBestHistoryAll_Func_No%d_NP%d_D%d_MaxG%d_Extime%d_median.csv", Func_No, NP, D, MaxGrnration, EXTIME);
	fp = fopen(filename, "a");

	for (i = 0; i < MaxGrnration; i++) {
		for (j = 0; j < DeAlgorithmNo ; j++) {
			fprintf(fp, "%.11lf,", gBestHistoryMedian[Func_No - 1][j][i]);

//			printf("%.11lf\n", gBestHistoryMedian[Func_No - 1][j][i]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}

void Output_To_File7(void)
{

	int i, j;					//繰返し用変数.
	double gBestHistorySum = 0;
	FILE *fp;					//ファイルポインタ
	char filename[120];			//ファイル名
	time_t timer;				//時間計測用
	struct tm *t_st;			//時間計測用
	time(&timer);				//時間の取得
	t_st = localtime(&timer);	//時間の変換
	sprintf_s(filename, "ave_stv_Func_No%d_NP%d_D%d_MaxG%d_Extime%d.csv", Func_No, NP, D, MaxGrnration, EXTIME);
	fp = fopen(filename, "a");
	for (j = 0; j < DeAlgorithmNo - 1; j++) {
		fprintf(fp, "func_no%d_DE_No=%d,", Func_No, j);
		fprintf(fp, "ave=%.11lf,", ave[Func_No - 1][j]);
		fprintf(fp, "stv=%.11lf,", stv[Func_No - 1][j]);
		fprintf(fp, "\n");
	}

	fclose(fp);
}




//------------------------------------------------------------
//メイン関数
//------------------------------------------------------------
int main(void)
{
	int pop;
	int J_rand = 0;
	int j = 0, i = 0, k = 0,l=0;
	int count_data = 0;
	p_best = (int)(NP*P_BEST); //上位p%の設定
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 6; j++) {
			for (k = 0; k < MAXGRNRATION; k++) {
				gBestHistoryAll[i][j][k] = 0;
				gBestHistoryMedian[i][j][k] = 0;
				for (l = 0; l < EXTIME;l++) 	gBestHistoryMed[k][l] = 0;
			}

			stv[i][j] = 0;
			ave[i][j] = 0;
		}
	}
	Uf_best_History[EXTIME][MAXGRNRATION], Uf_rand_History[EXTIME][MAXGRNRATION];
	for (Func_No = 1; Func_No <= 4; Func_No++) {
		for (DeAlgorithmNo = 1; DeAlgorithmNo <= 6; DeAlgorithmNo++) {
			//if(DeAlgorithmNo ==5)	Output_To_File4();
			printf("DeAlgorithmNo=%d\nFunc_No=%d\n", DeAlgorithmNo, Func_No);
			pop = 0, J_rand = 0, j = 0, i = 0;
			vDEParameter();//DEのパラメータの設定
			best_Initialize();
			init_genrand((unsigned)time(NULL));	//MTの初期化
			for (iteration = 0; iteration < EXTIME; iteration++) {//試行回数
				Initialize();		//初期化
				while (episode < MaxGrnration) {
					Select_Elite_Vector(iteration, episode);
					//Calc_Diversity(iteration, episode);
					for (pop = 0; pop < Np; pop++) {
						Select_pVector(pop);
						DE_Operation(pop, episode);
						Evaluate_New_Vector(pop);
						//printf("%d Ucr%lf Uf_best%lf Uf_rand%lf \n", episode, Ucr, Uf_best,Uf_rand);
					}

					Parameter_Format(DeAlgorithmNo);
					Compare_Vector();

					if (DeAlgorithmNo == 5 || DeAlgorithmNo == 6) {
						bubbleSort(nFitness);
					}
					episode++;
					if (Terminate > gBestFitness) break;
				}

				gTable[iteration] = episode;
				if (gBestFitness < Terminate)sRate[iteration] = 1;
				else sRate[iteration] = 0;
				gBestTable[iteration] = gBestFitness;

			}
			Calc_Diversity2();
		Output_To_File1();
			//Output_To_File2();
			//Output_To_File3();
/*			if (DeAlgorithmNo == 5 || DeAlgorithmNo==6|| DeAlgorithmNo==7) {
				Output_To_File4();
				Output_To_File5();
			}
*/
		}
		Output_To_File5();
		Output_To_File6();
		Output_To_File7();
	}

	return 0;
}


/*
for (i = 0; i < Np; i++) {
for (j = 0; j < d; j++) {
printf("変化前nVect[%d][%d]=%0.10f\n", i, j, nVect[i][j]);
}
printf("\n");
}
printf("\n\n");
for (i = 0; i < Np; i++) {
printf("変化前nFitness[%d]=%0.10f\n", i, nFitness[i]);
}
*/

//QuickSort();

/*
printf("\n");
for (i = 0; i < Np; i++) {
printf("変化後nFitness[%d]=%0.10f\n", i, nFitness[i]);
}
printf("\n");
for (i = 0; i < Np; i++) {
for (j = 0; j < d; j++) {
printf("変化後nVect[%d][%d]=%0.10f\n", i, j, nVect[i][j]);
}
printf("\n");
}
printf("\n\n");
*/

//					printf("Uf=%5.2f \n Ucr=%5.2f\n", Uf, Ucr);




/*
for (i = 0; i < Np; i++) {
printf("変化後nFitness[%d]=%0.10f\n", i, nFitness[i]);
}
printf("\n");
for (i = 0; i < Np; i++) {
for (j = 0; j < d; j++) {
printf("変化後nVect[%d][%d]=%0.10f\n", i, j, nVect[i][j]);
}
printf("\n");
}
*/
//printf("CR%lf\tF%lf\n", average_CR, average_F);