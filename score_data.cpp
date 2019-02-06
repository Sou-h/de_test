#define NP	100	//最大個体数
#define D	30	//最大次元数（問題の次元数）
#define FUNC_NO	1	//最適化問題の種類
#define RANGE	5.12	//最適化問題の定義域
#define MAXGRNRATION	5000	//最大繰り返し回数
#define DE_ALGORITHM_NO	1	//DEのアルゴリズム
#define EXTIME	100	//試行回数
#define Terminate	1.0e-11	//終了条件
#define Fl 0.1	//jDEのパラメータ
#define Fu 0.9	//jDEのパラメータ
#define Tf 0.1	//jDEのパラメータ
#define Tcr 0.1	//jDEのパラメータ
#define P_BEST 0.10	//標準は0.05(5%)
#define P_BEST_COUNT	6
#define _USE_MATH_DEFINES

double sum_data = 0;	//合計
double ave[4][6];	//平均
double stv[4][6];	//標準偏差


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>

#include "random.h"
#include "func_.h"
// 以下はメルセンヌツイスターで使用する
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
//ファイル出力1　最良値の推移
//------------------------------------------------------------
void Output_To_File1(void)
{
int i, j;	//繰返し用変数.
FILE *fp;	//ファイルポインタ
char filename[80];	//ファイル名
time_t timer;		//時間計測用
struct tm *t_st;	//時間計測用
time(&timer);		//時間の取得
t_st = localtime(&timer);//時間の変換
sprintf_s(filename, "DE_gBestHistory%04d%02d%02d%02d%02d_DE_NO%d_FUNC_NO%d_NP%d_D%d_C%.2lf.csv",
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


//分散の抽出
void Calc_Diversity2()
{
int i, j;					//繰返し用変数.
double div = 0;	//分散

//合計
for (i = 0; i < EXTIME; i++) {
	sum_data += gBestHistory[i][MaxGrnration - 1];
}
if (sum_data != 0) {
	//平均
	ave[Func_No - 1][DeAlgorithmNo - 1] = sum_data / (EXTIME);

	//分散
	for (i = 0; i < EXTIME; i++) {
		div += pow((ave[Func_No - 1][DeAlgorithmNo - 1] 
			- gBestHistory[i][MaxGrnration - 1]), 2.0);
	}
	stv[Func_No - 1][DeAlgorithmNo - 1] += sqrt(div / (EXTIME));
}
else {
	stv[Func_No - 1][DeAlgorithmNo - 1] = 0;
}
for (i = 0; i < MaxGrnration;i++) {
	bubbleSort_cnter(gBestHistoryMed,i);

}
for (i = 0; i < MaxGrnration; i++) {
	gBestHistoryMedian[Func_No - 1][DeAlgorithmNo - 1][i] = 
		gBestHistoryMed[i][(int)((EXTIME/2)+0.5)-1];

	for (j = 0; j < EXTIME; j++) {
		gBestHistoryMed[i][j] = 0;
	}
}

sum_data = 0;

}
//平均値の出力
void Output_To_File5(void)
{
int i, j;	//繰返し用変数.
FILE *fp;	//ファイルポインタ
char filename[100];	//ファイル名
time_t timer;		//時間計測用
struct tm *t_st;	//時間計測用
time(&timer);		//時間の取得
t_st = localtime(&timer);	//時間の変換
sprintf_s(filename, "gBestHistoryAll_Func_No%d_NP%d_D%d_MaxG%d_Extime%d_ave.csv"
	, Func_No, NP, D, MaxGrnration, EXTIME);
fp = fopen(filename, "a");

for (i = 0; i < MaxGrnration; i++) {
	for (j = 0; j < DeAlgorithmNo; j++) {
		fprintf(fp, "%.11lf,", gBestHistoryAll[Func_No - 1][j][i] / EXTIME);
	}
	fprintf(fp, "\n");
}

fclose(fp);
}
//中央値の出力
void Output_To_File6(void)
{
int i, j;	//繰返し用変数.
FILE *fp;	//ファイルポインタ
char filename[100];	//ファイル名
time_t timer;		//時間計測用
struct tm *t_st;	//時間計測用
time(&timer);		//時間の取得
t_st = localtime(&timer);	//時間の変換
sprintf_s(filename, "gBestHistoryAll_Func_No%d_NP%d_D%d_MaxG%d_Extime%d_median.csv"
	, Func_No, NP, D, MaxGrnration, EXTIME);
fp = fopen(filename, "a");

for (i = 0; i < MaxGrnration; i++) {
	for (j = 0; j < DeAlgorithmNo ; j++) {
		fprintf(fp, "%.11lf,", gBestHistoryMedian[Func_No - 1][j][i]);

	}
	fprintf(fp, "\n");
}

fclose(fp);
}
//平均値と標準偏差の出力
void Output_To_File7(void)
{
int  i;		//繰返し用変数.
FILE *fp;	//ファイルポインタ
char filename[120];	//ファイル名
time_t timer;		//時間計測用
struct tm *t_st;	//時間計測用
time(&timer);		//時間の取得
t_st = localtime(&timer);	//時間の変換
sprintf_s(filename, "ave_stv_Func_No%d_NP%d_D%d_MaxG%d_Extime%d.csv"
	, Func_No, NP, D, MaxGrnration, EXTIME);
fp = fopen(filename, "a");
for (i = 0; i < DeAlgorithmNo - 1; i++) {
	fprintf(fp, "func_no%d_DE_No=%d,", Func_No, i);
	fprintf(fp, "ave=%.11lf,", ave[Func_No - 1][i]);
	fprintf(fp, "stv=%.11lf,", stv[Func_No - 1][i]);
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
p_best = (int)(NP*P_BEST); //上位p%の設定JADE,jDEで使用
First_time_only_Initialize();
for (Func_No = 1; Func_No <= 1; Func_No++) {
	for (DeAlgorithmNo = 5; DeAlgorithmNo <= 5; DeAlgorithmNo++) {
		printf("DeAlgorithmNo=%d\nFunc_No=%d\n", DeAlgorithmNo, Func_No);
		vDEParameter();//DEのパラメータの設定
		best_Initialize();//gBestHistoryの初期化
		for (iteration = 0; iteration < EXTIME; iteration++) {//試行回数
			Initialize();	//初期化
			while (episode < MaxGrnration) {//世代数
				Select_Elite_Vector(iteration, episode);//最良値の格納
				for (pop = 0; pop < Np; pop++) {
					Select_pVector(pop);//p1,p2,p3の選択
					DE_Operation(pop, episode);//交叉,突然変異
					Evaluate_New_Vector(pop);//親子の比較
				}
				//JADE,jDEの際にパラメータの更新
				Parameter_Format(DeAlgorithmNo);
				Compare_Vector();
				//JADE,jDEの際にソートを実行
				bubbleSort(nFitness);
				episode++;
				if (Terminate > gBestFitness) break; //目標値との比較
			}
		}
		Calc_Diversity2();//分散の抽出
		Output_To_File1();//最良値の出力
	}
	Output_To_File5();//平均値の出力
	Output_To_File6();//中央値の出力
	Output_To_File7();//平均値と標準偏差の出力
}
return 0;
}
