//------------------------------------------------------------
// 標準的な差分進化のアルゴリズムを実装したプログラム
//------------------------------------------------------------
//------------------------------------------------------------
//「error C2065: 'M_PI' : 定義されていない識別子です」への対応
//　VisualStudioの場合は以下のコメントを外してコンパイル
//------------------------------------------------------------
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
// 以下はメルセンヌツイスターで使用する
#include "random.h"
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

#define NP	30				//最大個体数
#define D		30				//最大次元数（問題の次元数）
#define FUNC_NO		1					//最適化問題の種類
#define RANGE		5.12				//最適化問題の定義域
#define MRATE			0.8				//突然変異率
#define CRATE			0.6				//交叉率
#define MAXGRNRATION	5000		//最大繰り返し回数
#define DE_ALGORITHM_NO	1	//DEのアルゴリズム
#define EXTIME			50		//試行回数
#define Terminate		1.0e-55			//終了条件

//-----------------------------------------------------------
int Np = NP;				//最大個体数
int d = D;		//最大次元数（問題の次元数）
int Func_No = FUNC_NO;					//最適化問題の種類
double Range = RANGE;				//最適化問題の定義域
int MaxGrnration = MAXGRNRATION;		//最大繰り返し回数
int DeAlgorithmNo = DE_ALGORITHM_NO;	//DEのアルゴリズム
//------------------------------------------------------------
double cVect[NP][D];				//時刻Tの個体（ベクトル）
double cFitness[NP];				//時刻Tの個体の評価値
double pBestVector[NP][D];	//時刻Tにおける各個体の最良解の履歴
double pBestFitness[NP];		//時刻Tにおける各個体の最良解評価値の履歴
double gBestVector[D];			//時刻Tにおける集団全体での最良解の履歴
double gBestFitness;				//時刻Tにおける集団全体での最良解評価値の履歴
double pVect1[D];					//親ベクトル１
double pVect2[D];					//親ベクトル２
double pVect3[D];					//親ベクトル３
double nVect[NP][D];				//時刻T+1の個体（ベクトル）
double nFitness[NP];				//時刻T+1の個体の評価値
double gBestHistory[EXTIME][MAXGRNRATION];	//各試行におけるgBestの履歴
double pDiversity[EXTIME][MAXGRNRATION];		//各試行における集団の多様性
int gTable[EXTIME];					//最適解発見時の世代数
int sRate[EXTIME];						//最適解の発見率
double gBestTable[EXTIME];		//終了時点でのgBestの値
double C = 0.1;							//JADE のパラメータ
double Ucr = 0.5;
double Uf = 0.5;
double Sf[NP];
double Scr[NP];
double CR[NP];
double F[NP];
double average_F = 0, average_CR = 0, average_Sn = 0;
FILE *fp_4;					//ファイルポインタ
char filename_4[100];			//ファイル名


//------------------------------------------------------------
//------------------------------------------------------------
// [min, max]の一様乱数
//------------------------------------------------------------
double uniform(double min, double max) {
	return min + (max - min)*genrand_real1();
}
//------------------------------------------------------------
// テスト関数
// F1 Sphere関数
// [-5.12, +5.12] x=(0,..,0) F1=0 
//------------------------------------------------------------
double Sphere(double *x) {
	int i;
	double sum = 0.0;
	for (i = 0; i<d; i++) {
		sum += x[i] * x[i];
	}
	return sum;
}
//------------------------------------------------------------
// F2 Rosenbrock関数 Chain型
// [-2.048, +2.048] x=(1,..,1) F2=0 
//------------------------------------------------------------
double Rosenbrock(double *x) {
	int i;
	double sum = 0.0;
	/*
	for (i = 0; i<d - 1; i++){
	sum += 100 * (x[i + 1] - x[i] * x[i])*(x[i + 1] - x[i] * x[i]) + (x[i] - 1.0)*(x[i] - 1.0);
	}
	*/
	// UNDXの文献の式に統一
	// Rosenbrock関数 Star型
	for (i = 2; i<d; i++) {
		sum += (100 * (x[0] - x[i] * x[i])*(x[0] - x[i] * x[i]) + (x[1] - 1.0)*(x[1] - 1.0));
	}
	return sum;
}
//------------------------------------------------------------
// F3 Rastrigin関数
//------------------------------------------------------------
double Rastrigin(double *x) {
	int i;
	double sum1 = 0, sum2 = 0;
	for (i = 0; i<d; i++) {
		sum1 += x[i] * x[i] - 10 * cos(2 * M_PI*x[i]);
	}
	return (10 * d + sum1);
}
//------------------------------------------------------------
// F4 Griewank関数
//------------------------------------------------------------
double Griewank(double *x) {
	int i;
	double sum1 = 0.0, sum2 = 1.0;
	for (i = 0; i<d; i++) {
		sum1 += x[i] * x[i];
		sum2 *= cos(x[i] / sqrt((i + 1)*1.0));
	}
	return ((sum1 / 4000) - sum2 + 1);
}
//------------------------------------------------------------
// F5 Ackley関数
//------------------------------------------------------------
double Ackley(double *x) {
	int i;
	double sum1, sum2;
	for (i = 0, sum1 = 0.0; i<d; i++)sum1 += x[i] * x[i];
	for (i = 0, sum2 = 0.0; i<d; i++)sum2 += cos(2.0*M_PI*d);
	return (20.0 + M_E - 20.0 * exp(-0.2*sqrt(sum1 / d)) - exp(sum2 / d));
}
//------------------------------------------------------------
// F6 Schwefel関数
//------------------------------------------------------------
double Schwefel(double *x) {
	int i;
	double sum;
	for (i = 0, sum = 0.0; i<d; i++) {
		sum += x[i] * sin(sqrt(fabs(x[i])));
	}
	return (sum + 418.9828872724338 * d);
}
//------------------------------------------------------------
// F7 Ridge関数
//------------------------------------------------------------
double Ridge(double *x) {
	int i, j;
	double sum1, sum2;
	for (i = 0, sum1 = 0; i<d; i++) {
		for (j = 0, sum2 = 0; j <= i; j++) {
			sum2 += x[j];
		}
		sum1 += sum2 * sum2;
	}
	return sum1;
}
//------------------------------------------------------------
// F8 Bohachevsky関数
//------------------------------------------------------------
double Bohachevsky(double *x) {
	int i;
	double sum;
	for (i = 0, sum = 0; i<d - 1; i++) {
		sum += x[i] * x[i] + 2 * x[i + 1] * x[i + 1] - 0.3*cos(3 * M_PI*x[i]) - 0.4*cos(4 * M_PI*x[i + 1]) + 0.7;
	}
	return sum;
}
//------------------------------------------------------------
// F9 Schaffer関数
//------------------------------------------------------------
double Schaffer(double *x) {
	int i;
	double tp1 = 0, tp2 = 0, tp3 = 0, tp4 = 0, tp5 = 0;
	for (i = 0, tp5 = 0; i<d - 1; i++) {
		tp1 = x[i] * x[i] + x[i + 1] * x[i + 1];
		tp2 = pow(tp1, 0.25);
		tp3 = pow(tp1, 0.1);
		tp4 = sin(50 * tp3) * sin(50 * tp3);
		tp5 += tp2 * (tp4 + 1.0);
	}
	return tp5;
}
//------------------------------------------------------------
// F10 Ellipsoid関数
//------------------------------------------------------------
double Ellipsoid(double *x) {
	int i;
	double sum;
	for (i = 0, sum = 0; i<d; i++) {
		sum += (pow(1000, (double)i / (double)(d - 1)) * x[i]) * (pow(1000, (double)i / (double)(d - 1)) * x[i]);
	}
	return sum;
}
//------------------------------------------------------------
// F11 k-tablet関数
//------------------------------------------------------------
double K_Tablet(double *x) {
	int i;
	double sum1, sum2;
	for (i = 0, sum1 = 0; i<(int)d / 2; i++) {
		sum1 += x[i] * x[i];
	}
	for (i = (int)d / 2 + 1, sum2 = 0; i<d; i++) {
		sum2 += (100 * x[i]) * (100 * x[i]);
	}
	return (sum1 + sum2);
}
//------------------------------------------------------------
// F12 Shifted-Rastrigin関数
//------------------------------------------------------------
double Shifted_Rastrigin(double *x) {
	int i;
	double sum1 = 0.0, sum2 = 0.0;
	for (i = 0; i<d; i++) sum1 += (x[i] - 1) * (x[i] - 1) - 10 * cos(2 * M_PI*(x[i] - 1));
	sum2 = 10 * d + sum1;
	return sum2;
}
//------------------------------------------------------------
// F13 Cigar関数
//------------------------------------------------------------
double Cigar(double *x) {
	int i;
	double sum;
	for (i = 1, sum = 0; i<d; i++) {
		sum += (1000 * x[i]) * (1000 * x[i]);
	}
	return (x[0] * x[0] + sum);
}
//------------------------------------------------------------
// F13 2n-minima関数
//------------------------------------------------------------
double N2_Minima(double *x) {
	int i;
	double sum;
	for (i = 0, sum = 0; i<d; i++) {
		sum += (pow(x[i], 4) - 16 * pow(x[i], 2) + 5 * x[i]);
	}
	return (sum / 2);
}
//------------------------------------------------------------
// 目的関数値を計算
//------------------------------------------------------------
double Calc_Objective_Function(double *x)
{
	if (Func_No == 1)return  Sphere(x);
	else if (Func_No == 2)return  Rosenbrock(x);
	else if (Func_No == 3)return  Rastrigin(x);
	else if (Func_No == 4)return  Griewank(x);
	else if (Func_No == 5)return  Ackley(x);
	else if (Func_No == 6)return  Schwefel(x);
	else if (Func_No == 7)return  Ridge(x);
	else if (Func_No == 8)return  Bohachevsky(x);
	else if (Func_No == 9)return  Schaffer(x);
	else if (Func_No == 10)return Ellipsoid(x);
	else if (Func_No == 11)return K_Tablet(x);
	else if (Func_No == 12)return Shifted_Rastrigin(x);
	else if (Func_No == 13)return Cigar(x);
	else if (Func_No == 14)return N2_Minima(x);
	else exit(0);
}
//------------------------------------------------------------
//初期個体群の生成
//------------------------------------------------------------
void Init_Vector(void)
{
	int i, j;		//繰り返し用変数
	double r;		//定義域用変数
	void vRange();
	r = Range;
	for (i = 0; i<Np; i++) {
		for (j = 0; j<d; j++) {
			cVect[i][j] = Range * (genrand_real1() * 2 - 1);
		}
	}
	//初期ベクトルをpBestVectorに保存する
	for (i = 0; i<Np; i++) {
		for (j = 0; j<d; j++) {
			pBestVector[i][j] = cVect[i][j];
		}
	}
}

void vRange() {
	if (Func_No == 1)		Range = 100;
	else if (Func_No == 2) Range=30;
	else if (Func_No == 3) Range=5.12;
	else if (Func_No == 4) Range=600;
	else if (Func_No == 5) Range=32.768;
	else if (Func_No == 6) Range=500;
	else if (Func_No == 9) Range=100;
	else if (Func_No == 10) Range =5.12;
	else {
		printf("Rangeにあうパラメータが設定されていません\nFunc_No=%d",Func_No);
		exit(0);
	}

}
//------------------------------------------------------------
// 初期集団の評価（目的関数値を適応度として利用）
//------------------------------------------------------------
void Evaluate_Init_Vector(void) {
	int i;
	for (i = 0; i<Np; i++) {
		cFitness[i] = Calc_Objective_Function(cVect[i]);
		//初期値を初期pBestとして保存
		pBestFitness[i] = cFitness[i];
	}
}
//------------------------------------------------------------
//新規ベクトルの生成のための親ベクトルの選択
//------------------------------------------------------------
void Select_pVector(int pop1)
{
	register int i;
	int pop2, pop3, pop4;
	do {
		pop2 = (int)Np*genrand_real1();
		pop3 = (int)Np*genrand_real1();
		pop4 = (int)Np*genrand_real1();
	} while (pop1 == pop2 || pop1 == pop3 || pop1 == pop4 || pop2 == pop3 || pop2 == pop4 || pop3 == pop4);
	for (i = 0; i<d; i++)pVect1[i] = cVect[pop2][i];
	for (i = 0; i<d; i++)pVect2[i] = cVect[pop3][i];
	for (i = 0; i<d; i++)pVect3[i] = cVect[pop4][i];
}

//パラメータの更新	JADE 作成中
void New_parameter() {
	int i, j;
	double Sn = 0.0, Sf = 0.0, Sf2 = 0.0, Scr = 0.0;
	double hoge_fp = 0;
	for (i = 0; i < Np; i++) {
		if (nFitness[i] > cFitness[i]) {
			//		if (nFitness[i] < gBestFitness) {
			Sn += 1;
			Sf += F[i];
			Sf2 += F[i] * F[i];
			Scr += CR[i];
//			if (cFitness[i] < gBestFitness) hoge_fp += 1;
//			printf("nfit%lf,gBest%f\n", cFitness[i], gBestFitness);
		}
	}

//	Pameter_Filter(Sf, Scr);
//	printf("hoge_fp=%lf\n", hoge_fp);
//	if (Sn > 10) {
	if (Sn != 0) {
//		Ucr = Sn / Np;
//		Uf = Sn / Np;
		Ucr = (1 - C)*Ucr + C * ((Scr / Sn));
		Uf = (1 - C)*Uf + C * (Sf2 / Sf);
//		Ucr = (1 - C)*Ucr + C * (Sf2 / Sf);
//		Uf = (1 - C)*Uf + C * ((Scr / Sn));
//		printf("Sf2=%lf Sf=%lf\n", Sf2, Sf);
//		printf("cScr=%lf cSf2=%lf\n",C* Sf2/Sf,C*Scr/Sn);
//		printf("Ucr%lf Uf%lf\n", Ucr, Uf);
//		printf("Sn=%.0lf\n", Sn);

	}else {
//		printf("何もなかった\t%d\n",Func_No);
//		printf("Sn=%.0lf\n", Sn);
//		getchar();
	}
//	Ucr = 0.3;
//	Uf = 0.3;
}
//パラメータの更新
void Parameter_Format() {
	int i = 0;
	double sum_F = 0, sum_CR = 0;

	//	printf("F=%lf\n",F[i]);
	//	printf("CR=%lf\n", CR[i]);
	for (i = 0; i < Np; i++) {
		do {
			F[i] = rand_cauchy(Uf, 0.1);
			if (F[i] > 1) F[i] = 1;
		} while (F[i] < 0);
		sum_F += F[i];
		do {
			CR[i] = rand_normal(Ucr, 0.1);
			if (CR[i] > 1) CR[i] = 1;
		} while (CR[i] < 0);
		sum_CR += CR[i];

	}

/*
			F[i] = rand_cauchy(Uf, 0.1);
			if (F[i] > 1) F[i] = 1;
			if (F[i] < 0) F[i] = 1.0e-10;

			CR[i] =rand_normal(Ucr, 0.1);
			if (CR[i] > 1) CR[i] = 1;
			if (CR[i] < 0) CR[i] = 1.0e-10;

	}
*/
	average_F = sum_F / Np;
	average_CR = sum_CR / Np;

	//	if (DeAlgorithmNo == 5) fprintf(fp_4, "%lf\t%lf\n", average_F, average_CR);

}
//パラメータの初期化
void Parameter_Initialization() {
	int i;
	for (i = 0; i < Np; i++) {
		F[i] = 0.5;
		CR[i] = 0.5;
	}
}

//------------------------------------------------------------
//DEの操作
//------------------------------------------------------------
void DE_Operation(int i_Np, int g_GSIZE)
{
	register int i;
	int N = 0, L = 0;
	//DE/rand/1/exp
	if (DeAlgorithmNo == 1) {
		for (i = 0; i<d; i++) nVect[i_Np][i] = cVect[i_Np][i];
		N = (int)(genrand_real1()*d);
		L = 0;
		do {
			nVect[i_Np][N] = pVect1[N] + MRATE * (pVect2[N] - pVect3[N]);
			if (nVect[i_Np][N] < -Range) nVect[i_Np][N] = pVect1[N] + genrand_real1() * (-Range - pVect1[N]);
			if (nVect[i_Np][N] >  Range) nVect[i_Np][N] = pVect1[N] + genrand_real1() * (Range - pVect1[N]);
			N = (N + 1) % d;
			L++;
		} while (genrand_real1() < CRATE && L < d);
	}
	//DE/best/1/exp
	else if (DeAlgorithmNo == 2) {
		for (i = 0; i<d; i++) nVect[i_Np][i] = cVect[i_Np][i];
		N = (int)(genrand_real1()*d);
		L = 0;
		do {
			nVect[i_Np][N] = gBestVector[N] + MRATE * (pVect2[N] - pVect3[N]);
			if (nVect[i_Np][N] < -Range) nVect[i_Np][N] = pVect1[N] + genrand_real1() * (-Range - pVect1[N]);
			if (nVect[i_Np][N] >  Range) nVect[i_Np][N] = pVect1[N] + genrand_real1() * (Range - pVect1[N]);
			N = (N + 1) % d;
			L++;
		} while (genrand_real1() < CR[i] && L < d);
	}
	//DE/rand/1/bin
	else if (DeAlgorithmNo == 3) {
		for (i = 0; i<d; i++) nVect[i_Np][i] = cVect[i_Np][i];
		N = (int)(genrand_real1()*d);
		for (L = 0; L<d; L++) {
			if (L == 0 || genrand_real1() < CRATE) {
				nVect[i_Np][N] = pVect1[N] + MRATE * (pVect2[N] - pVect3[N]);
				if (nVect[i_Np][N] < -Range) nVect[i_Np][N] = pVect1[N] + genrand_real1() * (-Range - pVect1[N]);
				if (nVect[i_Np][N] >  Range) nVect[i_Np][N] = pVect1[N] + genrand_real1() * (Range - pVect1[N]);
			}
			else {
				nVect[i_Np][N] = cVect[i_Np][N];
			}
			N = (N + 1) % d;
		}
	}
	//DE/best/1/bin
	else if (DeAlgorithmNo == 4) {
		for (i = 0; i<d; i++) nVect[i_Np][i] = cVect[i_Np][i];
		N = (int)(genrand_real1()*d);
		for (L = 0; L<d; L++) {
			if (L == 0 || genrand_real1() < CRATE) {
				nVect[i_Np][N] = gBestVector[N] + MRATE * (pVect2[N] - pVect3[N]);
				if (nVect[i_Np][N] < -Range) nVect[i_Np][N] = pVect1[N] + genrand_real1() * (-Range - pVect1[N]);
				if (nVect[i_Np][N] >  Range) nVect[i_Np][N] = pVect1[N] + genrand_real1() * (Range - pVect1[N]);
			}
			else {
				nVect[i_Np][N] = cVect[i_Np][N];
			}
			N = (N + 1) % d;
		}

	}
	//JADE/rand/1/bin
	else if (DeAlgorithmNo == 5) {
		for (i = 0; i < d; i++) nVect[i_Np][i] = cVect[i_Np][i];
		for (L = 0; L < d; L++) {
			if (L == 0 || genrand_real1() > CR[i]) {
				nVect[i_Np][N] = cVect[i_Np][N] + F[i_Np] * (gBestVector[N] - cVect[i_Np][N]) + F[i_Np] * (pVect1[N] - pVect2[N]);
//				nVect[i_Np][N] = pVect1[N] + MRATE * (pVect2[N] - pVect3[N]);
				if (nVect[i_Np][N] < -Range) nVect[i_Np][N] = pVect1[N] + genrand_real1() * (-Range - pVect1[N]);
				if (nVect[i_Np][N] > Range) nVect[i_Np][N] = pVect1[N] + genrand_real1() * (Range - pVect1[N]);
			}else {
				nVect[i_Np][N] = cVect[i_Np][N];
//				nVect[i_Np][N] = nVect[i_Np][N] = pVect1[N] + F[i_Np] * (pVect2[N] - pVect3[N]);
			}
			N = (N + 1) % d;
		}
		New_parameter();
	}
	//JADE	DE/rand/exp
	else if (DeAlgorithmNo == 6) {
		L = 0;
		for (i = 0; i < d; i++) nVect[i_Np][i] = cVect[i_Np][i];
		N = (int)(genrand_real1()*d);
			do{
				//nVect[i_Np][N] = pVect1[i_Np] + Mrate[i_Np] * (gBestVector[N] - pVect1[N]) + Mrate[i_Np] * (pVect1[i_Np] - pVect2[i_Np]);	//変更中
				nVect[i_Np][N] = cVect[i_Np][N] + F[i_Np] * (gBestVector[N] - cVect[i_Np][N]) + F[i_Np] * (pVect1[i_Np] - pVect2[i_Np]);
				if (nVect[i_Np][N] < -Range) nVect[i_Np][N] = pVect1[N] + genrand_real1() * (-Range - pVect1[N]);
				if (nVect[i_Np][N] > Range) nVect[i_Np][N] = pVect1[N] + genrand_real1() * (Range - pVect1[N]);
				N = (N + 1) % d;
				L++;
			}while (genrand_real1() < CR[i] && L < d);
			New_parameter();
	}
		//JADE	DE/curreny-to-pbest/1
	else if (DeAlgorithmNo == 7) {
		for (i = 0; i<d; i++) nVect[i_Np][i] = cVect[i_Np][i];
		N = (int)(genrand_real1()*d);
		L = 0;
		i = 0;
		do {
			nVect[i_Np][N] = pVect1[N] + F[i] * (pVect2[N] - pVect3[N]);
			if (nVect[i_Np][N] < -Range) nVect[i_Np][N] = pVect1[N] + genrand_real1() * (-Range - pVect1[N]);
			if (nVect[i_Np][N] >  Range) nVect[i_Np][N] = pVect1[N] + genrand_real1() * (Range - pVect1[N]);
			N = (N + 1) % d;
			L++;
			i++;
		} while (genrand_real1() < CR[i] && L < d);
		New_parameter();
	}
	else exit(0);
}

//------------------------------------------------------------
//新しいベクトルの評価
//------------------------------------------------------------
void Evaluate_New_Vector(int pop1)
{
	nFitness[pop1] = Calc_Objective_Function(nVect[pop1]);
}
//------------------------------------------------------------
//ベクトルの比較
//------------------------------------------------------------
void Compare_Vector(void)
{
	int i, j;				//繰返し用変数.
	for (i = 0; i<Np; i++) {
		//新しいベクトルが良ければ置き換え操作を行う
		if (nFitness[i] < cFitness[i]) {
			cFitness[i] = nFitness[i];
			for (j = 0; j<d; j++)cVect[i][j] = nVect[i][j];
		}
		else continue;
	}
}
//------------------------------------------------------------
//エリート選択
//------------------------------------------------------------
void Select_Elite_Vector(int itime, int gtime)
{
	int i;					//繰返し用変数.
	int num;			//添字
	double best;		//一時保存用
	int Ucr_fg = 0;
	for (i = 0, num = 0, best = cFitness[0]; i<Np; i++) {
		if (cFitness[i]<best) {
			best = cFitness[i];
			num = i;
			Ucr_fg = 1;
		}
	}
	if (Ucr_fg == 0) {
//		Ucr = 0.3;
//		Uf = 0.9;
//		getchar();
//		printf("失敗");
	}
	for (i = 0; i<d; i++) gBestVector[i] = cVect[num][i];
	gBestFitness = cFitness[num];
	gBestHistory[itime][gtime] = gBestFitness;
	//	printf("%20.10lf\n", gBestFitness);
}
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
	sprintf_s(filename, "DE_gBestHistory%04d%02d%02d%02d%02d_DE_NO%d_FUNC_NO%d_C%lf.txt",
		t_st->tm_year + 1900, t_st->tm_mon + 1,
		t_st->tm_mday, t_st->tm_hour, t_st->tm_min, DeAlgorithmNo, Func_No,C);
	fp = fopen(filename, "a");
	for (i = 0; i<MaxGrnration; i++) {
		for (j = 0; j<EXTIME; j++) {
			fprintf(fp, "%20.60lf\t", gBestHistory[j][i]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
	fclose(fp);

}
//------------------------------------------------------------
//集団の多様性の評価
//------------------------------------------------------------
void Calc_Diversity(int itime, int gtime)
{
	int i, j;					//繰返し用変数.
	double sum[D];		//合計
	double ave[D];		//平均
	double div;					//分散
								//初期化
	for (i = 0; i<d; i++) {
		sum[i] = 0.0;
		ave[i] = 0.0;
	}
	//合計
	for (i = 0; i<Np; i++) {
		for (j = 0; j<d; j++) {
			sum[j] += cVect[i][j];
		}
	}
	//平均
	for (i = 0; i<d; i++) {
		ave[i] = sum[i] / Np;
	}
	//分散
	for (i = 0, div = 0.0; i<Np; i++) {
		for (j = 0; j<d; j++) {
			div += (ave[j] - cVect[i][j])*(ave[j] - cVect[i][j]);
		}
	}
	pDiversity[itime][gtime] = div / Np;
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
	for (i = 0; i<MaxGrnration; i++) {
		for (j = 0; j<EXTIME; j++) {
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
		t_st->tm_mday, t_st->tm_hour, t_st->tm_min, DeAlgorithmNo, Func_No,C);
	fp = fopen(filename, "a");
	for (i = 0; i<EXTIME; i++) {
		fprintf(fp, "%6d %d\n", gTable[i], sRate[i]);
	}
	fclose(fp);
}
void Output_To_File4(void)
{
	int i, j;					//繰返し用変数.
	time_t timer;				//時間計測用
	struct tm *t_st;			//時間計測用
	time(&timer);				//時間の取得
	t_st = localtime(&timer);	//時間の変換
	sprintf_s(filename_4, "DE_CR%04d%02d%02d%02d%02d_DE_NO%d_FUNC_NO%dC%lf.txt",
		t_st->tm_year + 1900, t_st->tm_mon + 1,
		t_st->tm_mday, t_st->tm_hour, t_st->tm_min, DeAlgorithmNo, Func_No,C);
	fp_4 = fopen(filename_4, "a");


}


//------------------------------------------------------------
//メイン関数
//------------------------------------------------------------
int main(void)
{
	for (DeAlgorithmNo = 2; DeAlgorithmNo <= 5; DeAlgorithmNo++) {
		for (Func_No =1; Func_No <= 4; Func_No++) {
	//		if(DeAlgorithmNo ==5)	Output_To_File4();
			printf("DeAlgorithmNo=%d\nFunc_No=%d\n", DeAlgorithmNo, Func_No);
			int episode;
			int iteration;
			int pop;
			int J_rand = 0;
			int j = 0;
			init_genrand((unsigned)time(NULL));	//MTの初期化
			for (iteration = 0; iteration < EXTIME; iteration++) {
				episode = 0;
				Init_Vector();
				Evaluate_Init_Vector();
				Parameter_Initialization();
				while (episode < MaxGrnration) {
					Select_Elite_Vector(iteration, episode);
					Calc_Diversity(iteration, episode);
					for (pop = 0; pop < Np; pop++) {
						Parameter_Format();
						Select_pVector(pop);
						DE_Operation(pop, episode);
						Evaluate_New_Vector(pop);
					}
//					printf("%d\nUcr%lf Uf%lf\n", episode, Ucr, Uf);
					Compare_Vector();
//					printf("CR%lf\tF%lf\n", average_CR, average_F);
					episode++;
				}

				gTable[iteration] = episode;
				if (gBestFitness < Terminate)sRate[iteration] = 1;
				else sRate[iteration] = 0;
				gBestTable[iteration] = gBestFitness;
			}
			Output_To_File1();
			//Output_To_File2();
			//Output_To_File3();
//			if (DeAlgorithmNo == 5) fclose(fp_4);
		}

	}
//	printf("end\n");
//	getchar();
	return 0;
}