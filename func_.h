#pragma once
//-----------------------------------------------------------
int Np = NP;
int d = D;		
int Func_No = FUNC_NO;
double Range = RANGE;	
int MaxGrnration = MAXGRNRATION;
int DeAlgorithmNo = DE_ALGORITHM_NO;
//------------------------------------------------------------
double cVect[NP][D];	//時刻Tの個体（ベクトル）
double cVect_best[NP][D];
double cFitness[NP];	//時刻Tの個体の評価値
double pBestVector[NP][D];//時刻Tにおける各個体の最良解の履歴
double pBestFitness[NP];	//時刻Tにおける各個体の最良解評価値の履歴
double gBestVector[D];	//時刻Tにおける集団全体での最良解の履歴
double gBestFitness;	//時刻Tにおける集団全体での最良解評価値の履歴
double cBestFitness;	//時刻T-1における全体の最良値
double pVect1[D];	//親ベクトル１
double pVect2[D];	//親ベクトル２
double pVect3[D];	//親ベクトル３
double nVect[NP][D];	//時刻T+1の個体（ベクトル）
double nFitness[NP];	//時刻T+1の個体の評価値
double gBestHistory[EXTIME][MAXGRNRATION];	//各試行におけるgBestの履歴
double gBestHistoryMed[MAXGRNRATION][EXTIME];	//各試行におけるgBestの履歴
double gBestHistoryAll[4][6][MAXGRNRATION];
double gBestHistoryMedian[4][6][MAXGRNRATION];
double pDiversity[EXTIME][MAXGRNRATION];	//各試行における集団の多様性

double gBestTable[EXTIME];	//終了時点でのgBestの値
double MRATE = 0.3;	//突然変異率0.9
double CRATE = 0.5;	//交叉率
double C = 0.1;		//JADE のパラメータ
double Ucr = 0.5;	//JADEのパラメータ
double Uf = 0.5;	//JADEのパラメータ
double CR[NP];		//JADE,jDEのパラメータ
double F[NP];		//JADE,jDEのパラメータ
int episode;
int iteration;
int p_best;

void vDEParameter() {
if (Func_No == 1) {
	if (DeAlgorithmNo==1) {
		MRATE = 0.2;
		CRATE = 0.8;

	}
	else if(DeAlgorithmNo==2){
		MRATE = 0.4;
		CRATE = 0.8;
	}
	else if (DeAlgorithmNo == 3) {
		MRATE = 0.2;
		CRATE = 0.4;

	}
	else if (DeAlgorithmNo == 4) {
		MRATE = 0.5;
		CRATE = 0.4;

	}

}
else if (Func_No == 2) {
	if (DeAlgorithmNo == 1) {
		MRATE = 0.6;
		CRATE = 0.9;

	}
	else if (DeAlgorithmNo == 2) {
		MRATE = 0.9;
		CRATE = 0.9;
	}
	else if (DeAlgorithmNo == 3) {
		MRATE = 0.5;
		CRATE = 0.9;

	}
	else if (DeAlgorithmNo == 4) {
		MRATE = 0.8;
		CRATE = 0.8;
	}

}
else if (Func_No==3) {
	if (DeAlgorithmNo == 1) {
		MRATE = 0.2;
		CRATE = 0.2;

	}
	else if (DeAlgorithmNo == 2) {
		MRATE = 0.5;
		CRATE = 0.1;
	}
	else if (DeAlgorithmNo == 3) {
		MRATE = 0.1;
		CRATE = 0.1;

	}
	else if (DeAlgorithmNo == 4) {
		MRATE = 0.5;
		CRATE = 0.1;

	}

}
else if (Func_No == 4) {
	if (DeAlgorithmNo == 1) {
		MRATE = 0.3;
		CRATE = 0.7;

	}
	else if (DeAlgorithmNo == 2) {
		MRATE = 0.5;
		CRATE = 0.1;
	}
	else if (DeAlgorithmNo == 3) {
		MRATE = 0.3;
		CRATE = 0.3;

	}
	else if (DeAlgorithmNo == 4) {
		MRATE = 0.7;
		CRATE = 0.1;

	}
}
}


void vRange() {
if (Func_No == 1) {
	Range = 100;
	MaxGrnration = 1000;
}
else if (Func_No == 2) {
	Range = 30;
	MaxGrnration = 5000;//20000
}
else if (Func_No == 3) {
	Range = 5.12;
	MaxGrnration = 4000;
}
else if (Func_No == 4) {
	Range = 600;
	MaxGrnration = 4000;
}
else {
	printf("Rangeにあうパラメータが設定されていません\nFunc_No=%d", Func_No);
	exit(0);
}

}


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
for (i = 0; i < d; i++) {
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

for (i = 0; i < d - 1; i++) {
	sum += 100 * (x[i + 1] - x[i] * x[i])*(x[i + 1] - x[i] * x[i]) 
		+ (x[i] - 1.0)*(x[i] - 1.0);
}

// UNDXの文献の式に統一
// Rosenbrock関数 Star型
/*
for (i = 2; i<d; i++) {
sum += (100 * (x[0] - x[i] * x[i])*(x[0] - x[i] * x[i]) 
             + (x[1] - 1.0)*(x[1] - 1.0));
}
*/
return sum;
}
//------------------------------------------------------------
// F3 Rastrigin関数
//------------------------------------------------------------
double Rastrigin(double *x) {
int i;
double sum1 = 0.0, sum2 = 0.0;
for (i = 0; i < d; i++) {
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
for (i = 0; i < d; i++) {
	sum1 += x[i] * x[i];
	sum2 *= cos(x[i] / sqrt((i + 1)*1.0));
}
return ((sum1 / 4000) - sum2 + 1);
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
else {
	getchar();
	exit(0);
}
}

//------------------------------------------------------------
//初期個体群の生成
//------------------------------------------------------------
void Init_Vector(void)
{
int i, j;		//繰り返し用変数
double r;		//定義域用変数
vRange();
r = Range;
for (i = 0; i < Np; i++) {
	for (j = 0; j < d; j++) {
		cVect[i][j] = Range * (genrand_real1() * 2 - 1);
	}
}
//初期ベクトルをpBestVectorに保存する
for (i = 0; i < Np; i++) {
	for (j = 0; j < d; j++) {
		pBestVector[i][j] = cVect[i][j];
	}
}

}

//------------------------------------------------------------
// 初期集団の評価（目的関数値を適応度として利用）
//------------------------------------------------------------
void Evaluate_Init_Vector(void) {
int i;
for (i = 0; i < Np; i++) {
	cFitness[i] = Calc_Objective_Function(cVect[i]);
	//初期値を初期pBestとして保存
	pBestFitness[i] = cFitness[i];
}
}

void Select_pVector1(int pop1)
{
register int i;
int pop2, pop3;
do {
	pop2 = (int)Np*genrand_real1();
	pop3 = (int)Np*genrand_real1();
} while (pop1 == pop2 || pop1 == pop3 || pop2 == pop3);
for (i = 0; i < d; i++)pVect1[i] = cVect[pop2][i];
for (i = 0; i < d; i++)pVect2[i] = cVect[pop3][i];

}

//------------------------------------------------------------
//新規ベクトルの生成のための親ベクトルの選択
//------------------------------------------------------------
void Select_pVector(int pop1)
{
if (DeAlgorithmNo == 5 || DeAlgorithmNo == 6) {
	Select_pVector1(pop1);
}
else {
	register int i;
	int pop2, pop3, pop4;
	do {
		pop2 = (int)Np*genrand_real1();
		pop3 = (int)Np*genrand_real1();
		pop4 = (int)Np*genrand_real1();
	} while (pop1 == pop2 || pop1 == pop3 || pop1 == pop4 
		|| pop2 == pop3 || pop2 == pop4 || pop3 == pop4);
	for (i = 0; i < d; i++)pVect1[i] = cVect[pop2][i];
	for (i = 0; i < d; i++)pVect2[i] = cVect[pop3][i];
	for (i = 0; i < d; i++)pVect3[i] = cVect[pop4][i];
}
}

void First_time_only_Initialize() {
init_genrand((unsigned)time(NULL));	//MTの初期化
int j = 0, i = 0, k = 0, l = 0;
for (i = 0; i < 4; i++) {
	for (j = 0; j < 6; j++) {
		for (k = 0; k < MAXGRNRATION; k++) {
			gBestHistoryAll[i][j][k] = 0;
			gBestHistoryMedian[i][j][k] = 0;
			for (l = 0; l < EXTIME; l++) 	gBestHistoryMed[k][l] = 0;
		}

		stv[i][j] = 0;
		ave[i][j] = 0;
	}
}

}

//bestの配列のデータの更新
void best_Initialize() {
int i, j;
for (i = 0; i < MaxGrnration; i++) {
	for (j = 0; j < EXTIME; j++) {
		gBestHistory[j][i] = 0;
	}
}
}

//パラメータの初期化
void Parameter_Initialization() {
int i;
if (DeAlgorithmNo == 5) {
	for (i = 0; i < Np; i++) {
		F[i] = 0.5;
		CR[i] = 0.5;
	}
}else if(DeAlgorithmNo==6){
	for (i = 0; i < Np; i++) {
		F[i] = 0.5;
		CR[i] = 0.9;
	}
}
Uf = 0.5;
Ucr = 0.5;
}

void swap(double *a, double *b) {
double c;
c = *a;
*a = *b;
*b = c;
}

void vnVect(int arry_num1, int arry_num2) {
int i;
double arry_data[D];
for (i = 0; i < D; i++) {
	arry_data[i] = nVect[arry_num1][i];
	nVect[arry_num1][i] = nVect[arry_num2][i];
	nVect[arry_num2][i] = arry_data[i];
}
}

void vcVect(int arry_num1, int arry_num2) {
int i;
double arry_data[D];
for (i = 0; i < D; i++) {
	arry_data[i] = cVect[arry_num1][i];
	cVect[arry_num1][i] = cVect[arry_num2][i];
	cVect[arry_num2][i] = arry_data[i];
}
}

void vcVect_best(int arry_num1, int arry_num2) {
int i;
double arry_data[D];
for (i = 0; i < D; i++) {
	arry_data[i] = cVect_best[arry_num1][i];
	cVect_best[arry_num1][i] = cVect_best[arry_num2][i];
	cVect_best[arry_num2][i] = arry_data[i];
}
}

void  bubbleSort(double *N) {
if (DeAlgorithmNo == 5 || DeAlgorithmNo == 6) {
	int i, j;
	for (i = 0; i < Np; i++) {
		for (j = Np - 1; j > i; j--) {
			if (N[j] < N[j - 1]) {
				swap(&N[j], &N[j - 1]);
				swap(&cFitness[j], &cFitness[j - 1]);
				vnVect(j, j - 1);
				vcVect(j, j - 1);
				swap(&F[j], &F[j - 1]);
				swap(&CR[j], &CR[j - 1]);
			}
		}
	}
}
}

void  bubbleSort_cnter(double N[MAXGRNRATION][EXTIME],int Nom) {
int i, j;
for (i = 0; i < EXTIME; i++) {
	for (j = EXTIME - 1; j > i; j--) {
		if (N[Nom][j] < N[Nom][j - 1]) {
			swap(&N[Nom][j], &N[Nom][j - 1]);
		}
	}
}
}

void  Initialize_bubbleSort(double *N) {
int i, j;
for (i = 0; i < Np; i++) {
	for (j = Np - 1; j > i; j--) {
		if (N[j] < N[j - 1]) {
			swap(&N[j], &N[j - 1]);
			vcVect(j, j - 1);
		}
	}
}
}

//------------------------------------------------------------
//DEの操作
//------------------------------------------------------------
void DE_Operation(int i_Np, int g_GSIZE)
{
register int i;
int N = 0, L = 0;
int best_rand_N;
//DE/rand/1/exp
if (DeAlgorithmNo == 1) {
	for (i = 0; i < d; i++) nVect[i_Np][i] = cVect[i_Np][i];
	N = (int)(genrand_real2()*d);
	L = 0;
	do {
		nVect[i_Np][N] = pVect1[N] + MRATE * (pVect2[N] - pVect3[N]);
		if (nVect[i_Np][N] < -Range) 	nVect[i_Np][N] = pVect1[N] 
		+ genrand_real1() * (-Range - pVect1[N]);
		if (nVect[i_Np][N] > Range) 	nVect[i_Np][N] = pVect1[N] 
		+ genrand_real1() * (Range - pVect1[N]);
		N = (N + 1) % d;
		L++;
	} while (genrand_real1() < CRATE && L < d);
}
//DE/best/1/exp
else if (DeAlgorithmNo == 2) {
	for (i = 0; i < d; i++) nVect[i_Np][i] = cVect[i_Np][i];
	N = (int)(genrand_real2()*d);
	L = 0;
	do {
	 nVect[i_Np][N] = gBestVector[N] + MRATE * (pVect2[N] - pVect3[N]);
	 if (nVect[i_Np][N] < -Range) nVect[i_Np][N] = pVect1[N] 
		+ genrand_real1() * (-Range - pVect1[N]);
	 if (nVect[i_Np][N] > Range) nVect[i_Np][N] = pVect1[N] 
		+ genrand_real1() * (Range - pVect1[N]);
	 N = (N + 1) % d;
	 L++;
	} while (genrand_real1() < CRATE && L < d);
}
//DE/rand/1/bin
else if (DeAlgorithmNo == 3) {
	for (i = 0; i < d; i++) nVect[i_Np][i] = cVect[i_Np][i];
	N = (int)(genrand_real2()*d);
	for (L = 0; L < d; L++) {
	 if (L == 0 || genrand_real1() < CRATE) {
	 	nVect[i_Np][N] = pVect1[N] + MRATE * (pVect2[N] - pVect3[N]);
	 	if (nVect[i_Np][N] < -Range) nVect[i_Np][N] = pVect1[N] 
			+ genrand_real1() * (-Range - pVect1[N]);
	 	if (nVect[i_Np][N] > Range) nVect[i_Np][N] = pVect1[N] 
			+ genrand_real1() * (Range - pVect1[N]);
	 }
	 else {
	 	nVect[i_Np][N] = cVect[i_Np][N];
	 }
	 N = (N + 1) % d;
	}
}
//DE/best/1/bin
else if (DeAlgorithmNo == 4) {
	for (i = 0; i < d; i++) nVect[i_Np][i] = cVect[i_Np][i];
	N = (int)(genrand_real2()*d);
	for (L = 0; L < d; L++) {
	 if (L == 0 || genrand_real1() < CRATE) {
	 	nVect[i_Np][N] = gBestVector[N] + MRATE * (pVect2[N] - pVect3[N]);
	 	if (nVect[i_Np][N] < -Range)	nVect[i_Np][N] = pVect1[N] 
			+ genrand_real1() * (-Range - pVect1[N]);
	 	if (nVect[i_Np][N] > Range) nVect[i_Np][N] = pVect1[N] 
			+ genrand_real1() * (Range - pVect1[N]);
	 }
	 else {
	 	nVect[i_Np][N] = cVect[i_Np][N];
	 }
	 N = (N + 1) % d;
	}
}

//JADE
else if (DeAlgorithmNo == 5) {
	L = 0;
	for (i = 0; i < d; i++) nVect[i_Np][i] = cVect[i_Np][i];
	N = (int)(genrand_real1()*d);
	best_rand_N = (int)(genrand_real2()*p_best);
	do {
	 if (genrand_real1() < CR[i_Np] || L == N) {
	 	nVect[i_Np][N] = cVect[i_Np][N] + (F[i_Np])* (cVect[best_rand_N][N] - cVect[i_Np][N] 
			+ (F[i_Np])* (pVect1[N] - pVect2[N]));
	 	if (nVect[i_Np][N] < -Range) 	nVect[i_Np][N] = pVect1[N] 
		+ genrand_real1() * (-Range - pVect1[N]);
	 	if (nVect[i_Np][N] > Range) nVect[i_Np][N] = pVect1[N] 
		+ genrand_real1() * (Range - pVect1[N]);
	 }
	 else {
	 	nVect[i_Np][N] = cVect[i_Np][N];
	 }
	 N = (N + 1) % d;
	 L++;
 } while (L < d);

}
//jde
else if (DeAlgorithmNo == 6) {
	L = 0;
	for (i = 0; i < d; i++) nVect[i_Np][i] = cVect[i_Np][i];
	N = (int)(genrand_real1()*d);
	best_rand_N = (int)(genrand_real2()*p_best);
	do {
 	if (genrand_real1() < CR[i_Np] || L == N) {
 		nVect[i_Np][N] = cVect[i_Np][N] + (F[i_Np])* (cVect[best_rand_N][N] - cVect[i_Np][N])
			+ (F[i_Np])* (pVect1[N] - pVect2[N]);
 		if (nVect[i_Np][N] < -Range) 	nVect[i_Np][N] = pVect1[N] 
			+ genrand_real1() * (-Range - pVect1[N]);
 		if (nVect[i_Np][N] > Range) nVect[i_Np][N] = pVect1[N] 
			+ genrand_real1() * (Range - pVect1[N]);
 	}
 	else {
 		nVect[i_Np][N] = cVect[i_Np][N];
 	}
 	N = (N + 1) % d;
 	L++;
 } while (L < d);
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
for (i = 0; i < Np; i++) {
	//新しいベクトルが良ければ置き換え操作を行う
	if (nFitness[i] < cFitness[i]) {
		cFitness[i] = nFitness[i];
		for (j = 0; j < d; j++)	cVect[i][j] = nVect[i][j];
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
cBestFitness = gBestFitness;
for (i = 0, num = 0, best = cFitness[0]; i < Np; i++) {
	if (cFitness[i] < best) {
		best = cFitness[i];
		num = i;
	}
}
for (i = 0; i < d; i++) gBestVector[i] = cVect[num][i];
gBestFitness = cFitness[num];
gBestHistory[itime][gtime] = gBestFitness;
gBestHistoryMed[gtime][itime] = gBestFitness;
gBestHistoryAll[Func_No-1][DeAlgorithmNo-1][episode]+=gBestFitness;

//	printf("%20.10lf\n", gBestFitness);
}


//初期化
void Initialize() {
episode = 0;
Init_Vector();
Evaluate_Init_Vector();
Parameter_Initialization();
Initialize_bubbleSort(cFitness);
}

//	JADEのパラメータ更新
void New_parameter() {
int i;
double Sn = 0.0, Sf = 0.0, Sf2 = 0.0, Scr = 0.0;
double hoge_fp = 0;
for (i = 0; i < Np; i++) {
	if (nFitness[i] < cFitness[i]) {
		Sn += 1;
		Sf += F[i];
		Sf2 += F[i] * F[i];
		Scr += CR[i];

	}
}
if (Sn != 0) {
	Ucr = (1 - C)*Ucr + C * ((Scr / Sn));
	Uf = (1 - C)*Uf + C * (Sf2 / Sf);
}

}


void vJade_Parameter_Format() {
int i = 0;
double sum_F = 0;
for (i = 0; i < Np; i++) {
	do {
		F[i] = rand_cauchy(Uf, 0.1);
		if (F[i] > 1) F[i] = 1;
	} while (F[i] < 0);
	do {
		CR[i] = rand_normal(Ucr, 0.1);
		if (CR[i] > 1) CR[i] = 1;
	} while (CR[i] < 0);
}
}


void vjde_Parameter_Format() {
int i;
for (i = 0; i < Np; i++) {
	if (genrand_real2() < Tf) {
		F[i] = Fl + genrand_real2()*Fu;
	}
	if (genrand_real2() < Tcr) {
		CR[i] = genrand_real2();
	}
}
}

//パラメータの更新
void Parameter_Format(int iDe_nomber) {

if (iDe_nomber == 5) {
	New_parameter();
	vJade_Parameter_Format();
}
else if (iDe_nomber == 6) {
	vjde_Parameter_Format();

}
}
