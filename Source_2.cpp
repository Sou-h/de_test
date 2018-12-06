/*

���݁C���삳���Ă���̂�DeAlgorithmNo==6
Best�̒l�̍X�V���ƑO��̌̂Ƃ̍X�V���Ńp�����[�^�𕪗����l����D

Schwefel function�œ���s�ǁ@�œK���𒴂����ŏ��l���Ƃ�

*/


//------------------------------------------------------------
// �W���I�ȍ����i���̃A���S���Y�������������v���O����
//------------------------------------------------------------
//------------------------------------------------------------
//�uerror C2065: 'M_PI' : ��`����Ă��Ȃ����ʎq�ł��v�ւ̑Ή�
//�@VisualStudio�̏ꍇ�͈ȉ��̃R�����g���O���ăR���p�C��
//------------------------------------------------------------
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
// �ȉ��̓����Z���k�c�C�X�^�[�Ŏg�p����
#include "random.h"
//------------------------------------------------------------
// �����Z���k�c�C�X�^�[�Ŏg�p����֐��̊T��
// genrand_int32() //�����Ȃ�32�r�b�g������
// genrand_int31() //�����Ȃ�31�r�b�g������
// genrand_real1() //��l������[0,1] (32�r�b�g���x)
// genrand_real2() //��l������[0,1) (32�r�b�g���x)
// genrand_real3() //��l������(0,1) (32�r�b�g���x)
// genrand_res53() //��l������[0,1) (53�r�b�g���x)

#define _CRT_SECURE_NO_DEPRECATE
#pragma warning(disable:4996)
//------------------------------------------------------------
// �ݒ�p�����[�^ �������ɂ��̕���������������
//------------------------------------------------------------

#define NP	100				//�ő�̐�
#define D		10				//�ő原�����i���̎������j
#define FUNC_NO		1					//�œK�����̎��
#define MRATE			0.8				//�ˑR�ψٗ�0.9
#define CRATE			0.5				//������
#define MAXGRNRATION	1000		//�ő�J��Ԃ���
#define DE_ALGORITHM_NO	1	//DE�̃A���S���Y��
#define EXTIME			10		//���s��
#define Terminate		1.0e-20			//�I������

//-----------------------------------------------------------
int Np = NP;				//�ő�̐�
int d = D;					//�ő原�����i���̎������j
int Func_No = FUNC_NO;		//�œK�����̎��
double Range ;						//�œK�����̒�`��
int MaxGrnration = MAXGRNRATION;			//�ő�J��Ԃ���
int DeAlgorithmNo = DE_ALGORITHM_NO;	//DE�̃A���S���Y��
//------------------------------------------------------------
double cVect[NP][D];				//����T�̌́i�x�N�g���j
double cFitness[NP];				//����T�̌̂̕]���l
double pBestVector[NP][D];	//����T�ɂ�����e�̂̍ŗǉ��̗���
double pBestFitness[NP];		//����T�ɂ�����e�̂̍ŗǉ�]���l�̗���
double gBestVector[D];			//����T�ɂ�����W�c�S�̂ł̍ŗǉ��̗���
double gBestFitness;				//����T�ɂ�����W�c�S�̂ł̍ŗǉ�]���l�̗���
double cBestFitness;				//����T-1�ɂ�����S�̂̍ŗǒl
double pVect1[D];					//�e�x�N�g���P
double pVect2[D];					//�e�x�N�g���Q
double pVect3[D];					//�e�x�N�g���R
double nVect[NP][D];			//����T+1�̌́i�x�N�g���j
double nFitness[NP];				//����T+1�̌̂̕]���l
double gBestHistory[EXTIME][MAXGRNRATION];	//�e���s�ɂ�����gBest�̗���
double pDiversity[EXTIME][MAXGRNRATION];		//�e���s�ɂ�����W�c�̑��l��
int gTable[EXTIME];					//�œK�𔭌����̐��㐔
int sRate[EXTIME];					//�œK���̔�����
double gBestTable[EXTIME];		//�I�����_�ł�gBest�̒l
double C = 0.1;							//JADE �̃p�����[�^
double Ucr = 0.5;
double Uf = 0.5;
double Uf_best = 0.9,Uf_rand=0.9;
double CR[NP];
double F[NP];
double F_best[NP], F_rand[NP];
double Uf_best_History[EXTIME][MAXGRNRATION], Uf_rand_History[EXTIME][MAXGRNRATION];
double average_F = 0, average_CR = 0, average_Sn = 0;
int episode;
int iteration;


//------------------------------------------------------------
//------------------------------------------------------------
// [min, max]�̈�l����
//------------------------------------------------------------
double uniform(double min, double max) {
	return min + (max - min)*genrand_real1();
}
//------------------------------------------------------------
// �e�X�g�֐�
// F1 Sphere�֐�
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
// F2 Rosenbrock�֐� Chain�^
// [-2.048, +2.048] x=(1,..,1) F2=0 
//------------------------------------------------------------
double Rosenbrock(double *x) {
	int i;
	double sum = 0.0;

	for (i = 0; i<d - 1; i++){
	sum += 100 * (x[i + 1] - x[i] * x[i])*(x[i + 1] - x[i] * x[i]) + (x[i] - 1.0)*(x[i] - 1.0);
	}

	// UNDX�̕����̎��ɓ���
	// Rosenbrock�֐� Star�^
	/*
	for (i = 2; i<d; i++) {
		sum += (100 * (x[0] - x[i] * x[i])*(x[0] - x[i] * x[i]) + (x[1] - 1.0)*(x[1] - 1.0));
	}
	*/
	return sum;
}
//------------------------------------------------------------
// F3 Rastrigin�֐�
//------------------------------------------------------------
double Rastrigin(double *x) {
	int i;
	double sum1 = 0.0, sum2 = 0.0;
	for (i = 0; i<d; i++) {
		sum1 += x[i] * x[i] - 10 * cos(2 * M_PI*x[i]);
	}
	return (10 * d + sum1);
}
//------------------------------------------------------------
// F4 Griewank�֐�
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
// F5 Ackley�֐�
//------------------------------------------------------------
double Ackley(double *x) {
	int i;
	double sum1, sum2;
	for (i = 0, sum1 = 0.0; i<d; i++)sum1 += x[i] * x[i];
	for (i = 0, sum2 = 0.0; i<d; i++)sum2 += cos(2.0*M_PI*d);
	return (20.0 + M_E - 20.0 * exp(-0.2*sqrt(sum1 / d)) - exp(sum2 / d));
}
//------------------------------------------------------------
// F6 Schwefel�֐�
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
// �ړI�֐��l���v�Z
//------------------------------------------------------------
double Calc_Objective_Function(double *x)
{
	if (Func_No == 1)return  Sphere(x);
	else if (Func_No == 2)return  Rosenbrock(x);
	else if (Func_No == 3)return  Rastrigin(x);
	else if (Func_No == 4)return  Griewank(x);
	else if (Func_No == 5)return  Ackley(x);
	else if (Func_No == 6)return  Schwefel(x);
	else {
		getchar();
		exit(0);
	}
}
void vRange() {
	if (Func_No == 1)		Range = 100;
	else if (Func_No == 2) Range = 30;
	else if (Func_No == 3) Range = 5.12;
	else if (Func_No == 4) Range = 600;
	else if (Func_No == 5) Range = 32.768;
	else if (Func_No == 6) Range = 500;
	else if (Func_No == 9) Range = 100;
	else if (Func_No == 10) Range = 5.12;
	else {
		printf("Range�ɂ����p�����[�^���ݒ肳��Ă��܂���\nFunc_No=%d", Func_No);
		getchar();
		exit(0);
	}

}
//------------------------------------------------------------
//�����̌Q�̐���
//------------------------------------------------------------
void Init_Vector(void)
{
	int i, j;		//�J��Ԃ��p�ϐ�
	double r;		//��`��p�ϐ�
	vRange();
	r = Range;
	for (i = 0; i<Np; i++) {
		for (j = 0; j<d; j++) {
			cVect[i][j] = Range * (genrand_real1() * 2 - 1);
		}
	}
	//�����x�N�g����pBestVector�ɕۑ�����
	for (i = 0; i<Np; i++) {
		for (j = 0; j<d; j++) {
			pBestVector[i][j] = cVect[i][j];
		}
	}

}


//------------------------------------------------------------
// �����W�c�̕]���i�ړI�֐��l��K���x�Ƃ��ė��p�j
//------------------------------------------------------------
void Evaluate_Init_Vector(void) {
	int i;
	for (i = 0; i<Np; i++) {
		cFitness[i] = Calc_Objective_Function(cVect[i]);
		//�����l������pBest�Ƃ��ĕۑ�
		pBestFitness[i] = cFitness[i];
	}
}
//------------------------------------------------------------
//�V�K�x�N�g���̐����̂��߂̐e�x�N�g���̑I��
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

//�p�����[�^�̍X�V	JADE �쐬��
void New_parameter() {
	int i;
	double Sn = 0.0, Sf = 0.0, Sf2 = 0.0, Scr = 0.0;
	double hoge_fp = 0;
	for (i = 0; i < Np; i++) {
		if (nFitness[i] > cFitness[i]) {
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


//�p�����[�^�̍X�V	JADE �쐬��
void New_parameter_2() {
	int i;
	double Sn = 0.0,Scr = 0.0,Sn_best=0;
	double Sf_rand = 0,Sf_rand_2=0;
	double Sf_best = 0, Sf_best_2=0;
	bool hoge_fp = 0;
	for (i = 0; i < Np; i++) {
		if (cBestFitness > cFitness[i]) {
			Sf_best += F_best[i];
			Sf_best_2 += F_best[i] * F_best[i];
			hoge_fp = 1;
			Sn_best += 1;
		}
		if (nFitness[i] > cFitness[i]) {
			Sn += 1;
			Sf_rand += F_rand[i];
			Sf_rand_2 += F_rand[i] * F_rand[i];
			Scr += CR[i];
		}
	}

	if (Sn != 0) {
		Ucr = (1 - C)*Ucr + C * ((Scr / Sn));
		Uf_rand = (1 - C)*Uf_rand + C * (Sf_rand_2 / Sf_rand);
	}


	if (hoge_fp != 0) {
		Uf_best = (1 - C)*Uf_best + C * (Sf_best_2 / Sf_best);
	}
	if (Uf_rand > 0.9) {
		Uf_rand = 0.9;
	}
	if (Uf_best > 0.9) {
		Uf_best = 0.9;
	}
	else if (Uf_best < 0.1) {
		Uf_best = 0.1;
	}

}

//�p�����[�^�̍X�V
void Parameter_Format() {
	int i = 0;
	double sum_F = 0, sum_CR = 0;
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

}

void Parameter_Format_2() {
	int i = 0;
	double sum_F = 0, sum_CR = 0;
	for (i = 0; i < Np; i++) {

		F_best[i] = rand_cauchy(Uf_best, 1);
		if (F_best[i] > 0.9) {
			F_best[i] = 0.9;
		}
		else if (F_best[i] < 0.1) {
			F_best[i] = 0.1;
			F_rand[i] = rand_cauchy(Uf_rand, 1);
			if (F_rand[i] > 0.9) {
				F_rand[i] = 0.9;
			}
			else 	if (F_rand[i] < 0.1) {
				F_rand[i] = 0.1;
			}

			CR[i] = rand_normal(Ucr, 0.2);
			if (CR[i] > 0.9) {
				CR[i] = 0.9;
			}
			else if (CR[i] < 0.1) {
				CR[i] = 0.1;
			}
		}
	}

}

//�p�����[�^�̏�����
void Parameter_Initialization() {
	int i;
	for (i = 0; i < Np; i++) {
		F[i] = 0.85;
		CR[i] = 0.5;
		F_best[i] = 0.85;
		F_rand[i] = 0.5;
	}
	Uf = 0.85;
	Ucr = 0.5;
	Uf_rand = 0.5;
	Uf_best = 0.5;
}

//------------------------------------------------------------
//DE�̑���
//------------------------------------------------------------
void DE_Operation(int i_Np, int g_GSIZE)
{
	register int i;
	int N = 0, L = 0;
	int d_rand = 0;

	//DE/rand/1/exp
	if (DeAlgorithmNo == 1) {
		for (i = 0; i<d; i++) nVect[i_Np][i] = cVect[i_Np][i];
		N = (int)(genrand_real1()*d);
		L = 0;
		do {
			nVect[i_Np][N] = pVect1[N] + MRATE * (pVect2[N] - pVect3[N]);
			if (nVect[i_Np][N] < -Range) 	nVect[i_Np][N] = pVect1[N] + genrand_real1() * (-Range - pVect1[N]);
			if (nVect[i_Np][N] > Range) 	nVect[i_Np][N] = pVect1[N] + genrand_real1() * (Range - pVect1[N]);
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
		} while (genrand_real1() < CRATE && L < d);
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

	//JADE	
	else if (DeAlgorithmNo == 5) {

		L = 0;
		for (i = 0; i < d; i++) nVect[i_Np][i] = cVect[i_Np][i];
		N = (int)(genrand_real1()*d);
		do {
			if (genrand_real1()<CR[i_Np] || L == 0) {
				nVect[i_Np][N] = cVect[i_Np][N] + F[i_Np] * (gBestVector[N] - cVect[i_Np][N]) + F[i_Np] * (pVect1[N] - pVect2[N]);
				if (nVect[i_Np][N] < -Range) nVect[i_Np][N] = genrand_real1() * (-Range);
				if (nVect[i_Np][N] > Range) nVect[i_Np][N] =  genrand_real1() * Range ;
			}
			else {
				nVect[i_Np][N] = cVect[i_Np][N];
			}
			N = (N + 1) % d;
			L++;
		} while (L < d);

	}

	//JADE 2 
	else if (DeAlgorithmNo == 6) {
		L = 0;
		for (i = 0; i < d; i++) nVect[i_Np][i] = cVect[i_Np][i];
		N = (int)(genrand_real1()*d);
		N = L;
		d_rand=(int)(genrand_real1()*d);
		do {
			if (genrand_real1() < CR[i_Np] || L == 0) {
				nVect[i_Np][N] = pVect1[N] + (F_best[i_Np]) * (gBestVector[N] - cVect[i_Np][N]) + (F_rand[i_Np]) * (pVect3[N] - pVect2[N]);
//				nVect[i_Np][N] = pVect1[i_Np] + (F_best[i_Np]) * (gBestVector[i_Np] - pVect1[i_Np]) + (1 - F_rand[i_Np]) * (pVect1[i_Np] - pVect2[i_Np]);

				if (nVect[i_Np][N] < -Range) 	nVect[i_Np][N] =  genrand_real3() * (-Range);
				if (nVect[i_Np][N] > Range)	nVect[i_Np][N] = genrand_real3() * (Range );
			}
			else {
//				nVect[i_Np][d_rand] = cVect[i_Np][d_rand];
			}
			N = L;
			L++;

		} while (L < d);

	}

	else exit(0);
}

//------------------------------------------------------------
//�V�����x�N�g���̕]��
//------------------------------------------------------------
void Evaluate_New_Vector(int pop1)
{
	nFitness[pop1] = Calc_Objective_Function(nVect[pop1]);
}
//------------------------------------------------------------
//�x�N�g���̔�r
//------------------------------------------------------------
void Compare_Vector(void)
{
	int i, j;				//�J�Ԃ��p�ϐ�.
	for (i = 0; i<Np; i++) {
		//�V�����x�N�g�����ǂ���Βu������������s��
		if (nFitness[i] < cFitness[i]) {
			cFitness[i] = nFitness[i];
			for (j = 0; j<d; j++)cVect[i][j] = nVect[i][j];
		}
		else continue;
	}
}
//------------------------------------------------------------
//�G���[�g�I��
//------------------------------------------------------------
void Select_Elite_Vector(int itime, int gtime)
{
	int i;					//�J�Ԃ��p�ϐ�.
	int num;			//�Y��
	double best;		//�ꎞ�ۑ��p
	cBestFitness = gBestFitness;
	for (i = 0, num = 0, best = cFitness[0]; i<Np; i++) {
		if (cFitness[i]<best) {
			best = cFitness[i];
			num = i;
		}

	}

	for (i = 0; i<d; i++) gBestVector[i] = cVect[num][i];
	gBestFitness = cFitness[num];
	gBestHistory[itime][gtime] = gBestFitness;

}
//------------------------------------------------------------
//�t�@�C���o��1�@�ŗǒl�̐���
//------------------------------------------------------------
void Output_To_File1(void)
{
	int i, j;					//�J�Ԃ��p�ϐ�.
	FILE *fp;					//�t�@�C���|�C���^
	char filename[80];			//�t�@�C����
	time_t timer;				//���Ԍv���p
	struct tm *t_st;			//���Ԍv���p
	time(&timer);				//���Ԃ̎擾
	t_st = localtime(&timer);	//���Ԃ̕ϊ�
	sprintf_s(filename, "DE_gBestHistory%04d%02d%02d%02d%02d_DE_No%d_FUNC_No%d.txt",
		t_st->tm_year + 1900, t_st->tm_mon + 1,
		t_st->tm_mday, t_st->tm_hour, t_st->tm_min, DeAlgorithmNo, Func_No);
	fp = fopen(filename, "a");
	for (i = 0; i<MaxGrnration; i++) {
		for (j = 0; j < EXTIME; j++) {
			fprintf(fp, "%20.20lf\t", gBestHistory[j][i]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);

}


void Output_To_File4(void)
{
	int i, j;					//�J�Ԃ��p�ϐ�.
	FILE *fp;					//�t�@�C���|�C���^
	char filename[100];			//�t�@�C����
	time_t timer;				//���Ԍv���p
	struct tm *t_st;			//���Ԍv���p
	time(&timer);				//���Ԃ̎擾
	t_st = localtime(&timer);	//���Ԃ̕ϊ�
	sprintf_s(filename, "DE_Uf_rand_DE_NO%d_FUNC_NO%d_NP%d_D%d_ave.txt",
		DeAlgorithmNo, Func_No, NP, D);
	fp = fopen(filename, "a");

	for (i = 0; i <MaxGrnration; i++) {
		for (j = 0; j <EXTIME; j++) {
			fprintf(fp, "%.5lf\t", Uf_rand_History[j][i]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}

void Output_To_File5(void)
{
	int i, j;					//�J�Ԃ��p�ϐ�.
	FILE *fp;					//�t�@�C���|�C���^
	char filename[100];			//�t�@�C����
	time_t timer;				//���Ԍv���p
	struct tm *t_st;			//���Ԍv���p
	time(&timer);				//���Ԃ̎擾
	t_st = localtime(&timer);	//���Ԃ̕ϊ�
	sprintf_s(filename, "DE_Uf_best_DE_NO%d_FUNC_NO%d_NP%d_D%d_ave.txt",
		DeAlgorithmNo, Func_No, NP, D);
	fp = fopen(filename, "a");

	for (i = 0; i <MaxGrnration; i++) {
		for (j = 0; j <EXTIME; j++) {
			fprintf(fp, "%.5lf\t", Uf_best_History[j][i]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}


//------------------------------------------------------------
//���C���֐�
//------------------------------------------------------------
int main(void)
{
	for (Func_No = 1; Func_No <= 4; Func_No++) {
		for (DeAlgorithmNo = 6; DeAlgorithmNo <= 6; DeAlgorithmNo++) {
			printf("DeAlgorithmNo=%d\nFunc_No=%d\n", DeAlgorithmNo, Func_No);
			int pop;
			int J_rand = 0;
			int j = 0,i=0;
			for (i = 0; i < MaxGrnration; i++) {
				for (j = 0; j < EXTIME; j++) {
					gBestHistory[j][i] = 0;
					Uf_best_History[j][i] = 0;
					Uf_rand_History[j][i] = 0;
				}
			}
			init_genrand((unsigned)time(NULL));	//MT�̏�����
			for (iteration = 0; iteration < EXTIME; iteration++) {
				episode = 0;
				Init_Vector();
				Evaluate_Init_Vector();
				Parameter_Initialization();
				while (episode < MaxGrnration) {
					Select_Elite_Vector(iteration, episode);
					for (pop = 0; pop < Np; pop++) {
						Select_pVector(pop);
						DE_Operation(pop, episode);
						Evaluate_New_Vector(pop);
					}
					Uf_best_History[iteration][episode] = Uf_best;
					Uf_rand_History[iteration][episode] = Uf_rand;
					Compare_Vector();

					if(DeAlgorithmNo==5){ //�p�����[�^�̍X�V
						Parameter_Format();
						New_parameter();
					}else if (DeAlgorithmNo == 6) {
						Parameter_Format_2();
						New_parameter_2();
					}

					episode++;
					if (Terminate > gBestFitness) break;	//�I�������𖞂����ΏI��
				}
				gTable[iteration] = episode;
				if (gBestFitness < Terminate)sRate[iteration] = 1;
				else sRate[iteration] = 0;
				gBestTable[iteration] = gBestFitness;
			}
			Output_To_File1();
			if (DeAlgorithmNo ==6) {
				Output_To_File4();
				Output_To_File5();
			}

		}

	}
	return 0;
}