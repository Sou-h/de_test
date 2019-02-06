#define NP	100	//�ő�̐�
#define D	30	//�ő原�����i���̎������j
#define FUNC_NO	1	//�œK�����̎��
#define RANGE	5.12	//�œK�����̒�`��
#define MAXGRNRATION	5000	//�ő�J��Ԃ���
#define DE_ALGORITHM_NO	1	//DE�̃A���S���Y��
#define EXTIME	100	//���s��
#define Terminate	1.0e-11	//�I������
#define Fl 0.1	//jDE�̃p�����[�^
#define Fu 0.9	//jDE�̃p�����[�^
#define Tf 0.1	//jDE�̃p�����[�^
#define Tcr 0.1	//jDE�̃p�����[�^
#define P_BEST 0.10	//�W����0.05(5%)
#define P_BEST_COUNT	6
#define _USE_MATH_DEFINES

double sum_data = 0;	//���v
double ave[4][6];	//����
double stv[4][6];	//�W���΍�


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>

#include "random.h"
#include "func_.h"
// �ȉ��̓����Z���k�c�C�X�^�[�Ŏg�p����
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
//�t�@�C���o��1�@�ŗǒl�̐���
//------------------------------------------------------------
void Output_To_File1(void)
{
int i, j;	//�J�Ԃ��p�ϐ�.
FILE *fp;	//�t�@�C���|�C���^
char filename[80];	//�t�@�C����
time_t timer;		//���Ԍv���p
struct tm *t_st;	//���Ԍv���p
time(&timer);		//���Ԃ̎擾
t_st = localtime(&timer);//���Ԃ̕ϊ�
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


//���U�̒��o
void Calc_Diversity2()
{
int i, j;					//�J�Ԃ��p�ϐ�.
double div = 0;	//���U

//���v
for (i = 0; i < EXTIME; i++) {
	sum_data += gBestHistory[i][MaxGrnration - 1];
}
if (sum_data != 0) {
	//����
	ave[Func_No - 1][DeAlgorithmNo - 1] = sum_data / (EXTIME);

	//���U
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
//���ϒl�̏o��
void Output_To_File5(void)
{
int i, j;	//�J�Ԃ��p�ϐ�.
FILE *fp;	//�t�@�C���|�C���^
char filename[100];	//�t�@�C����
time_t timer;		//���Ԍv���p
struct tm *t_st;	//���Ԍv���p
time(&timer);		//���Ԃ̎擾
t_st = localtime(&timer);	//���Ԃ̕ϊ�
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
//�����l�̏o��
void Output_To_File6(void)
{
int i, j;	//�J�Ԃ��p�ϐ�.
FILE *fp;	//�t�@�C���|�C���^
char filename[100];	//�t�@�C����
time_t timer;		//���Ԍv���p
struct tm *t_st;	//���Ԍv���p
time(&timer);		//���Ԃ̎擾
t_st = localtime(&timer);	//���Ԃ̕ϊ�
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
//���ϒl�ƕW���΍��̏o��
void Output_To_File7(void)
{
int  i;		//�J�Ԃ��p�ϐ�.
FILE *fp;	//�t�@�C���|�C���^
char filename[120];	//�t�@�C����
time_t timer;		//���Ԍv���p
struct tm *t_st;	//���Ԍv���p
time(&timer);		//���Ԃ̎擾
t_st = localtime(&timer);	//���Ԃ̕ϊ�
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
//���C���֐�
//------------------------------------------------------------
int main(void)
{
int pop;
p_best = (int)(NP*P_BEST); //���p%�̐ݒ�JADE,jDE�Ŏg�p
First_time_only_Initialize();
for (Func_No = 1; Func_No <= 1; Func_No++) {
	for (DeAlgorithmNo = 5; DeAlgorithmNo <= 5; DeAlgorithmNo++) {
		printf("DeAlgorithmNo=%d\nFunc_No=%d\n", DeAlgorithmNo, Func_No);
		vDEParameter();//DE�̃p�����[�^�̐ݒ�
		best_Initialize();//gBestHistory�̏�����
		for (iteration = 0; iteration < EXTIME; iteration++) {//���s��
			Initialize();	//������
			while (episode < MaxGrnration) {//���㐔
				Select_Elite_Vector(iteration, episode);//�ŗǒl�̊i�[
				for (pop = 0; pop < Np; pop++) {
					Select_pVector(pop);//p1,p2,p3�̑I��
					DE_Operation(pop, episode);//����,�ˑR�ψ�
					Evaluate_New_Vector(pop);//�e�q�̔�r
				}
				//JADE,jDE�̍ۂɃp�����[�^�̍X�V
				Parameter_Format(DeAlgorithmNo);
				Compare_Vector();
				//JADE,jDE�̍ۂɃ\�[�g�����s
				bubbleSort(nFitness);
				episode++;
				if (Terminate > gBestFitness) break; //�ڕW�l�Ƃ̔�r
			}
		}
		Calc_Diversity2();//���U�̒��o
		Output_To_File1();//�ŗǒl�̏o��
	}
	Output_To_File5();//���ϒl�̏o��
	Output_To_File6();//�����l�̏o��
	Output_To_File7();//���ϒl�ƕW���΍��̏o��
}
return 0;
}
