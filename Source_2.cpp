//------------------------------------------------------------
// �W���I�ȍ����i���̃A���S���Y�������������v���O����
//------------------------------------------------------------
//------------------------------------------------------------
//�uerror C2065: 'M_PI' : ��`����Ă��Ȃ����ʎq�ł��v�ւ̑Ή�
//�@VisualStudio�̏ꍇ�͈ȉ��̃R�����g���O���ăR���p�C��

//best���~���ɕ��ёւ���100p%�̊m��(p=0.05�����5�ʂ܂�)�Ń����_����best��I�ԁD
//���Ȃ���΂Ȃ�Ȃ����Ƃ�best�̕��ёւ��Cbest�����100p%�ɂȂ�悤�ɂ���D

// �ȉ��̓����Z���k�c�C�X�^�[�Ŏg�p����
#define NP	100	//�ő�̐�
#define D		30				//�ő原�����i���̎������j
#define FUNC_NO		1					//�œK�����̎��
#define RANGE			5.12				//�œK�����̒�`��
#define MRATE			0.3				//�ˑR�ψٗ�0.9
#define CRATE			0.5				//������
#define MAXGRNRATION	150000		//�ő�J��Ԃ���
#define DE_ALGORITHM_NO	1	//DE�̃A���S���Y��
#define EXTIME			10		//���s��
#define Terminate		1.0e-11			//�I������
#define Fl 0.1
#define Fu 0.9
#define Tf 0.29
#define Tcr 0.24
#define P_BEST 0.05		//�W����0.05(5%)

//------------------------------------------------------------

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <vector>

// �ȉ��̓����Z���k�c�C�X�^�[�Ŏg�p����
#include "random.h"
#include "hoge.h"
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
//�W�c�̑��l���̕]��
//------------------------------------------------------------
void Calc_Diversity(int itime, int gtime)
{
	int i, j;					//�J�Ԃ��p�ϐ�.
	double sum[D];		//���v
	double ave[D];		//����
	double div;					//���U
								//������
	for (i = 0; i < d; i++) {
		sum[i] = 0.0;
		ave[i] = 0.0;
	}
	//���v
	for (i = 0; i < Np; i++) {
		for (j = 0; j < d; j++) {
			sum[j] += cVect[i][j];
		}
	}
	//����
	for (i = 0; i < d; i++) {
		ave[i] = sum[i] / Np;
	}
	//���U
	for (i = 0, div = 0.0; i < Np; i++) {
		for (j = 0; j < d; j++) {
			div += (ave[j] - cVect[i][j])*(ave[j] - cVect[i][j]);
		}
	}
	pDiversity[itime][gtime] = div / Np;
}
//------------------------------------------------------------
//�t�@�C���o��2�@���U�l�̐���
//------------------------------------------------------------
void Output_To_File2(void)
{
	int i, j;					//�J�Ԃ��p�ϐ�.
	FILE *fp;					//�t�@�C���|�C���^
	char filename[50];			//�t�@�C����
	time_t timer;				//���Ԍv���p
	struct tm *t_st;			//���Ԍv���p
	time(&timer);				//���Ԃ̎擾
	t_st = localtime(&timer);	//���Ԃ̕ϊ�
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
//�t�@�C���o��3�@�œK�𔭌�����C������
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
	int i, j;					//�J�Ԃ��p�ϐ�.
	double gBestHistorySum = 0;
	FILE *fp;					//�t�@�C���|�C���^
	char filename[100];			//�t�@�C����
	time_t timer;				//���Ԍv���p
	struct tm *t_st;			//���Ԍv���p
	time(&timer);				//���Ԃ̎擾
	t_st = localtime(&timer);	//���Ԃ̕ϊ�
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
	int i, j;					//�J�Ԃ��p�ϐ�.
	double gBestHistorySum = 0;
	FILE *fp;					//�t�@�C���|�C���^
	char filename[100];			//�t�@�C����
	time_t timer;				//���Ԍv���p
	struct tm *t_st;			//���Ԍv���p
	time(&timer);				//���Ԃ̎擾
	t_st = localtime(&timer);	//���Ԃ̕ϊ�
	sprintf_s(filename, "DE_Uf_best_DE_NO%d_FUNC_NO%d_NP%d_D%d_ave.txt",
		DeAlgorithmNo, Func_No, NP, D);
	fp = fopen(filename, "a");

	for (i = 0; i < MaxGrnration; i++) {
		for (j = 0; j < EXTIME; j++) {
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
	int pop;
	int J_rand = 0;
	int j = 0, i = 0;
	p_best = (int)(NP*P_BEST);
	for (Func_No = 1; Func_No <= 4; Func_No++) {
		for (DeAlgorithmNo = 1; DeAlgorithmNo <= 5; DeAlgorithmNo++) {
			//if(DeAlgorithmNo ==5)	Output_To_File4();
			printf("DeAlgorithmNo=%d\nFunc_No=%d\n", DeAlgorithmNo, Func_No);
			pop = 0, J_rand = 0, j = 0, i = 0;
			best_Initialize();
			init_genrand((unsigned)time(NULL));	//MT�̏�����
			for (iteration = 0; iteration < EXTIME; iteration++) {//���s��
				Initialize();		//������

				while (episode < MaxGrnration) {
					Select_Elite_Vector(iteration, episode);
					//					Calc_Diversity(iteration, episode);
					for (pop = 0; pop < Np; pop++) {
						Select_pVector(pop);
						DE_Operation(pop, episode);
						Evaluate_New_Vector(pop);
						//printf("%d Ucr%lf Uf_best%lf Uf_rand%lf \n", episode, Ucr, Uf_best,Uf_rand);
					}

					if (DeAlgorithmNo == 6) {
						Uf_best_History[iteration][episode] = Uf_best;
						Uf_rand_History[iteration][episode] = Uf_rand;
						printf("Uf_best=%5.2f\nUf_ranf=%5.2f\n", Uf_best, Uf_rand);

					}

					Parameter_Format(DeAlgorithmNo);
					Compare_Vector();
					bubbleSort(nFitness);

					episode++;
					if (Terminate > gBestFitness) break;
				}

				gTable[iteration] = episode;
				if (gBestFitness < Terminate)sRate[iteration] = 1;
				else sRate[iteration] = 0;
				gBestTable[iteration] = gBestFitness;
			}
			Output_To_File1();
			//Output_To_File2();
			//Output_To_File3();
/*			if (DeAlgorithmNo == 5 || DeAlgorithmNo==6|| DeAlgorithmNo==7) {
				Output_To_File4();
				Output_To_File5();
			}
*/
		}

	}

	return 0;
}


/*
for (i = 0; i < Np; i++) {
for (j = 0; j < d; j++) {
printf("�ω��OnVect[%d][%d]=%0.10f\n", i, j, nVect[i][j]);
}
printf("\n");
}
printf("\n\n");
for (i = 0; i < Np; i++) {
printf("�ω��OnFitness[%d]=%0.10f\n", i, nFitness[i]);
}
*/

//QuickSort();

/*
printf("\n");
for (i = 0; i < Np; i++) {
printf("�ω���nFitness[%d]=%0.10f\n", i, nFitness[i]);
}
printf("\n");
for (i = 0; i < Np; i++) {
for (j = 0; j < d; j++) {
printf("�ω���nVect[%d][%d]=%0.10f\n", i, j, nVect[i][j]);
}
printf("\n");
}
printf("\n\n");
*/

//					printf("Uf=%5.2f \n Ucr=%5.2f\n", Uf, Ucr);




/*
for (i = 0; i < Np; i++) {
printf("�ω���nFitness[%d]=%0.10f\n", i, nFitness[i]);
}
printf("\n");
for (i = 0; i < Np; i++) {
for (j = 0; j < d; j++) {
printf("�ω���nVect[%d][%d]=%0.10f\n", i, j, nVect[i][j]);
}
printf("\n");
}
*/
//printf("CR%lf\tF%lf\n", average_CR, average_F);