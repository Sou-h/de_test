#pragma once

//-----------------------------------------------------------
int Np = NP;				//�ő�̐�
int d = D;		//�ő原�����i���̎������j
int Func_No = FUNC_NO;					//�œK�����̎��
double Range = RANGE;				//�œK�����̒�`��
int MaxGrnration = MAXGRNRATION;		//�ő�J��Ԃ���
int DeAlgorithmNo = DE_ALGORITHM_NO;	//DE�̃A���S���Y��
										//------------------------------------------------------------
double cVect[NP][D];				//����T�̌́i�x�N�g���j
double cFitness[NP];				//����T�̌̂̕]���l
double pBestVector[NP][D];	//����T�ɂ�����e�̂̍ŗǉ��̗���
double pBestFitness[NP];		//����T�ɂ�����e�̂̍ŗǉ�]���l�̗���
double gBestVector[D];			//����T�ɂ�����W�c�S�̂ł̍ŗǉ��̗���
double gBestFitness;				//����T�ɂ�����W�c�S�̂ł̍ŗǉ�]���l�̗���
double cBestFitness;					//����T-1�ɂ�����S�̂̍ŗǒl
double pVect1[D];					//�e�x�N�g���P
double pVect2[D];					//�e�x�N�g���Q
double pVect3[D];					//�e�x�N�g���R
double nVect[NP][D];				//����T+1�̌́i�x�N�g���j
double nFitness[NP];				//����T+1�̌̂̕]���l
double gBestHistory[EXTIME][MAXGRNRATION];	//�e���s�ɂ�����gBest�̗���
double pDiversity[EXTIME][MAXGRNRATION];		//�e���s�ɂ�����W�c�̑��l��
int gTable[EXTIME];					//�œK�𔭌����̐��㐔
int sRate[EXTIME];						//�œK���̔�����
double gBestTable[EXTIME];		//�I�����_�ł�gBest�̒l
double C = 0.1;							//JADE �̃p�����[�^
double Ucr = 0.5;
double Uf = 0.5;
double Uf_best = 0.9, Uf_rand = 0.9;
//double Sf_best[NP],Sf_rand[NP];
//double Scr[NP];
double CR[NP];
double F[NP];
double F_best[NP], F_rand[NP];
double Uf_best_History[EXTIME][MAXGRNRATION], Uf_rand_History[EXTIME][MAXGRNRATION];
double average_F = 0, average_CR = 0, average_Sn = 0;
FILE *fp_4;					//�t�@�C���|�C���^
char filename_4[100];			//�t�@�C����
int No_best_sum = 1;
int episode;
int iteration;
int p_best;

void vRange() {
	if (Func_No == 1) {
		Range = 100;
		MaxGrnration = 1000;
	}
	else if (Func_No == 2) {
		Range = 30;
		MaxGrnration = 20000;
	}
	else if (Func_No == 3) {
		Range = 5.12;
		MaxGrnration = 4000;
	}
	else if (Func_No == 4) {
		Range = 600;
		MaxGrnration = 5000;
	}
	else if (Func_No == 5) {
		Range = 32.768;
		MaxGrnration = 1000;
	}
	else if (Func_No == 6) {
		Range = 500;
		MaxGrnration = 1000;
	}
	else if (Func_No == 9) {
		Range = 100;
		MaxGrnration = 1000;
	}
	else if (Func_No == 10) {
		Range = 5.12;
		MaxGrnration = 1000;
	}
	else {
		printf("Range�ɂ����p�����[�^���ݒ肳��Ă��܂���\nFunc_No=%d", Func_No);
		exit(0);
	}

}


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

	for (i = 0; i<d - 1; i++) {
		sum += 100 * (x[i + 1] - x[i] * x[i])*(x[i + 1] - x[i] * x[i]) + (x[i] - 1.0)*(x[i] - 1.0);
	}

	// UNDX�̕����̎��ɓ���
	// Rosenbrock�֐� Star�^
	/*
	for (i = 2; i<d; i++) {
	sum += (100 * (x[0] - x[i] * x[i])*(x[0] - x[i] * x[i]) + (x[1] - 1.0)*(x[1] - 1.0));
	}
	*/
	return sum - 3.98662385;
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
	printf("%lf\n", sum);
	return (418.9828872724338 * d - sum);
}
//------------------------------------------------------------
// F7 Ridge�֐�
//------------------------------------------------------------
double Ridge(double *x) {
	int i, j;
	double sum1, sum2;
	for (i = 0, sum1 = 0.0; i<d; i++) {
		for (j = 0, sum2 = 0; j <= i; j++) {
			sum2 += x[j];
		}
		sum1 += sum2 * sum2;
	}
	return sum1;
}
//------------------------------------------------------------
// F8 Bohachevsky�֐�
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
// F9 Schaffer�֐�
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
// F10 Ellipsoid�֐�
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
// F11 k-tablet�֐�
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
// F12 Shifted-Rastrigin�֐�
//------------------------------------------------------------
double Shifted_Rastrigin(double *x) {
	int i;
	double sum1 = 0.0, sum2 = 0.0;
	for (i = 0; i<d; i++) sum1 += (x[i] - 1) * (x[i] - 1) - 10 * cos(2 * M_PI*(x[i] - 1));
	sum2 = 10 * d + sum1;
	return sum2;
}
//------------------------------------------------------------
// F13 Cigar�֐�
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
// F13 2n-minima�֐�
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
	else if (Func_No == 7)return  Ridge(x);
	else if (Func_No == 8)return  Bohachevsky(x);
	else if (Func_No == 9)return  Schaffer(x);
	else if (Func_No == 10)return Ellipsoid(x);
	else if (Func_No == 11)return K_Tablet(x);
	else if (Func_No == 12)return Shifted_Rastrigin(x);
	else if (Func_No == 13)return Cigar(x);
	else if (Func_No == 14)return N2_Minima(x);
	else {
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


void Select_pVector1(int pop1)
{
	register int i;
	int pop2, pop3;
	do {
		pop2 = (int)Np*genrand_real1();
		pop3 = (int)Np*genrand_real1();
	} while (pop1 == pop2 || pop1 == pop3 || pop2 == pop3);
	for (i = 0; i<d; i++)pVect1[i] = cVect[pop2][i];
	for (i = 0; i<d; i++)pVect2[i] = cVect[pop3][i];

}
//best�̔z��̃f�[�^�̍X�V
void best_Initialize() {
	int i, j;
	for (i = 0; i < MaxGrnration; i++) {
		for (j = 0; j < EXTIME; j++) {
			Uf_best_History[j][i] = 0;
			Uf_rand_History[j][i] = 0;
			gBestHistory[j][i] = 0;
		}
	}


}


//�p�����[�^�̍X�V	JADE �쐬��
void New_parameter() {
	int i;
	double Sn = 0.0, Sf = 0.0, Sf2 = 0.0, Scr = 0.0;
	double hoge_fp = 0;
	for (i = 0; i < Np; i++) {
		if (nFitness[i] < cFitness[i]) {
			//		if (nFitness[i] < gBestFitness) {
			Sn += 1;
			Sf += F[i];
			Sf2 += F[i] * F[i];
			Scr += CR[i];
			//			printf("F=%lf\n",F[i]);
			//			if (cFitness[i] < gBestFitness) hoge_fp += 1;
			//			printf("nfit%lf,gBest%f\n", cFitness[i], gBestFitness);
		}
	}


	if (Sn != 0) {

		Ucr = (1 - C)*Ucr + C * ((Scr / Sn));
		Uf = (1 - C)*Uf + C * (Sf2 / Sf);



	}
	else {
		//		Uf += 0.05;
		//		Ucr += -0.01;
	}
}


//�p�����[�^�̍X�V	JADE �쐬��
void New_parameter_2() {
	int i;
	double Sn = 0.0, Scr = 0.0, Sn_best = 0;
	double Sf_rand = 0, Sf_rand_2 = 0;
	double Sf_best = 0, Sf_best_2 = 0;
	int hoge_fp = 0;
	for (i = 0; i < Np; i++) {
		if (cBestFitness > cFitness[i]) {
			Sf_best += F_best[i];
			Sf_best_2 += F_best[i] * F_best[i];
			hoge_fp += 1;
			Sn_best += 1;
		}
		if (nFitness[i] > cFitness[i]) {
			Sn += 1;
			Sf_rand += F_rand[i];
			Sf_rand_2 += F_rand[i] * F_rand[i];
			Scr += CR[i];
		}
	}

	if (Sn > 5) {
		Ucr = (1 - C)*Ucr + C * ((Scr / Sn));
		Uf_rand = (1 - C)*Uf_rand + C * (Sf_rand_2 / Sf_rand);
	}


	if (hoge_fp > 0) {
		Uf_best = (1 - C)*Uf_best + C * (Sf_best_2 / Sf_best);
		//		Uf_best = (1 - C)*Uf_best + C * (Sf_best/Sn_best);
	}
	else {
//		printf("UF_best=%f\n", Uf_best);
		Uf_best -= 0.01;
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

void New_parameter3() {
	int i;
	for (i = 0; i < Np; i++) {
		if (nFitness[i] < cFitness[i]) {
			if (genrand_real2() < Tf) {
				F[i] = Fl + genrand_real2()*Fu;
			}
			if (genrand_real2() < Tcr) {
				CR[i] = genrand_real2();
			}

		}
	}
}




//�p�����[�^�̍X�V
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
}

void Parameter_Format_2() {
	int i = 0;
	double sum_F = 0, sum_CR = 0;

	//	printf("F=%lf\n",F[i]);
	//	printf("CR=%lf\n", CR[i]);
	//	printf("Uf_best=%lf\n", Uf_best);
	for (i = 0; i < Np; i++) {
		do {
			F_best[i] = rand_cauchy(Uf_best, 0.5);
		} while (F_best[i] < 0.1);
		if (F_best[i] > 0.9) {
			F_best[i] = 0.9;
		}
		//			if (F_best[i] == 0.000) F_best[i] = 0.001;

		do {
			F_rand[i] = rand_cauchy(Uf_rand, 0.5);
		} while (F_rand[i] < 0.1);
		if (F_rand[i] > 0.9) {
			F_rand[i] = 0.9;
		}


		do {
			CR[i] = rand_cauchy(Ucr, 0.5);
		} while (CR[i] < 0.1);
		if (CR[i] > 0.9) {
			CR[i] = 0.9;
		}
	}

}

//�p�����[�^�̏�����
void Parameter_Initialization() {
	int i;
	for (i = 0; i < Np; i++) {
		F[i] = 0.5;
		CR[i] = 0.5;
		//		F_best[i] = 0.85;
		//		F_rand[i] = 0.5;
	}
	Uf = 0.5;
	Ucr = 0.5;
	Uf_rand = 0.5;
	Uf_best = 0.5;
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
	/*	for (i = 0; i < d; i++) {
	printf("�ω��OcVect[%d][%d]=%0.10f\n", arry_num1, i, cVect[arry_num1][i]);
	}
	for (i = 0; i < d; i++) {
	printf("�ω��OcVect[%d][%d]=%0.10f\n", arry_num2, i, cVect[arry_num2][i]);
	}
	printf("\n");
	*/
	for (i = 0; i < D; i++) {
		arry_data[i] = cVect[arry_num1][i];
		cVect[arry_num1][i] = cVect[arry_num2][i];
		cVect[arry_num2][i] = arry_data[i];
	}
	/*
	for (i = 0; i < d; i++) {
	printf("�ω���cVect[%d][%d]=%0.10f\n", arry_num1, i, cVect[arry_num1][i]);
	}
	for (i = 0; i < d; i++) {
	printf("�ω���cVect[%d][%d]=%0.10f\n", arry_num2, i, cVect[arry_num2][i]);
	}
	printf("\n\n");
	*/
}

/* �z�� a �� [left, right) ���\�[�g���܂� */
void QuickSort(double *a, int left, int right) {
	if (right - left <= 1) return;

	int pivot_index = (left + right) / 2;  // �K���ɂ����ł͒��_�Ƃ��܂�
	int pivot = a[pivot_index];
	swap(&a[pivot_index], &a[right - 1]);    // pivot �ƉE�[�� swap
	vnVect(pivot_index, right - 1);
	vcVect(pivot_index, right - 1);

	int i = left; // iterator
	for (int j = left; j < right - 1; ++j) { // j �͑S�̂𒭂߂�
		if (a[j] < pivot) { // pivot �����̂��������獶�ɋl�߂Ă���
			swap(&a[i++], &a[j]);
			vnVect((i + 1), j);
			vcVect((i + 1), j);
		}
	}
	swap(&a[i], &a[right - 1]); // pivot ��K�؂ȏꏊ�ɑ}��
	vnVect(i, (right - 1));
	vcVect(i, (right - 1));
	/* �ċA�I�ɉ��� */
	QuickSort(a, left, i);    // ������ (pivot ����)
	QuickSort(a, i + 1, right); // �E���� (pivot �ȏ�)
}


/* �z�� a �� [left, right) ���\�[�g���܂� */
void QuickSort() {

	int left = 0, right = Np;
	if (right - left <= 1) return;

	int pivot_index = (left + right) / 2;  // �K���ɂ����ł͒��_�Ƃ��܂�
	double pivot = nFitness[pivot_index];
	swap(&nFitness[pivot_index], &nFitness[right - 1]);    // pivot �ƉE�[�� swap
	vnVect(pivot_index, (right - 1));
	vcVect(pivot_index, (right - 1));
	int i = left; // iterator
	for (int j = left; j < right - 1; ++j) { // j �͑S�̂𒭂߂�
		if (nFitness[j] < pivot) { // pivot �����̂��������獶�ɋl�߂Ă���
			swap(&nFitness[i++], &nFitness[j]);
			vnVect(i + 1, j);
			vcVect(i + 1, j);
		}
	}
	swap(&nFitness[i], &nFitness[right - 1]); // pivot ��K�؂ȏꏊ�ɑ}��
	vnVect(i + 1, right - 1);
	vcVect(i + 1, right - 1);
	/* �ċA�I�ɉ��� */
	QuickSort(nFitness, left, i);    // ������ (pivot ����)
	QuickSort(nFitness, i + 1, right); // �E���� (pivot �ȏ�)
}


void  bubbleSort(double *N) {
	int i, j;
	double temp, c_temp;
	double F_temp, CR_temp;
	for (i = 0; i<Np; i++) {
		for (j = Np - 1; j>i; j--) {
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

void  Initialize_bubbleSort(double *N) {
	int i, j;
	double temp;
	for (i = 0; i<Np; i++) {
		for (j = Np - 1; j>i; j--) {
			if (N[j] < N[j - 1]) {
				swap(&N[j], &N[j - 1]);
				vcVect(j, j - 1);
			}
		}
	}
}






//------------------------------------------------------------
//DE�̑���
//------------------------------------------------------------
void DE_Operation(int i_Np, int g_GSIZE)
{
	register int i;
	int N = 0, L = 0;
	int best_rand_N;
	//DE/rand/1/exp
	if (DeAlgorithmNo == 1) {
		for (i = 0; i<d; i++) nVect[i_Np][i] = cVect[i_Np][i];
		N = (int)(genrand_real2()*d);
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
		N = (int)(genrand_real2()*d);
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
		N = (int)(genrand_real2()*d);
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
		N = (int)(genrand_real2()*d);
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


	//JADE	DE/rand/exp
	else if (DeAlgorithmNo == 5) {
		L = 0;
		for (i = 0; i < d; i++) nVect[i_Np][i] = cVect[i_Np][i];
//		N = (int)(genrand_real1()*d);
		N = 0;
		//		printf("best_P%d\n",best_rand_N);
		best_rand_N = (int)(genrand_real2()*p_best);
		do {
			if (genrand_real1()<CR[i_Np] || L == 0) {
				nVect[i_Np][N] = cVect[i_Np][N] + (F[i_Np])* (cVect[best_rand_N][N] - cVect[i_Np][N]) + (F[i_Np])* (pVect1[N] - pVect2[N]);
				//nVect[i_Np][N] = cVect[i_Np][N] + (F[i_Np] * (gBestVector[N] - cVect[i_Np][N]) + F[i_Np] * (pVect1[N] - pVect2[N]))/2;
				if (nVect[N][N] < -Range) 	nVect[i_Np][N] = pVect1[N] + genrand_real1() * (-Range - pVect1[N]);
				if (nVect[i_Np][N] > Range) nVect[i_Np][N] = pVect1[N] + genrand_real1() * (Range - pVect1[N]);
			}
			else {
				nVect[i_Np][N] = cVect[i_Np][N];
				//nVect[i_Np][N] = nVect[i_Np][N] = pVect1[N] + F[i_Np] * (pVect2[N] - pVect3[N]);
			}

//			N = (N + 1) % d;
			L++;
			N = L;
		} while (L < d);

	}
	//JADE/exp
	else if (DeAlgorithmNo == 8) {

		L = 0;
		for (i = 0; i < d; i++) nVect[i_Np][i] = cVect[i_Np][i];
		N = (int)(genrand_real1()*d);
		do {
			if (genrand_real1() < 0.2 || L == 0) {
				//				nVect[i_Np][N] = pVect1[i_Np] + (F_best[i_Np]) * (gBestVector[N] - cVect[i_Np][N]) + (1-F_rand[i_Np]) * (pVect1[i_Np] - pVect2[i_Np]);	//�ύX��
				nVect[i_Np][N] = cVect[i_Np][N] + 0.4*(gBestVector[N] - pVect1[N]) + 0.0* (pVect2[N] - pVect3[N]);
				if (nVect[i_Np][N] < -Range) 	nVect[i_Np][N] = pVect1[i_Np] + genrand_real1() * (-Range - pVect1[i_Np]);
				if (nVect[i_Np][N] > Range)		nVect[i_Np][N] = pVect1[i_Np] + genrand_real1() * (Range - pVect1[i_Np]);
			}
			else {
				//				nVect[i_Np][N] = cVect[i_Np][N];
			}
			N = (N + 1) % d;
			L++;

		} while (L < d);

	}
	else if (DeAlgorithmNo == 7) {

		L = 0;
		for (i = 0; i < d; i++) nVect[i_Np][i] = cVect[i_Np][i];
		N = (int)(genrand_real1()*d);
		do {
			if (genrand_real1() < CR[i_Np] || L == 0) {
				//				nVect[i_Np][N] = cVect[i_Np][N] + (F_best[i_Np]) * (gBestVector[N] - pVect1[N]) + (F_rand[N]) * (pVect3[N] - pVect2[N]);
				nVect[i_Np][N] = pVect1[N] + (F_best[i_Np]) * (gBestVector[N] - pVect1[N]) + (F_rand[N]) * (pVect3[N] - pVect2[N]);
				if (nVect[i_Np][N] < -Range) 	nVect[i_Np][N] = pVect1[i_Np] + genrand_real1() * (-Range - pVect1[N]);
				if (nVect[i_Np][N] > Range)	nVect[i_Np][N] = pVect1[i_Np] + genrand_real1() * (Range - pVect1[N]);
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
		/*
		L = 0;
		for (i = 0; i < d; i++) nVect[i_Np][i] = cVect[i_Np][i];
		N = (int)(genrand_real1()*d);
		do {
		if (genrand_real1() < CR[i_Np] || L == 0) {
		nVect[i_Np][N] = pVect1[N] + F[i_Np] * (pVect2[N] - pVect3[N]);
		if (nVect[i_Np][N] < -Range) 	nVect[i_Np][N] = pVect1[i_Np] + genrand_real1() * (-Range - pVect1[N]);
		if (nVect[i_Np][N] > Range)	nVect[i_Np][N] = pVect1[i_Np] + genrand_real1() * (Range - pVect1[N]);
		}
		else {
		nVect[i_Np][N] = cVect[i_Np][N];
		}
		N = (N + 1) % d;
		L++;
		} while (L < d);
		*/
		L = 0;
		for (i = 0; i < d; i++) nVect[i_Np][i] = cVect[i_Np][i];
		N = (int)(genrand_real1()*d);
		//printf("best_P%d\n",best_rand_N);
		do {
			if (genrand_real1()<CR[i_Np] || L == 0) {
				best_rand_N = (int)(genrand_real2()*p_best);
				nVect[i_Np][N] = cVect[i_Np][N] + (F[i_Np])* (cVect[best_rand_N][N] - cVect[i_Np][N] + (F[i_Np])* (pVect1[N] - pVect2[N]));
				//nVect[i_Np][N] = cVect[i_Np][N] + (F[i_Np] * (gBestVector[N] - cVect[i_Np][N]) + F[i_Np] * (pVect1[N] - pVect2[N]))/2;
				if (nVect[N][N] < -Range) 	nVect[i_Np][N] = pVect1[N] + genrand_real1() * (-Range - pVect1[N]);
				if (nVect[i_Np][N] > Range) nVect[i_Np][N] = pVect1[N] + genrand_real1() * (Range - pVect1[N]);
			}
			else {
				nVect[i_Np][N] = cVect[i_Np][N];
				//nVect[i_Np][N] = nVect[i_Np][N] = pVect1[N] + F[i_Np] * (pVect2[N] - pVect3[N]);
			}

			N = (N + 1) % d;
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
			for (j = 0; j<d; j++)	cVect[i][j] = nVect[i][j];
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
	//	printf("%20.10lf\n", gBestFitness);
}


//������
void Initialize() {
	episode = 0;
	No_best_sum = 0;
	Init_Vector();
	Evaluate_Init_Vector();
	Parameter_Initialization();
	if(DeAlgorithmNo==5) Initialize_bubbleSort(cFitness);
}

void vJade_Parameter_Format() {
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


}

void vJade_Parameter_Format2() {
	int i;
	double Sn = 0.0, Scr = 0.0, Sn_best = 0;
	double Sf_rand = 0, Sf_rand_2 = 0;
	double Sf_best = 0, Sf_best_2 = 0;
	int hoge_fp = 0;
	for (i = 0; i < Np; i++) {
		if (cBestFitness > cFitness[i]) {
			Sf_best += F_best[i];
			Sf_best_2 += F_best[i] * F_best[i];
			hoge_fp += 1;
			Sn_best += 1;
		}
		if (nFitness[i] > cFitness[i]) {
			Sn += 1;
			Sf_rand += F_rand[i];
			Sf_rand_2 += F_rand[i] * F_rand[i];
			Scr += CR[i];
		}
	}

	if (Sn > 5) {
		Ucr = (1 - C)*Ucr + C * ((Scr / Sn));
		Uf_rand = (1 - C)*Uf_rand + C * (Sf_rand_2 / Sf_rand);
	}


	if (hoge_fp > 0) {
		Uf_best = (1 - C)*Uf_best + C * (Sf_best_2 / Sf_best);
		//		Uf_best = (1 - C)*Uf_best + C * (Sf_best/Sn_best);
	}
	else {
//		printf("UF_best=%f\n", Uf_best);
		Uf_best -= 0.01;
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

void vJade_Parameter_Format3() {
	int i;
	double Sn = 0.0, Scr = 0.0, Sn_best = 0;
	double Sf_rand = 0, Sf_rand_2 = 0;
	double Sf_best = 0, Sf_best_2 = 0;
	int hoge_fp = 0;
	for (i = 0; i < Np; i++) {
		if (cBestFitness > cFitness[i]) {
			Sf_best += F_best[i];
			Sf_best_2 += F_best[i] * F_best[i];
			hoge_fp += 1;
			Sn_best += 1;
		}
		if (nFitness[i] > cFitness[i]) {
			Sn += 1;
			Sf_rand += F_rand[i];
			Sf_rand_2 += F_rand[i] * F_rand[i];
			Scr += CR[i];
		}
	}

	if (Sn > 5) {
		Ucr = (1 - C)*Ucr + C * ((Scr / Sn));
		Uf_rand = (1 - C)*Uf_rand + C * (Sf_rand_2 / Sf_rand);
	}


	if (hoge_fp > 0) {
		Uf_best = (1 - C)*Uf_best + C * (Sf_best_2 / Sf_best);
		//		Uf_best = (1 - C)*Uf_best + C * (Sf_best/Sn_best);
	}
	else {
//		printf("UF_best=%f\n", Uf_best);
		Uf_best -= 0.01;
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

void vjde_Parameter_Format() {
	int i;
	for (i = 0; i < Np; i++) {
		if (nFitness[i] < cFitness[i]) {
			if (genrand_real2() < Tf) {
				F[i] = Fl + genrand_real2()*Fu;
			}
			if (genrand_real2() < Tcr) {
				CR[i] = genrand_real2();
			}

		}
	}
}

//�p�����[�^�̍X�V
void Parameter_Format(int iDe_nomber) {

	if (iDe_nomber == 5) {
		New_parameter();
		vJade_Parameter_Format();
	}
	else if (iDe_nomber == 6 || iDe_nomber == 7) {
		New_parameter_2();
		vJade_Parameter_Format2();
	}
	else if (iDe_nomber == 8) {
		vjde_Parameter_Format();

	}


}





