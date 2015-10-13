



void multiplyMatrix(int n, int m, int k, float* A, float* B, float* C, const char* tA, const char* tB) {
	float ONE = 1.0f;
	float ZERO = 0.0f;
	int lda = *tA == 'n' ? n : k;
	int ldb = *tB == 'n' ? m : k;
	F77_CALL(gemm)(tA,tB, &n, &m, &k, &ONE, A, &lda, B, &ldb, &ZERO, C, &n);
}



void multiplySymmetricMatrix(int n, int m, float* a, float* b, float* c) {
	float ONE = 1.0f;
	float ZERO = 0.0f;
	F77_CALL(symm)("l", "l", &n,&m,&ONE,a,&n, b,&n,&ZERO,c,&n);
}

void multiplyMatrixVector(int k, int n, float* a, float* b, float* c) {
	float ONE = 1.0f;
	float ZERO = 0.0f;
	int ONE_INT = 1;
	F77_CALL(gemv)("n", &k, &n, &ONE, a, &k, b, &ONE_INT, &ZERO, c, &ONE_INT);
}


void rankOneUpdate(int n, int k, float* a, float* b, float* c) {
	float ONE = 1.0f;
	int ONE_INT = 1;
	F77_CALL(ger)(&n, &k, &ONE, a, &ONE_INT,  b, &ONE_INT, c, &n);
}
