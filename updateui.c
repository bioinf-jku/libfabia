
#include <stdio.h>

// very basic updateUI, used only when compiling libfabia.so as standalone library
void updateUI(int iter, float elapsedTime, int k, int n, int l, float* L, float* Z, float* Psi, float* lapla){
	printf("Elapsed time: %7.1fs Cycle: %d\n", elapsedTime, iter);
}
