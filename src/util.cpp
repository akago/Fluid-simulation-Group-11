#include "util.h"

void printVector(int n, double* num) {
	for (int ii = 0; ii < n; ii++) {
		printf("%.8f  ", num[ii]);
	}
	printf("\n");
}

