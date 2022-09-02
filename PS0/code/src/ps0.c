# include <stdio.h>
# include <stdlib.h>
// Task 1

int N = 200;

float * allocateArray(int n) {
    return malloc(n * sizeof(float));
}

float square(float x) {
    return x * x;
}

int main() {
    float* function_vals = allocateArray(N);
    float* derivative_vals = allocateArray(N);
    return 1;
}