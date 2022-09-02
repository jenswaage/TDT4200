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

void printValues(float* array, int size) {
    for(int i = 0; i < size; i = i + 1) {
        printf("%.2f\n", array[i]);
    }
}

int main() {
    float* function_vals = allocateArray(N);
    float* derivative_vals = allocateArray(N);

    for(int i = 0; i < N; i = i + 1) {
        function_vals[i] = square(i);
    }
    printValues(function_vals, N);

    return 1;
}