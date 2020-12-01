#define infile "data.dat"
#include <iostream>
#include <fstream>
#include <string.h>
#include <memory>
#include <cmath>
#include <memory.h>
using std::cout;
using std::endl;
using std::ifstream;

void normalizeVectors(double* &matrix, int n) {
    int counter = 0;
    double* magnitudes;
    for (int i = 0; i < (sizeof(matrix) / sizeof(double)); i += n) {
        for (int j = i; j < n + i; j++, counter++) {
            magnitudes[i / n] += (matrix[counter] * matrix[counter]);
        }
        magnitudes[i / n] = sqrt(magnitudes[i / n]);
    }
    for (int i = 0; i < (sizeof(matrix) / sizeof(double)); i += n) {
        for (int j = i; j < n + i; j++, counter++) {
            matrix[counter] /= magnitudes[i / n];
        }
    }
    delete[] magnitudes;
}
double dotProduct(double* &matrix, int n, int v1, int v2) {
    double total = 0;
    for (int i = (v1*n); i < ((v1+1)*n); i++) {
        total += (matrix[i]*matrix[i+((v2-v1)*n)]);
    }
    return total;
}

void GramSchmidt(double* &matrix, int n) {
    int counter = 0;
    normalizeVectors(matrix, n);
    for (int i = 0; i < (sizeof(matrix) / sizeof(double)); i+=n) {
        for (int j = i; j < n + i; j++, counter++) {
            for (int k = 0; k < i; k++) {
                matrix[counter] -= dotProduct(matrix, n, i, k);
            }
        }
    }
    normalizeVectors(matrix, n);
}

void printMatrix(double* &matrix, int n) {
    int counter = 0;
    for (int i = 0; i < (sizeof(matrix) / sizeof(double)); i += n) {
        for (int j = i; j < n + i; j++, counter++) {
            cout << matrix[counter];
            if (j == (n + i - 1))
                cout << endl;
            else
                cout << "/t";
        }
    }
}

int main() {
    double* matrix;
    for (int i = 0; i < 9; i++) {
        matrix[i] = 2*(double)i - 1;
    }
    GramSchmidt(matrix, 3);
    printMatrix(matrix, 3);
}

