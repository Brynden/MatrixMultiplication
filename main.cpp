#include <iostream>
#include <vector>

using namespace std;


typedef vector<int> Row;
typedef vector<Row> Matrix;


void printMatrix(Matrix, int, int ,int );
Matrix multiply(Matrix, Matrix);

int main() {

    //int matrix2[4][4] = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
    int num = 4;
    Matrix A(num, Row(num));
    Matrix B(num, Row(num));
    Matrix C(num, Row(num));
    A ={{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
    B ={{13, 14, 15, 16},{5, 6, 7, 8},{9, 10, 11, 12},{1, 2, 3, 4}};
    //printMatrix(A, 0, 0, 4);
    C = multiply(A, B);
    printMatrix(C, 0, 0, num);
    return  0;
}


void printMatrix(Matrix matrix, int fila, int columna, int largo) {
    int n = matrix.size();
    for (int i = fila; i < fila+largo; i++) {
        for (int j = columna; j < columna+largo; j++)
            cout << matrix[i][j] << ' ';
        cout << endl;
    }
}


Matrix multiply(Matrix A, Matrix B){
    int n = A.size();
    Matrix C(n, Row(n, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
    return C;
}