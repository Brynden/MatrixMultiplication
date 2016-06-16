#include <iostream>
#include <vector>

using namespace std;


typedef vector<int> Row;
typedef vector<Row> Matrix;


void printMatrix(Matrix, int, int ,int );

int main() {

    //int matrix2[4][4] = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
    int num = 4;
    Matrix my_matrix(num, Row(num));
    my_matrix ={{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
    printMatrix(my_matrix, 0, 0, 4);

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


