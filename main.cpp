#include <iostream>
#include <vector>
#include <ctime>
#include <algorithm>
#include <fstream>
#include <climits>

using namespace std;


typedef vector<int> Row;
typedef vector<Row> Matrix;


void PrintMatrix(const Matrix &, int);
void MultiplyNaive(const Matrix &, const Matrix &, Matrix &);
void MultiplyNaive2(const Matrix &, const Matrix &, Matrix &, int, int, int, int, size_t);

void MultiplyRecursive(const Matrix &, const Matrix &, Matrix &, int, int, int, int, size_t, int);
void MultiplyRecursive2(const Matrix &, const Matrix &, Matrix &, int, int, int, int, size_t);


void SumMatrix(const Matrix &, const Matrix &, Matrix &, size_t);
void SumMatrix2(const Matrix &, const Matrix &, Matrix &, int, int, int, int, size_t);
void SubMatrix(const Matrix &A, const Matrix &B, Matrix &, size_t);
void SubMatrix2(const Matrix &A, const Matrix &B, Matrix &, int, int, int, int, size_t);

void MultiplyStrassen(const Matrix &, const Matrix &, Matrix &, size_t);
void MultiplyStrassen2(const Matrix &, const Matrix &, Matrix &, size_t);
void MultiplyStrassen3(const Matrix &, const Matrix &, Matrix &, int, int, int, int, size_t);

int main() {

    size_t num = 256;
    Matrix A(num, Row(num));
    Matrix B(num, Row(num));
    Matrix C(num, Row(num,0));
    Matrix D(num, Row(num,0));
    Matrix E(num, Row(num,0));
    Matrix F(num, Row(num,0));
    srand(time(0));
    for (int i = 0; i < num; ++i) {
        for (int j = 0; j < num; ++j) {
            //A[i][j] = rand() % (int)(INT_MAX / sqrt(num));
            //B[i][j] = rand() % (int)(INT_MAX / sqrt(num));
            A[i][j] = rand() % 10;
            B[i][j] = rand() % 10;
        }
    }

    /*
    clock_t begin = clock();
    MultiplyNaive(A, B, D);

    clock_t end = clock();
    double diff = end - begin;

    cout<<"Naive: "<< diff/CLOCKS_PER_SEC<<endl;
     */
    /*
    ofstream myfile;
    myfile.open ("test2.txt");

    clock_t begin, end;
    double diff;
    double rec1=0, rec2=0, rec3=0;
    int run=4;
    int cutoff = 2;
    for (int l = 0; l < 7; ++l) {
        rec1 = 0;
        cout<<l<<endl;
        for (int k = 0; k < run; ++k) {
            begin = clock();
            MultiplyRecursive(A, B, C, 0, 0, 0, 0, num, cutoff);
            end = clock();
            diff = end - begin;
            rec1 += diff / CLOCKS_PER_SEC;
            if (cutoff == 2) {
                begin = clock();
                MultiplyRecursive2(A, B, D, 0, 0, 0, 0, num);
                end = clock();
                diff = end - begin;
                rec2 += diff / CLOCKS_PER_SEC;

                begin = clock();
                MultiplyStrassen(A, B, D, num);
                end = clock();
                diff = end - begin;
                rec3 += diff / CLOCKS_PER_SEC;
            }
        }
        if (cutoff == 2) {
            myfile << 0 << " "<< rec2/run <<"\n";
            myfile << 1 << " "<< rec3/run <<"\n";
        }
        myfile << cutoff << " "<< rec1/run <<"\n";
        cutoff *= 2;
    }

    myfile.close();
    cout<<"Recursive: "<< rec1/run <<endl;
    cout<<"Recursive 2: "<< rec2/run <<endl;
    */




    clock_t begin = clock();
    MultiplyRecursive2(A, B, D, 0,0,0,0,num);
    clock_t end = clock();
    double diff = end - begin;
    cout<<"Recursive2: "<< diff/CLOCKS_PER_SEC<<endl;

    begin = clock();
    MultiplyNaive(A, B, C);
    end = clock();
    diff = end - begin;
    cout<<"Naive: "<< diff/CLOCKS_PER_SEC<<endl;


    begin = clock();
    MultiplyStrassen3(A, B, F, 0,0,0,0, num);
    end = clock();
    diff = end - begin;
    cout<<"Strassen3: "<< diff/CLOCKS_PER_SEC<<endl;





    bool igual = (C==D);
    cout<<igual<<endl;

    igual = (C==F);
    cout<<igual<<endl;

    //PrintMatrix(C, num);
    //PrintMatrix(D, num);
    //PrintMatrix(E, num);
    return  0;
}


void PrintMatrix(const Matrix &matrix, int length) {
    cout<<endl;
    for (int i = 0; i < length; i++) {
        for (int j = 0; j < length; j++)
            cout << matrix[i][j] << ' ';
        cout << endl;
    }
    cout<<endl;
}


void MultiplyNaive(const Matrix &A, const Matrix &B, Matrix &C){
    size_t n = A.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void MultiplyNaive2(const Matrix &A, const Matrix &B, Matrix &C, int row_a, int column_a, int row_b, int column_b, size_t n){
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                C[i][j] += A[i+row_a][k+column_a] * B[k+row_b][j+column_b];
            }
        }
    }
}

void MultiplyRecursive(const Matrix &A, const Matrix &B, Matrix &C, int row_a, int column_a, int row_b, int column_b, size_t length, int cutoff) {
    if (length == cutoff) {
        MultiplyNaive2(A, B, C, row_a, column_a, row_b, column_b, length);
    }
    else {
        size_t new_length = length / 2;
        Matrix R(new_length, Row(new_length,0));
        Matrix S(new_length, Row(new_length,0));
        Matrix T(new_length, Row(new_length,0));
        Matrix U(new_length, Row(new_length,0));
        Matrix Tmp1(new_length, Row(new_length,0));
        Matrix Tmp2(new_length, Row(new_length,0));

        MultiplyRecursive(A, B, Tmp1, row_a, column_a, row_b, column_b, new_length, cutoff);
        MultiplyRecursive(A, B, Tmp2, row_a, column_a + new_length, row_b + new_length, column_b, new_length, cutoff);
        SumMatrix(Tmp1, Tmp2, R, new_length);

        MultiplyRecursive(A, B, Tmp1, row_a, column_a, row_b, column_b + new_length, new_length, cutoff);
        MultiplyRecursive(A, B, Tmp2, row_a, column_a + new_length, row_b + new_length, column_b + new_length, new_length, cutoff);
        SumMatrix(Tmp1, Tmp2, S, new_length);

        MultiplyRecursive(A, B, Tmp1, row_a + new_length, column_a, row_b, column_b, new_length, cutoff);
        MultiplyRecursive(A, B, Tmp2, row_a + new_length, column_a + new_length, row_b + new_length, column_b, new_length, cutoff);
        SumMatrix(Tmp1, Tmp2, T, new_length);

        MultiplyRecursive(A, B, Tmp1, row_a + new_length, column_a, row_b, column_b + new_length, new_length, cutoff);
        MultiplyRecursive(A, B, Tmp2, row_a + new_length, column_a + new_length, row_b + new_length, column_b + new_length, new_length, cutoff);
        SumMatrix(Tmp1, Tmp2, U, new_length);

        for (size_t i = 0; i < new_length ; i++) {
            for (size_t j = 0 ; j < new_length ; j++) {
                C[i][j] = R[i][j];
                C[i][j + new_length] = S[i][j];
                C[i + new_length][j] = T[i][j];
                C[i + new_length][j + new_length] = U[i][j];
            }
        }

    }
}

void MultiplyRecursive2(const Matrix &A, const Matrix &B, Matrix &C, int row_a, int column_a, int row_b, int column_b, size_t length) {
    if (length == 1) {
        C[0][0] = A[row_a][column_a]*B[row_b][column_b];
    }
    else {
        size_t new_length = length / 2;
        Matrix R(new_length, Row(new_length,0));
        Matrix S(new_length, Row(new_length,0));
        Matrix T(new_length, Row(new_length,0));
        Matrix U(new_length, Row(new_length,0));
        Matrix Tmp1(new_length, Row(new_length,0));
        Matrix Tmp2(new_length, Row(new_length,0));

        MultiplyRecursive2(A, B, Tmp1, row_a, column_a, row_b, column_b, new_length);
        MultiplyRecursive2(A, B, Tmp2, row_a, column_a + new_length, row_b + new_length, column_b, new_length);
        SumMatrix(Tmp1, Tmp2, R, new_length);

        MultiplyRecursive2(A, B, Tmp1, row_a, column_a, row_b, column_b + new_length, new_length);
        MultiplyRecursive2(A, B, Tmp2, row_a, column_a + new_length, row_b + new_length, column_b + new_length, new_length);
        SumMatrix(Tmp1, Tmp2, S, new_length);

        MultiplyRecursive2(A, B, Tmp1, row_a + new_length, column_a, row_b, column_b, new_length);
        MultiplyRecursive2(A, B, Tmp2, row_a + new_length, column_a + new_length, row_b + new_length, column_b, new_length);
        SumMatrix(Tmp1, Tmp2, T, new_length);

        MultiplyRecursive2(A, B, Tmp1, row_a + new_length, column_a, row_b, column_b + new_length, new_length);
        MultiplyRecursive2(A, B, Tmp2, row_a + new_length, column_a + new_length, row_b + new_length, column_b + new_length, new_length);
        SumMatrix(Tmp1, Tmp2, U, new_length);

        for (size_t i = 0; i < new_length ; i++) {
            for (size_t j = 0 ; j < new_length ; j++) {
                C[i][j] = R[i][j];
                C[i][j + new_length] = S[i][j];
                C[i + new_length][j] = T[i][j];
                C[i + new_length][j + new_length] = U[i][j];
            }
        }

    }
}


void MultiplyStrassen(const Matrix &A, const Matrix &B, Matrix &C, size_t length) {
    if (length == 8) {
        return MultiplyNaive(A, B, C);
    }
    else {
        int i, j;
        size_t new_length = (size_t ) length / 2;
        Matrix Ar(new_length, Row(new_length)), As(new_length, Row(new_length)),
                At(new_length, Row(new_length)), Au(new_length, Row(new_length)),
                Br(new_length, Row(new_length)), Bs(new_length, Row(new_length)),
                Bt(new_length, Row(new_length)), Bu(new_length, Row(new_length)),
                Cr(new_length, Row(new_length)), Cs(new_length, Row(new_length)),
                Ct(new_length, Row(new_length)), Cu(new_length, Row(new_length)),
                P1(new_length, Row(new_length)), P2(new_length, Row(new_length)),
                P3(new_length, Row(new_length)), P4(new_length, Row(new_length)),
                P5(new_length, Row(new_length)), P6(new_length, Row(new_length)),
                P7(new_length, Row(new_length)), Tmp1(new_length, Row(new_length)),
                Tmp2(new_length, Row(new_length));


        for (i = 0; i < new_length; i++) {
            for (j = 0; j < new_length; j++) {
                Ar[i][j] = A[i][j];
                As[i][j] = A[i][j + new_length];
                At[i][j] = A[i + new_length][j];
                Au[i][j] = A[i + new_length][j + new_length];

                Br[i][j] = B[i][j];
                Bs[i][j] = B[i][j + new_length];
                Bt[i][j] = B[i + new_length][j];
                Bu[i][j] = B[i + new_length][j + new_length];
            }
        }
        SumMatrix(Ar, Au, Tmp1, new_length); // a11 + a22
        SumMatrix(Br, Bu, Tmp2, new_length); // b11 + b22
        MultiplyStrassen(Tmp1, Tmp2, P1, new_length); // p1 = (a11+a22) * (b11+b22)

        SumMatrix(At, Au, Tmp1, new_length); // a21 + a22
        MultiplyStrassen(Tmp1, Br, P2, new_length); // p2 = (a21+a22) * (b11)

        SubMatrix(Bs, Bu, Tmp2, new_length); // b12 - b22
        MultiplyStrassen(Ar, Tmp2, P3, new_length); // p3 = (Ar) * (b12 - b22)

        SubMatrix(Bt, Br, Tmp2, new_length); // b21 - b11
        MultiplyStrassen(Au, Tmp2, P4, new_length); // p4 = (a22) * (b21 - b11)

        SumMatrix(Ar, As, Tmp1, new_length); // a11 + a12
        MultiplyStrassen(Tmp1, Bu, P5, new_length); // p5 = (a11+a12) * (b22)

        SubMatrix(At, Ar, Tmp1, new_length); // a21 - a11
        SumMatrix(Br, Bs, Tmp2, new_length); // b11 + b12
        MultiplyStrassen(Tmp1, Tmp2, P6, new_length); // p6 = (a21-a11) * (b11+b12)

        SubMatrix(As, Au, Tmp1, new_length); // a12 - a22
        SumMatrix(Bt, Bu, Tmp2, new_length); // b21 + b22
        MultiplyStrassen(Tmp1, Tmp2, P7, new_length); // p7 = (a12-a22) * (b21+b22)

        // calculating c21, c21, c11 e c22:

        SumMatrix(P3, P5, Cs, new_length); // c12 = p3 + p5
        SumMatrix(P2, P4, Ct, new_length); // c21 = p2 + p4

        SumMatrix(P1, P4, Tmp1, new_length); // p1 + p4
        SumMatrix(Tmp1, P7, Tmp2, new_length); // p1 + p4 + p7
        SubMatrix(Tmp2, P5, Cr, new_length); // c11 = p1 + p4 - p5 + p7

        SumMatrix(P1, P3, Tmp1, new_length); // p1 + p3
        SumMatrix(Tmp1, P6, Tmp2, new_length); // p1 + p3 + p6
        SubMatrix(Tmp2, P2, Cu, new_length); // c22 = p1 + p3 - p2 + p6
        
        for (i = 0; i < new_length ; i++) {
            for (j = 0 ; j < new_length ; j++) {
                C[i][j] = Cr[i][j];
                C[i][j + new_length] = Cs[i][j];
                C[i + new_length][j] = Ct[i][j];
                C[i + new_length][j + new_length] = Cu[i][j];
            }
        }
    }
}


void MultiplyStrassen2(const Matrix &A, const Matrix &B, Matrix &C, size_t length) {
    if (length == 1) {
        C[0][0] = A[0][0]*B[0][0];
    }
    else {
        int i, j;
        size_t new_length = (size_t ) length / 2;
        Matrix A11(new_length, Row(new_length)), A12(new_length, Row(new_length)),
                A21(new_length, Row(new_length)), A22(new_length, Row(new_length)),
                B11(new_length, Row(new_length)), B12(new_length, Row(new_length)),
                B21(new_length, Row(new_length)), B22(new_length, Row(new_length)),
                P1(new_length, Row(new_length)), P2(new_length, Row(new_length)),
                P3(new_length, Row(new_length)), P4(new_length, Row(new_length)),
                P5(new_length, Row(new_length)), P6(new_length, Row(new_length)),
                P7(new_length, Row(new_length)),
                S1(new_length, Row(new_length)), S2(new_length, Row(new_length)),
                S3(new_length, Row(new_length)), S4(new_length, Row(new_length)),
                T1(new_length, Row(new_length)), T2(new_length, Row(new_length)),
                T3(new_length, Row(new_length)), T4(new_length, Row(new_length)),
                U1(new_length, Row(new_length)), U2(new_length, Row(new_length)),
                U3(new_length, Row(new_length)), U4(new_length, Row(new_length)),
                U5(new_length, Row(new_length)), U6(new_length, Row(new_length)),
                U7(new_length, Row(new_length));


        for (i = 0; i < new_length; i++) {
            for (j = 0; j < new_length; j++) {
                A11[i][j] = A[i][j];
                A12[i][j] = A[i][j + new_length];
                A21[i][j] = A[i + new_length][j];
                A22[i][j] = A[i + new_length][j + new_length];

                B11[i][j] = B[i][j];
                B12[i][j] = B[i][j + new_length];
                B21[i][j] = B[i + new_length][j];
                B22[i][j] = B[i + new_length][j + new_length];
            }
        }
        SumMatrix(A21, A22, S1, new_length); // S1 = A21 + A22
        SubMatrix(S1, A11, S2, new_length); // S2 = S1 - A11
        SubMatrix(A11, A21, S3, new_length); // S3 = A11 - A21
        SubMatrix(A12, S2, S4, new_length); // S4 = A12 - S2
        SubMatrix(B12, B11, T1, new_length); // T1 = B12 - B11
        SubMatrix(B22, T1, T2, new_length); // T2 = B22 - T1
        SubMatrix(B22, B12, T3, new_length); // T3 = B22 - B12
        SubMatrix(T2, B21, T4, new_length); // T4 = T2 - B21


        MultiplyStrassen2(A11, B11, P1, new_length); // P1 = A11*B11
        MultiplyStrassen2(A12, B21, P2, new_length); // P2 = A12*B21
        MultiplyStrassen2(S4, B22, P3, new_length); // P3 = S4*B22
        MultiplyStrassen2(A22, T4, P4, new_length); // P4 = A22*T4
        MultiplyStrassen2(S1, T1, P5, new_length); // P5 = S1*T1
        MultiplyStrassen2(S2, T2, P6, new_length); // P6 = S2*T2
        MultiplyStrassen2(S3, T3, P7, new_length); // P7 = S3*T3

        SumMatrix(P1, P2, U1, new_length); // U1 = P1 + P2
        SumMatrix(P1, P6, U2, new_length); // U2 = P1 + P6
        SumMatrix(U2, P7, U3, new_length); // U3 = U2 + P7
        SumMatrix(U2, P5, U4, new_length); // U4 = U2 + P5
        SumMatrix(U4, P3, U5, new_length); // U5 = U4 + P3
        SubMatrix(U3, P4, U6, new_length); // U6 = U3 - P4
        SumMatrix(U3, P5, U7, new_length); // U7 = U3 + P5


        for (i = 0; i < new_length ; i++) {
            for (j = 0 ; j < new_length ; j++) {
                C[i][j] = U1[i][j];
                C[i][j + new_length] = U5[i][j];
                C[i + new_length][j] = U6[i][j];
                C[i + new_length][j + new_length] = U7[i][j];
            }
        }
    }
}

void MultiplyStrassen3(const Matrix &A, const Matrix &B, Matrix &C, int row_a, int column_a, int row_b, int column_b, size_t length) {
    if (length == 16) {
        MultiplyRecursive2(A, B, C, row_a,column_a,row_b,column_b,length);
    }
    else {
        int i, j;
        size_t new_length = (size_t ) length / 2;
        Matrix  P1(new_length, Row(new_length)), P2(new_length, Row(new_length)),
                P3(new_length, Row(new_length)), P4(new_length, Row(new_length)),
                P5(new_length, Row(new_length)), P6(new_length, Row(new_length)),
                P7(new_length, Row(new_length)),
                S1(new_length, Row(new_length)), S2(new_length, Row(new_length)),
                S3(new_length, Row(new_length)), S4(new_length, Row(new_length)),
                T1(new_length, Row(new_length)), T2(new_length, Row(new_length)),
                T3(new_length, Row(new_length)), T4(new_length, Row(new_length)),
                U1(new_length, Row(new_length)), U2(new_length, Row(new_length)),
                U3(new_length, Row(new_length)), U4(new_length, Row(new_length)),
                U5(new_length, Row(new_length)), U6(new_length, Row(new_length)),
                U7(new_length, Row(new_length));


        SumMatrix2(A, A, S1, row_a+new_length, column_a, row_a+new_length, column_a+new_length, new_length); // S1 = A21 + A22
        SubMatrix2(S1, A, S2, 0, 0, row_a, column_a, new_length); // S2 = S1 - A11
        SubMatrix2(A, A, S3, row_a, column_a, row_a+new_length, column_a, new_length); // S3 = A11 - A21
        SubMatrix2(A, S2, S4, row_a, column_a+new_length, 0, 0, new_length); // S4 = A12 - S2
        SubMatrix2(B, B, T1, row_b, column_b+new_length, row_b, column_b, new_length); // T1 = B12 - B11
        SubMatrix2(B, T1, T2, row_b+new_length, column_b+new_length, 0, 0, new_length); // T2 = B22 - T1
        SubMatrix2(B, B, T3, row_b+new_length, column_b+new_length, row_b, column_b+new_length, new_length); // T3 = B22 - B12
        SubMatrix2(T2, B, T4, 0, 0, row_b+new_length, column_b, new_length); // T4 = T2 - B21


        MultiplyStrassen3(A, B, P1, row_a, column_a, row_b, column_b, new_length); // P1 = A11*B11
        MultiplyStrassen3(A, B, P2, row_a, column_a+new_length, row_b+new_length, column_b, new_length); // P2 = A12*B21
        MultiplyStrassen3(S4, B, P3, 0, 0, row_b+new_length, column_b+new_length, new_length); // P3 = S4*B22
        MultiplyStrassen3(A, T4, P4, row_a+new_length, column_a+new_length, 0, 0, new_length); // P4 = A22*T4
        MultiplyStrassen3(S1, T1, P5, 0, 0, 0, 0, new_length); // P5 = S1*T1
        MultiplyStrassen3(S2, T2, P6, 0, 0, 0, 0, new_length); // P6 = S2*T2
        MultiplyStrassen3(S3, T3, P7, 0, 0, 0, 0, new_length); // P7 = S3*T3

        SumMatrix2(P1, P2, U1, 0, 0, 0, 0, new_length); // U1 = P1 + P2
        SumMatrix2(P1, P6, U2, 0, 0, 0, 0, new_length); // U2 = P1 + P6
        SumMatrix2(U2, P7, U3, 0, 0, 0, 0, new_length); // U3 = U2 + P7
        SumMatrix2(U2, P5, U4, 0, 0, 0, 0, new_length); // U4 = U2 + P5
        SumMatrix2(U4, P3, U5, 0, 0, 0, 0, new_length); // U5 = U4 + P3
        SubMatrix2(U3, P4, U6, 0, 0, 0, 0, new_length); // U6 = U3 - P4
        SumMatrix2(U3, P5, U7, 0, 0, 0, 0, new_length); // U7 = U3 + P5


        for (i = 0; i < new_length ; i++) {
            for (j = 0 ; j < new_length ; j++) {
                C[i][j] = U1[i][j];
                C[i][j + new_length] = U5[i][j];
                C[i + new_length][j] = U6[i][j];
                C[i + new_length][j + new_length] = U7[i][j];
            }
        }
    }
}

void SumMatrix(const Matrix &A, const Matrix &B, Matrix &C, size_t n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            C[i][j] += A[i][j] + B[i][j];
        }
    }
}

void SumMatrix2(const Matrix &A, const Matrix &B, Matrix &C, int row_a, int column_a, int row_b, int column_b, size_t n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            C[i][j] += A[i+row_a][j+column_a] + B[i+row_b][j+column_b];
        }
    }
}

void SubMatrix(const Matrix &A, const Matrix &B, Matrix &C, size_t n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            C[i][j] += A[i][j] - B[i][j];
        }
    }
}

void SubMatrix2(const Matrix &A, const Matrix &B, Matrix &C, int row_a, int column_a, int row_b, int column_b, size_t n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            C[i][j] += A[i+row_a][j+column_a] - B[i+row_b][j+column_b];
        }
    }
}






