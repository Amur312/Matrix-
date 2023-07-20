#ifndef SRC_S21_MATRIX_OOP_H_
#define SRC_S21_MATRIX_OOP_H_
#include <iostream>
#include <cmath>

class S21Matrix
{
private:
    int rows_ = 0;
    int cols_ = 0;
    double **matrix_ = nullptr;

public:
    S21Matrix();
    S21Matrix(const int rows, const int cols);
    S21Matrix(const S21Matrix &other);
    S21Matrix(S21Matrix &&other);
    ~S21Matrix();

    bool EqMatrix(const S21Matrix& other) const;
    void SumMatrix(const S21Matrix &other);
    void SubMatrix(const S21Matrix &other);
    void MulNumber(const double num);
    void MulMatrix(const S21Matrix &other);

    S21Matrix Transpose();
    S21Matrix CalcComplements();
    double Determinant();
    S21Matrix InverseMatrix();

     // перегруженные операторы
    S21Matrix operator+(const S21Matrix& other) const;
    S21Matrix operator-(const S21Matrix& other) const;
    S21Matrix operator*(const S21Matrix& other) const;
    S21Matrix operator*(double num) const;
    bool operator==(const S21Matrix& other) const;
    S21Matrix& operator=(const S21Matrix& other);
    S21Matrix& operator+=(const S21Matrix& other);
    S21Matrix& operator-=(const S21Matrix& other);
    S21Matrix& operator*=(const S21Matrix& other);
    double& operator()(int i, int j);
};

#endif // SRC_S21_MATRIX_OOP_H_