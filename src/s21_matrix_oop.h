#include <iostream>
#include <math.h>
#include <stdlib.h>

class S21Matrix
{
private:
    int rows_, cols_;
    double **matrix_;

public:
    S21Matrix() noexcept;
    S21Matrix(int rows, int cols);
    S21Matrix(const S21Matrix &other);
    S21Matrix(S21Matrix &&other) noexcept;
    ~S21Matrix() noexcept;
    int get_cols() const noexcept;
    int get_rows() const noexcept;
    void set_rows(int new_rows);
    void set_cols(int new_cols);
    void delete_matrix(S21Matrix &oth);
    bool EqMatrix(const S21Matrix &other) const;
    void SumMatrix(const S21Matrix &other);
    void SubMatrix(const S21Matrix &other);
    void MulNumber(const double number) noexcept;
    void MulMatrix(const S21Matrix &other);
    S21Matrix Transpose() const;
    S21Matrix GetMinorMatrix(int row, int col) const;
    S21Matrix CalcComplements() const;
    double Determinant() const;
    S21Matrix InverseMatrix() const;
    double &operator()(int row, int col) &;
    double &operator()(int row, int col) && = delete;
    const double &operator()(int row, int col) const &;
    const double &operator()(int row, int col) const && = delete;
    bool operator==(const S21Matrix &other) const;
    S21Matrix operator+(const S21Matrix &other) const;
    S21Matrix &operator+=(const S21Matrix &other);
    S21Matrix operator-(const S21Matrix &other) const;
    S21Matrix &operator-=(const S21Matrix &other);
    S21Matrix operator*(double number) const noexcept;
    friend S21Matrix operator*(double number, const S21Matrix &matrix) noexcept;
    S21Matrix &operator*=(double number);
    S21Matrix operator*(const S21Matrix &other) const;
    S21Matrix &operator*=(const S21Matrix &other);
    S21Matrix& operator=(const S21Matrix& other);
};