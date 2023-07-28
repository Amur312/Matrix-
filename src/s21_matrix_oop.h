#include <iostream>
#include <math.h>
#include <stdlib.h>

class S21Matrix
{
private:

  // Члены данных
  int rows_, cols_;
  double **matrix_;

  // Вспомогательные методы для работы с матрицами
  S21Matrix calc_minor(int excluded_row, int excluded_col);
  void copyMatrixContent(const S21Matrix &sourceMatrix);
  void allocateMemoryForMatrix(S21Matrix &matrixToAllocate);
  void delete_matrix(S21Matrix &matrixToDeallocate);

public:

  // Конструкторы, деструктор и операторы присваивания
  S21Matrix() noexcept;
  S21Matrix(int rows, int columns);
  ~S21Matrix();
  S21Matrix(const S21Matrix &other);
  S21Matrix(S21Matrix &&other) noexcept;

  // Геттеры и сеттеры
  int get_rows();
  void set_rows(int newRowCount);
  int get_cols();
  void set_cols(int newColCount);

  // Основные операции над матрицами
  bool EqMatrix(const S21Matrix &other);
  void SumMatrix(const S21Matrix &other);
  void SubMatrix(const S21Matrix &other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix &other);
  S21Matrix Transpose();
  double Determinant();
  S21Matrix CalcComplements();
  S21Matrix InverseMatrix();

  // Операторы для работы с матрицами
  bool operator==(S21Matrix &otherMatrix);
  double &operator()(int i, int j);

  S21Matrix operator+(S21Matrix &otherMatrix);
  S21Matrix operator-(S21Matrix &otherMatrix);
  S21Matrix operator*(S21Matrix &otherMatrix);
  S21Matrix operator*(double num);

  S21Matrix &operator=(S21Matrix otherMatrix);
  S21Matrix operator+=(S21Matrix &otherMatrix);
  S21Matrix operator-=(S21Matrix &otherMatrix);
  S21Matrix operator*=(S21Matrix &otherMatrix);
  S21Matrix operator*=(double num);
};
