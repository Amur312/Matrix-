#include <iostream>
#include <math.h>
#include <stdlib.h>

class S21Matrix {
private:
 
  int rows_, cols_;
  double **matrix_; 

  S21Matrix calc_minor(int a, int b);
  void copy_matrix(const S21Matrix &A);
  void create_matrix(S21Matrix &A);
  void delete_matrix(S21Matrix &A);
  void validate_rows(int new_rows);
  S21Matrix create_new_matrix(int new_rows);
  void update_data(S21Matrix& new_matrix);

public:
  S21Matrix(int rows, int columns);
  S21Matrix();  // Default constructor
  ~S21Matrix(); // Destructor
  S21Matrix(const S21Matrix &other);
  S21Matrix(S21Matrix &&A);

  int get_rows();
  void set_rows(int input);
  int get_cols();
  void set_cols(int input);

  bool EqMatrix(const S21Matrix &other);
  void SumMatrix(const S21Matrix &other);
  void SubMatrix(const S21Matrix &other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix &other);
  S21Matrix Transpose();
  double Determinant();
  S21Matrix CalcComplements();
  S21Matrix InverseMatrix();

  S21Matrix operator+(S21Matrix &A);
  S21Matrix operator-(S21Matrix &A);
  S21Matrix operator*(S21Matrix &A);
  S21Matrix operator*(double num);
  bool operator==(S21Matrix &A);
  S21Matrix &operator=(S21Matrix A);
  S21Matrix operator+=(S21Matrix &A);
  S21Matrix operator-=(S21Matrix &A);
  S21Matrix operator*=(S21Matrix &A);
  S21Matrix operator*=(double num);
  double &operator()(int i, int j);
};
