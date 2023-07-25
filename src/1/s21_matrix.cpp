
// #include "s21_matrix_oop.h"

// // constructors

// using namespace std;

// S21Matrix::S21Matrix() {
//   rows_ = 0;
//   cols_ = 0;
//   matrix_ = nullptr;
// }

// S21Matrix::S21Matrix(int rows, int columns) {
//   rows_ = rows;
//   cols_ = columns;
//   matrix_ = new double *[rows_];
//   for (int i = 0; i < rows_; i++) {
//     matrix_[i] = new double[cols_];
//   }
// }

// S21Matrix::S21Matrix(const S21Matrix &other) {
//   rows_ = other.rows_;
//   cols_ = other.cols_;
//   copy_matrix(other);
// }

// S21Matrix::S21Matrix(S21Matrix &&A) {
//   rows_ = A.rows_;
//   cols_ = A.cols_;
//   matrix_ = A.matrix_;
//   A.rows_ = 0;
//   A.cols_ = 0;
//   A.matrix_ = nullptr;
// }

// S21Matrix::~S21Matrix() { delete_matrix(*this); }

// // functions

// int S21Matrix::get_rows() { return rows_; }

// void S21Matrix::set_rows(int input) {
//   if (input <= 0)
//     throw std::out_of_range("Incorrect, not positive");
//   if (input != rows_) {
//     double **sol = new double *[input];
//     for (int i = 0; i < input; i++) {
//       sol[i] = new double[cols_];
//       for (int j = 0; j < cols_; j++) {
//         if (i < rows_) {
//           sol[i][j] = matrix_[i][j];
//         } else {
//           sol[i][j] = 0;
//         }
//       }
//     }
//     delete_matrix(*this);
//     rows_ = input;
//     matrix_ = sol;
//   }
// }

// int S21Matrix::get_cols() { return cols_; }

// void S21Matrix::set_cols(int input) {
//   if (input <= 0)
//     throw std::out_of_range("Incorrect, not positive");
//   if (input != cols_) {
//     double **sol = new double *[rows_];
//     for (int i = 0; i < rows_; i++) {
//       sol[i] = new double[input];
//       for (int j = 0; j < input; j++) {
//         if (j < cols_) {
//           sol[i][j] = matrix_[i][j];
//         } else {
//           sol[i][j] = 0;
//         }
//       }
//     }
//     delete_matrix(*this);
//     cols_ = input;
//     matrix_ = sol;
//   }
// }

// void S21Matrix::delete_matrix(S21Matrix &A) {
//   for (int i = 0; i < A.rows_; i++) {
//     delete[] A.matrix_[i];
//   }
//   delete[] A.matrix_;
// }

// void S21Matrix::create_matrix(S21Matrix &A) {
//   A.matrix_ = new double *[rows_]();
//   for (int i = 0; i < rows_; i++) {
//     A.matrix_[i] = new double[cols_]();
//   }
// }

// void S21Matrix::copy_matrix(const S21Matrix &A) {
//   create_matrix(*this);
//   cols_ = A.cols_;
//   rows_ = A.rows_;
//   for (int i = 0; i < rows_; i++) {
//     for (int j = 0; j < cols_; j++) {
//       matrix_[i][j] = A.matrix_[i][j];
//     }
//   }
// }

// bool S21Matrix::EqMatrix(const S21Matrix &other) {
//   if (other.cols_ != cols_ || other.rows_ != rows_)
//     throw std::out_of_range("Incorrect, not same size");

//   int flag = 0;
//   for (int i = 0; i < rows_; i++) {
//     for (int j = 0; j < cols_; j++) {
//       if (matrix_[i][j] == other.matrix_[i][j])
//         flag = 1;
//     }
//   }
//   return flag;
// }

// void S21Matrix::SumMatrix(const S21Matrix &other) {
//   if (other.cols_ != cols_ || other.rows_ != rows_)
//     throw std::out_of_range("Incorrect, not same size");

//   for (int i = 0; i < rows_; i++) {
//     for (int j = 0; j < cols_; j++) {
//       matrix_[i][j] += other.matrix_[i][j];
//     }
//   }
// }

// void S21Matrix::SubMatrix(const S21Matrix &other) {
//   if (other.cols_ != cols_ || other.rows_ != rows_)
//     throw std::out_of_range("Incorrect, not same size");

//   for (int i = 0; i < rows_; i++) {
//     for (int j = 0; j < cols_; j++) {
//       matrix_[i][j] -= other.matrix_[i][j];
//     }
//   }
// }

// void S21Matrix::MulNumber(const double num) {
//   for (int i = 0; i < rows_; i++) {
//     for (int j = 0; j < cols_; j++) {
//       matrix_[i][j] *= num;
//     }
//   }
// }

// void S21Matrix::MulMatrix(const S21Matrix &other) {
//   if (other.rows_ != cols_)
//     throw std::out_of_range("Incorrect input, not same size");
//   S21Matrix res(rows_, other.cols_);
//   for (int i = 0; i < rows_; i++) {
//     for (int j = 0; j < other.cols_; j++) {
//       res.matrix_[i][j] = 0;
//       for (int k = 0; k < cols_; k++) {
//         res.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
//       }
//     }
//   }
//   *this = res;
// }

// S21Matrix S21Matrix::Transpose() {
//   S21Matrix transp(cols_, rows_);
//   for (int i = 0; i < rows_; i++) {
//     for (int j = 0; j < cols_; j++) {
//       transp.matrix_[j][i] = matrix_[i][j];
//     }
//   }
//   return transp;
// }

// S21Matrix S21Matrix::calc_minor(int a, int b) {
//   int flag1 = 0, flag2 = 0;
//   S21Matrix result(rows_ - 1, cols_ - 1);
//   for (int i = 0; i < result.rows_; i++) {
//     if (i == a)
//       flag1 = 1;
//     for (int j = 0; j < result.cols_; j++) {
//       if (j == b)
//         flag2 = 1;
//       result.matrix_[i][j] = matrix_[i + flag1][j + flag2];
//     }
//     flag2 = 0;
//   }
//   return result;
// }

// double S21Matrix::Determinant() {
//   double result = 0;
//   if (rows_ != cols_)
//     throw std::out_of_range("Incorrect input, not same size");
//   if (rows_ == 1) {
//     result = matrix_[0][0];
//   } else {
//     for (int i = 0; i < rows_; i++) {
//       S21Matrix buff;
//       result += pow(-1, i) * matrix_[0][i] * calc_minor(0, i).Determinant();
//     }
//   }
//   return result;
// }

// S21Matrix S21Matrix::CalcComplements() {
//   S21Matrix buff(rows_, cols_);
//   if (rows_ != cols_)
//     throw std::out_of_range("Incorrect input, not same size");
//   for (int i = 0; i < rows_; i++) {
//     for (int j = 0; j < cols_; j++) {
//       buff.matrix_[i][j] = pow(-1, i + j) * calc_minor(i, j).Determinant();
//     }
//   }
//   return buff;
// }

// S21Matrix &S21Matrix::operator=(S21Matrix A) {
//   delete_matrix(*this);
//   copy_matrix(A);
//   return *this;
// }

// S21Matrix S21Matrix::InverseMatrix() {
//   double det = Determinant();
//   if (std::abs(det) < 1e-7) {
//     throw std::logic_error("Determinant must be non-zero");
//   }
//   S21Matrix inv(*this);
//   inv = inv.Transpose().CalcComplements();
//   for (int i = 0; i < inv.rows_; i++) {
//     for (int j = 0; j < inv.cols_; j++) {
//       inv.matrix_[i][j] /= det;
//     }
//   }
//   return inv;
// }

// S21Matrix S21Matrix::operator+(S21Matrix &A) {
//   S21Matrix res(*this);
//   res.SumMatrix(A);
//   return res;
// }

// S21Matrix S21Matrix::operator-(S21Matrix &A) {
//   S21Matrix res(*this);
//   res.SubMatrix(A);
//   return res;
// }

// S21Matrix S21Matrix::operator*(S21Matrix &A) {
//   S21Matrix res(*this);
//   res.MulMatrix(A);
//   return res;
// }

// S21Matrix S21Matrix::operator*(double num) {
//   S21Matrix res(*this);
//   res.MulNumber(num);
//   return res;
// }

// bool S21Matrix::operator==(S21Matrix &A) { return EqMatrix(A); }

// S21Matrix S21Matrix::operator+=(S21Matrix &A) {
//   SumMatrix(A);
//   return *this;
// }

// S21Matrix S21Matrix::operator-=(S21Matrix &A) {
//   SubMatrix(A);
//   return *this;
// }

// S21Matrix S21Matrix::operator*=(S21Matrix &A) {
//   MulMatrix(A);
//   return *this;
// }

// S21Matrix S21Matrix::operator*=(double num) {
//   MulNumber(num);
//   return *this;
// }

// double &S21Matrix::operator()(int i, int j) { return matrix_[i][j]; }