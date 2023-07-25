
#include "s21_matrix_oop.h"
#include <gtest/gtest.h>

void fill_matrix(S21Matrix &matrix) {
  int count = 0;
  int rows = matrix.get_rows();
  int cols = matrix.get_cols();
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      matrix(i, j) = count++;
    }
  }
}

TEST(Test, eq_matrix) {
  S21Matrix A(3, 3);
  S21Matrix B(3, 3);
  fill_matrix(A);
  if (A == B) {
    EXPECT_EQ(1, 1);
  }
}

TEST(Test, check_rows) {
  S21Matrix A(3, 3);
  EXPECT_EQ(3, A.get_rows());
}

TEST(Test, check_cols) {
  S21Matrix A(3, 3);
  EXPECT_EQ(3, A.get_cols());
}

TEST(Test, set_rows) {
  S21Matrix A(2, 2);
  A.set_rows(5);
  EXPECT_EQ(5, A.get_rows());
}

TEST(test, set_cols) {
  S21Matrix A(2, 2);
  A.set_cols(5);
  EXPECT_EQ(5, A.get_cols());
}

TEST(test, set_cols_less_1) {
  S21Matrix A(3, 3);
  EXPECT_THROW(A.set_cols(-1), std::logic_error);
}

TEST(test, get_rows) {
  S21Matrix A(5, 2);
  EXPECT_EQ(5, A.get_rows());
}

TEST(test, get_cols) {
  S21Matrix A(5, 2);
  EXPECT_EQ(2, A.get_cols());
}

TEST(test, def_construct) {
  S21Matrix A;
  EXPECT_EQ(0, A.get_rows());
  EXPECT_EQ(0, A.get_cols());
}

TEST(test, construct_param) {
  S21Matrix A(3, 5);
  EXPECT_EQ(3, A.get_rows());
  EXPECT_EQ(5, A.get_cols());
}
//
TEST(test, construct_other) {
  S21Matrix A(2, 6);
  S21Matrix B(A);
  EXPECT_EQ(2, B.get_rows());
  EXPECT_EQ(6, B.get_cols());
  EXPECT_EQ(1, A == B);
}

TEST(test, sum) {
  S21Matrix R(3, 2);
  S21Matrix CHECK(3, 2);
  S21Matrix A(3, 2);
  S21Matrix B(3, 2);
  A(0, 0) = 6;
  A(0, 1) = 6;
  B(0, 0) = 2;
  B(0, 1) = 10;
  CHECK(0, 0) = 8;
  CHECK(0, 1) = 16;
  R = A + B;
  EXPECT_EQ(1, CHECK == R);
}

TEST(Test, plus_throw) {
  S21Matrix A(3, 5);
  S21Matrix B(5, 4);
  EXPECT_THROW(A + B, std::logic_error);
}

TEST(Test, minus_throw) {
  S21Matrix A(3, 5);
  S21Matrix B(5, 4);
  EXPECT_THROW(A - B, std::logic_error);
}

TEST(test, minus) {
  S21Matrix R(3, 2);
  S21Matrix CHECK(3, 2);
  S21Matrix A(3, 2);
  S21Matrix B(3, 2);
  A(0, 0) = 6;
  A(0, 1) = 6;
  B(0, 0) = 2;
  B(0, 1) = 10;
  CHECK(0, 0) = 4;
  CHECK(0, 1) = -4;
  R = A - B;
  EXPECT_EQ(1, CHECK == R);
}

TEST(test, plus_eq) {
  S21Matrix A(1, 2);
  S21Matrix B(1, 2);
  B(0, 0) = 2;
  B(0, 1) = 10;
  A(0, 0) = 0;
  A(0, 1) = 0;
  A += B;
  EXPECT_EQ(1, A == B);
}
//
TEST(test, minus_eq) {
  S21Matrix A(3, 2);
  S21Matrix B(3, 2);
  S21Matrix CHECK(3, 2);
  CHECK(0, 0) = -2;
  CHECK(0, 1) = -10;
  B(0, 0) = 2;
  B(0, 1) = 10;
  A(0, 0) = 0;
  A(0, 1) = 0;
  A -= B;
  EXPECT_EQ(1, A == CHECK);
}

TEST(Test, mul_operator_num) {
  S21Matrix A(3, 4);
  S21Matrix B(3, 4);
  S21Matrix CHECK(3, 4);
  CHECK(0, 0) = 4;
  CHECK(0, 1) = 8;
  CHECK(0, 2) = 12;
  CHECK(0, 3) = 4;
  CHECK(1, 0) = 20;
  CHECK(1, 1) = 24;
  CHECK(1, 2) = 190;
  CHECK(1, 3) = 188;
  CHECK(2, 0) = 36;
  CHECK(2, 1) = 188;
  CHECK(2, 2) = 44;
  CHECK(2, 3) = 196;

  B(0, 0) = 2;
  B(0, 1) = 4;
  B(0, 2) = 6;
  B(0, 3) = 2;
  B(1, 0) = 10;
  B(1, 1) = 12;
  B(1, 2) = 95;
  B(1, 3) = 94;
  B(2, 0) = 18;
  B(2, 1) = 94;
  B(2, 2) = 22;
  B(2, 3) = 98;
  A = B * 2;
  EXPECT_EQ(1, A == CHECK);
}

TEST(Test, mul_operator_num_with_eq) {
  S21Matrix A(2, 2);
  S21Matrix CHECK(2, 2);
  A(0, 0) = 2;
  A(0, 1) = 4;
  CHECK(0, 0) = 10;
  CHECK(0, 1) = 20;
  A *= 5;
  EXPECT_EQ(1, A == CHECK);
}

TEST(Test, mul_matrix_throw) {
  S21Matrix A(3, 1);
  S21Matrix B(3, 3);
  EXPECT_THROW(A * B, std::logic_error);
}

TEST(Test, operation_func) {
  S21Matrix A(3, 3);
  A.operator()(0, 0) = 5;
  EXPECT_EQ(5, A.operator()(0, 0));
}

TEST(test_functional, determinant_zero) {
  const int size = 5;
  S21Matrix m(5, 5);
  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
      m(i, j) = j;
  EXPECT_NEAR(0, m.Determinant(), 1e-06);
}

//
TEST(test_functional, determinant_3x3) {
  const int size = 3;
  S21Matrix m(size, size);

  m(0, 0) = 2;
  m(0, 1) = 3;
  m(0, 2) = 1;
  m(1, 0) = 7;
  m(1, 1) = 4;
  m(1, 2) = 1;
  m(2, 0) = 9;
  m(2, 1) = -2;
  m(2, 2) = 1;

  double res = m.Determinant();
  EXPECT_NEAR(res, -32, 1e-6);
}
//
TEST(test_functional, determinant_2x2) {
  const int size = 2;
  S21Matrix m(size, size);

  m(0, 0) = -5;
  m(0, 1) = -4;
  m(1, 0) = -2;
  m(1, 1) = -3;

  double res = m.Determinant();
  EXPECT_NEAR(res, 7, 1e-6);
}
//
TEST(test_functional, complements_3x3) {
  const int rows = 3;
  const int cols = 3;

  S21Matrix given(rows, cols);
  S21Matrix expected(rows, cols);

  expected(0, 0) = 0;
  expected(0, 1) = 10;
  expected(0, 2) = -20;
  expected(1, 0) = 4;
  expected(1, 1) = -14;
  expected(1, 2) = 8;
  expected(2, 0) = -8;
  expected(2, 1) = -2;
  expected(2, 2) = 4;

  given(0, 0) = 1;
  given(0, 1) = 2;
  given(0, 2) = 3;
  given(1, 0) = 0;
  given(1, 1) = 4;
  given(1, 2) = 2;
  given(2, 0) = 5;
  given(2, 1) = 2;
  given(2, 2) = 1;

  S21Matrix res = given.CalcComplements();
  EXPECT_EQ(1, res == expected);
}
//
TEST(test_functional, complements_3x3_1) {
  const int rows = 3;
  const int cols = 3;

  S21Matrix given(rows, cols);
  S21Matrix expected(rows, cols);

  given(0, 0) = 5;
  given(0, 1) = -1;
  given(0, 2) = 1;
  given(1, 0) = 2;
  given(1, 1) = 3;
  given(1, 2) = 4;
  given(2, 0) = 1;
  given(2, 1) = 0;
  given(2, 2) = 3;

  expected(0, 0) = 9;
  expected(0, 1) = -2;
  expected(0, 2) = -3;
  expected(1, 0) = 3;
  expected(1, 1) = 14;
  expected(1, 2) = -1;
  expected(2, 0) = -7;
  expected(2, 1) = -18;
  expected(2, 2) = 17;

  S21Matrix res = given.CalcComplements();

  EXPECT_EQ(1, res == expected);
}
//
TEST(test_functional, complements_3x3_2) {
  const int rows = 3;
  const int cols = 3;

  S21Matrix given(rows, cols);
  S21Matrix expected(rows, cols);

  given(0, 0) = 1;
  given(0, 1) = 2;
  given(0, 2) = 3;
  given(1, 0) = 2;
  given(1, 1) = 4;
  given(1, 2) = 2;
  given(2, 0) = 5;
  given(2, 1) = 2;
  given(2, 2) = 1;

  expected(0, 0) = 0;
  expected(0, 1) = 10;
  expected(0, 2) = -20;
  expected(1, 0) = 4;
  expected(1, 1) = -14;
  expected(1, 2) = 8;
  expected(2, 0) = -8;
  expected(2, 1) = -2;
  expected(2, 2) = 4;

  S21Matrix res = given.CalcComplements();
  EXPECT_EQ(1, res == expected);
}

TEST(Test, trans) {
  S21Matrix A(3, 4);
  A.operator()(0, 0) = 1;
  A.operator()(0, 1) = 2;
  A.operator()(0, 2) = 3;
  A.operator()(0, 3) = 4;
  A.operator()(1, 0) = 5;
  A.operator()(1, 1) = 6;
  A.operator()(1, 2) = 7;
  A.operator()(1, 3) = 8;
  A.operator()(2, 0) = 9;
  A.operator()(2, 1) = 10;
  A.operator()(2, 2) = 11;
  A.operator()(2, 3) = 12;
  S21Matrix B(4, 3);
  B = A.Transpose();
  EXPECT_EQ(4, B.get_rows());
  EXPECT_EQ(3, B.get_cols());
}

TEST(test_functional, inverese_3x3_1) {
  const uint32_t size = 3;
  S21Matrix given(size, size);
  S21Matrix expected(size, size);
  expected(0, 0) = 1;
  expected(0, 1) = -1;
  expected(0, 2) = 1;
  expected(1, 0) = -38;
  expected(1, 1) = 41;
  expected(1, 2) = -34;
  expected(2, 0) = 27;
  expected(2, 1) = -29;
  expected(2, 2) = 24;

  given(0, 0) = 2.0;
  given(0, 1) = 5.0;
  given(0, 2) = 7.0;
  given(1, 0) = 6.0;
  given(1, 1) = 3.0;
  given(1, 2) = 4.0;
  given(2, 0) = 5.0;
  given(2, 1) = -2.0;
  given(2, 2) = -3.0;
  S21Matrix res = given.InverseMatrix();
  double deter_expect = -1.0;
  double deter_res = res.Determinant();
  EXPECT_NEAR(deter_res, deter_expect, 1e-6);
}

TEST(test_functional, inverse_throw) {
  S21Matrix m(2, 3);
  EXPECT_ANY_THROW(m.InverseMatrix());

  S21Matrix n(2, 2);
  EXPECT_ANY_THROW(n.InverseMatrix());
}

TEST(test, error_set_row) {
  S21Matrix A(1, 5);
  EXPECT_THROW(A.set_rows(0), std::logic_error);
}

TEST(test, error_set_col) {
  S21Matrix A(1, 5);
  EXPECT_THROW(A.set_cols(0), std::logic_error);
}

TEST(test, deter4x4) {
  S21Matrix A(4, 4);
  A(0, 0) = 3;
  A(0, 1) = -3;
  A(0, 2) = -5;
  A(0, 3) = 8;
  A(1, 0) = -3;
  A(1, 1) = 2;
  A(1, 2) = 4;
  A(1, 3) = -6;
  A(2, 0) = 2;
  A(2, 1) = -5;
  A(2, 2) = -7;
  A(2, 3) = 5;
  A(3, 0) = -4;
  A(3, 1) = 3;
  A(3, 2) = 5;
  A(3, 3) = -6;
  EXPECT_EQ(18.0, A.Determinant());
}

TEST(test, mult) {
  S21Matrix A(2, 2);
  S21Matrix B(2, 1);
  S21Matrix R(2, 1);
  A(0, 0) = 0;
  A(0, 1) = 1;
  A(1, 0) = 3;
  A(1, 1) = 4;

  B(0, 0) = 9;
  B(1, 0) = 6;

  R(0, 0) = 6;
  R(1, 0) = 51;

  S21Matrix res = A * B;
  EXPECT_EQ(1, R == res);
}

TEST(test, complem) {
  S21Matrix A(3, 3);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(0, 2) = 3;
  A(1, 0) = 0;
  A(1, 1) = 4;
  A(1, 2) = 2;
  A(2, 0) = 5;
  A(2, 1) = 2;
  A(2, 2) = 1;
  S21Matrix R(3, 3);
  R(0, 0) = 0;
  R(0, 1) = 10;
  R(0, 2) = -20;
  R(1, 0) = 4;
  R(1, 1) = -14;
  R(1, 2) = 8;
  R(2, 0) = -8;
  R(2, 1) = -2;
  R(2, 2) = 4;
  S21Matrix res = A.CalcComplements();
  EXPECT_EQ(1, res == R);
}

int main(int argc, char *argv[]) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
  return 0;
}