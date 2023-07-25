#include "s21_matrix_oop.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
S21Matrix::S21Matrix() noexcept : rows_(0), cols_(0), matrix_(nullptr) {}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols)
{
    if (rows_ < 0 || cols_ < 0)
    {
        throw std::length_error("Matrix size must be greater than or equal to 0");
    }
    matrix_ = new double *[rows_];
    for (int i = 0; i < rows_; ++i)
        matrix_[i] = new double[cols_]();
}

S21Matrix::S21Matrix(const S21Matrix &other) : rows_(other.rows_), cols_(other.cols_)
{
    matrix_ = new double *[rows_];
    for (int i = 0; i < rows_; ++i)
    {
        matrix_[i] = new double[cols_];
        for (int j = 0; j < cols_; ++j)
            matrix_[i][j] = other.matrix_[i][j];
    }
}

S21Matrix::S21Matrix(S21Matrix &&other) noexcept : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_)
{
    other.rows_ = 0;
    other.cols_ = 0;
    other.matrix_ = nullptr;
}

S21Matrix::~S21Matrix() noexcept
{
    for (int i = 0; i < rows_; ++i)
        delete[] matrix_[i];
    delete[] matrix_;
}

bool S21Matrix::EqMatrix(const S21Matrix &other) const
{
    if (rows_ != other.rows_ || cols_ != other.cols_)
        return false;
    for (int i = 0; i < rows_; i++)
    {
        for (int j = 0; j < cols_; j++)
        {
            if (matrix_[i][j] == other.matrix_[i][j])
                return false;
        }
    }
    return true;
}

void S21Matrix::SumMatrix(const S21Matrix &other)
{
    if (rows_ != other.rows_ || cols_ != other.cols_)
        throw std::invalid_argument("Matrix dimensions do not match for addition");
    for (int i = 0; i < rows_; i++)
    {
        for (int j = 0; j < cols_; j++)
        {
            matrix_[i][j] += other.matrix_[i][j];
        }
    }
}

void S21Matrix::SubMatrix(const S21Matrix &other)
{
    if (rows_ != other.rows_ || cols_ != other.cols_)
        throw std::invalid_argument("Matrix dimensions do not match for subtraction");

    for (int i = 0; i < rows_; ++i)
    {
        for (int j = 0; j < cols_; ++j)
        {
            matrix_[i][j] -= other.matrix_[i][j];
        }
    }
}

void S21Matrix::MulNumber(const double num) noexcept
{
    for (int i = 0; i < rows_; i++)
    {
        for (int j = 0; j < cols_; j++)
        {
            matrix_[i][j] *= num;
        }
    }
}

void S21Matrix::MulMatrix(const S21Matrix &other)
{
    if (cols_ != other.rows_)
        throw std::invalid_argument("Number of columns of the first matrix does not equal the number of rows of the second matrix");

    S21Matrix result(rows_, other.cols_);
    for (int i = 0; i < rows_; i++)
        for (int j = 0; j < other.cols_; j++)
            for (int k = 0; k < cols_; k++)
                result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];

    *this = std::move(result);
}

S21Matrix S21Matrix::Transpose() const
{
    S21Matrix result(cols_, rows_);
    for (int i = 0; i < rows_; i++)
        for (int j = 0; j < cols_; j++)
            result.matrix_[j][i] = matrix_[i][j];

    return result;
}
S21Matrix S21Matrix::GetMinorMatrix(int row, int col) const
{
    if (row >= rows_ || col >= cols_)
    {
        throw std::invalid_argument("Index outside matrix");
    }

    S21Matrix minor_matrix(rows_ - 1, cols_ - 1);

    int minor_row = 0;
    for (int i = 0; i < rows_; ++i)
    {
        if (i == row)
            continue;
        int minor_col = 0;
        for (int j = 0; j < cols_; ++j)
        {
            if (j == col)
                continue;
            minor_matrix.matrix_[minor_row][minor_col] = matrix_[i][j];
            minor_col++;
        }
        minor_row++;
    }
    return minor_matrix;
}

S21Matrix S21Matrix::CalcComplements() const
{
    if (rows_ != cols_)
    {
        throw std::logic_error("Incorrect matrix size for CalcComplements");
    }
    S21Matrix result(rows_, cols_);
    for (int i = 0; i < rows_; ++i)
    {
        for (int j = 0; j < cols_; ++j)
        {
            S21Matrix minor_matrix = this->GetMinorMatrix(i, j);
            result.matrix_[i][j] = minor_matrix.Determinant();

            // Изменение знака для каждого второго дополнения
            if ((i + j) % 2 != 0)
            {
                result.matrix_[i][j] *= -1;
            }
        }
    }
    return result;
}

double S21Matrix::Determinant() const
{
    if (rows_ != cols_)
        throw std::invalid_argument("The matrix is ​​not square");

    S21Matrix temp(*this);
    const double eps = 1e-10; // некоторое малое значение для проверки близости числа к нулю
    double det = 1.0;

    for (int i = 0; i < rows_; i++)
    {
        int k = i;
        for (int j = i + 1; j < rows_; j++)
        {
            if (abs(temp.matrix_[j][i]) > abs(temp.matrix_[k][i]))
                k = j;
        }

        if (abs(temp.matrix_[k][i]) < eps)
        {
            det = 0.0;
            break;
        }

        std::swap(temp.matrix_[i], temp.matrix_[k]); // меняем строки напрямую

        if (i != k)
            det = -det;

        det *= temp.matrix_[i][i];

        for (int j = i + 1; j < rows_; j++)
            temp.matrix_[i][j] /= temp.matrix_[i][i];

        for (int j = 0; j < rows_; j++)
        {
            if ((j != i) && (abs(temp.matrix_[j][i]) > eps))
            {
                for (int k = i + 1; k < rows_; k++)
                    temp.matrix_[j][k] -= temp.matrix_[i][k] * temp.matrix_[j][i];
            }
        }
    }

    return det;
}

S21Matrix S21Matrix::InverseMatrix() const
{
    const double epsilon = 1e-7;
    double det = Determinant();
    if (std::abs(det) < epsilon)
    {
        throw std::logic_error("Determinant must be non-zero");
    }

    S21Matrix inv = this->Transpose().CalcComplements();
    for (int i = 0; i < inv.rows_; i++)
    {
        for (int j = 0; j < inv.cols_; j++)
        {
            inv.matrix_[i][j] /= det;
        }
    }
    return inv;
}

int S21Matrix::get_rows() const noexcept { return rows_; }
int S21Matrix::get_cols() const noexcept { return cols_; }
void S21Matrix::set_rows(int new_rows)
{
    if (new_rows < 0)
    {
        throw std::length_error("matrix rows count must be non-negative");
    }

    if (new_rows != rows_)
    {
        double **sol = new double *[new_rows];
        for (int i = 0; i < new_rows; i++)
        {
            sol[i] = new double[cols_];
            for (int j = 0; j < cols_; j++)
            {
                if (i < rows_)
                {
                    sol[i][j] = matrix_[i][j];
                }
                else
                {
                    sol[i][j] = 0;
                }
            }
        }
        delete_matrix(*this);
        rows_ = new_rows;
        matrix_ = sol;
    }
}

void S21Matrix::set_cols(int new_cols)
{
    if (new_cols < 0)
    {
        throw std::length_error("matrix columns count must be non-negative");
    }

    if (new_cols != cols_)
    {
        double **sol = new double *[rows_];
        for (int i = 0; i < rows_; i++)
        {
            sol[i] = new double[new_cols];
            for (int j = 0; j < new_cols; j++)
            {
                if (j < cols_)
                {
                    sol[i][j] = matrix_[i][j];
                }
                else
                {
                    sol[i][j] = 0;
                }
            }
        }
        delete_matrix(*this);
        cols_ = new_cols;
        matrix_ = sol;
    }
}
void S21Matrix::delete_matrix(S21Matrix &oth)
{
    for (int i = 0; i < oth.rows_; i++)
    {
        delete[] oth.matrix_[i];
    }
    delete[] oth.matrix_;
}

double &S21Matrix::operator()(int row, int col) &
{
    if (row < 0 || row >= rows_ || col < 0 || col >= cols_)
    {
        throw std::out_of_range("Indexes are out of range");
    }
    return matrix_[row][col];
}

const double &S21Matrix::operator()(int row, int col) const &
{
    if (row < 0 || row >= rows_ || col < 0 || col >= cols_)
    {
        throw std::out_of_range("Indexes are out of range");
    }
    return matrix_[row][col];
}
S21Matrix S21Matrix::operator+(const S21Matrix &other) const
{
    S21Matrix result(*this);
    result.SumMatrix(other);
    return result;
}
S21Matrix &S21Matrix::operator+=(const S21Matrix &other) {
  SumMatrix(other);
  return *this;
}
S21Matrix S21Matrix::operator-(const S21Matrix &other) const
{
    S21Matrix result = *this;
    result.SubMatrix(other);
    return result;
}
S21Matrix &S21Matrix::operator-=(const S21Matrix &other) {
  SubMatrix(other);
  return *this;
}
S21Matrix S21Matrix::operator*(double number) const noexcept
{
    S21Matrix result = *this;
    result.MulNumber(number);
    return result;
}
S21Matrix operator*(double number, const S21Matrix &matrix) noexcept
{
    S21Matrix result = matrix * number;
    return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) const
{
    S21Matrix result = *this;
    result.MulMatrix(other);
    return result;
}

bool S21Matrix::operator==(const S21Matrix &other) const
{
    return this->EqMatrix(other);
}

S21Matrix &S21Matrix::operator*=(double number)
{
    MulNumber(number);
    return *this;
}

S21Matrix &S21Matrix::operator*=(const S21Matrix &other){
    if (cols_ != other.rows_)
        throw std::invalid_argument("Number of columns of the first matrix does not equal the number of rows of the second matrix");
    
    S21Matrix result(rows_, other.cols_);
    for (int i = 0; i < rows_; i++)
        for (int j = 0; j < other.cols_; j++)
            for (int k = 0; k < cols_; k++)
                result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];

    *this = std::move(result);
    return *this;
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other)
{
    if (this != &other) // проверка самоприсваивания
    {
        for (int i = 0; i < rows_; ++i) // освобождаем текущую матрицу
            delete[] matrix_[i];
        delete[] matrix_;

        rows_ = other.rows_;
        cols_ = other.cols_;
        
        matrix_ = new double*[rows_]; // создаем новую матрицу
        for (int i = 0; i < rows_; ++i)
        {
            matrix_[i] = new double[cols_];
            for (int j = 0; j < cols_; ++j)
                matrix_[i][j] = other.matrix_[i][j];
        }
    }
    return *this;
}


