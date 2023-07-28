
#include "s21_matrix_oop.h"

/**
 * @brief Конструктор по умолчанию для класса S21Matrix.
 *
 * Инициализирует пустую матрицу с 0 строками и столбцами.
 */
S21Matrix::S21Matrix() noexcept : rows_(0), cols_(0), matrix_(nullptr) {}

/**
 * @brief Конструктор с параметрами для класса S21Matrix.
 *
 * Создает матрицу заданного размера.
 *
 * @param rows Количество строк в новой матрице.
 * @param cols Количество столбцов в новой матрице.
 */
S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols)
{
  if (rows_ < 0 || cols_ < 0)
  {
    throw std::length_error("Matrix size must be greater or equal than 0");
  }
  matrix_ = new double *[rows_];
  for (int i = 0; i < rows_; i++)
  {
    matrix_[i] = new double[cols_];
  }
}

/**
 * @brief Конструктор копирования для класса S21Matrix.
 *
 * Создает копию другого объекта S21Matrix.
 *
 * @param other Матрица-источник, которую надо скопировать.
 */
S21Matrix::S21Matrix(const S21Matrix &other) : rows_(other.rows_), cols_(other.cols_)
{
  matrix_ = new double *[rows_];
  for (int i = 0; i < rows_; i++)
  {
    matrix_[i] = new double[cols_];
  }
  for (int i = 0; i < rows_; i++)
  {
    std::copy(other.matrix_[i], other.matrix_[i] + cols_, matrix_[i]);
  }
}

/**
 * @brief Конструктор перемещения для класса S21Matrix.
 *
 * Перемещает данные из другого объекта S21Matrix.
 *
 * @param other Матрица, данные которой надо переместить.
 */
S21Matrix::S21Matrix(S21Matrix &&other) noexcept
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_)
{

  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

/**
 * @brief Деструктор для класса S21Matrix.
 *
 * Освобождает память, выделенную под матрицу.
 */
S21Matrix::~S21Matrix() { delete_matrix(*this); }

/**
 * @brief Получение количества строк матрицы.
 *
 * @return int Количество строк матрицы.
 */
int S21Matrix::get_rows() { return rows_; }

/**
 * @brief Установка нового количества строк матрицы.
 *
 * Изменяет размер матрицы, добавляя или удаляя строки в конце.
 *
 * @param newRowCount Новое количество строк матрицы.
 */
void S21Matrix::set_rows(int newRowCount)
{
  if (newRowCount <= 0)
    throw std::out_of_range("Row count should be a positive integer.");

  if (newRowCount != rows_)
  {
    double **newMatrix = new double *[newRowCount];
    for (int i = 0; i < newRowCount; i++)
    {
      newMatrix[i] = new double[cols_];
      for (int j = 0; j < cols_; j++)
      {
        if (i < rows_)
        {
          newMatrix[i][j] = matrix_[i][j];
        }
        else
        {
          newMatrix[i][j] = 0;
        }
      }
    }
    delete_matrix(*this);
    rows_ = newRowCount;
    matrix_ = newMatrix;
  }
}

/**
 * @brief Получение количества столбцов матрицы.
 *
 * @return int Количество столбцов матрицы.
 */
int S21Matrix::get_cols() { return cols_; }

/**
 * @brief Установка нового количества столбцов матрицы.
 *
 * Изменяет размер матрицы, добавляя или удаляя столбцы в конце.
 *
 * @param newColCount Новое количество столбцов матрицы.
 */
void S21Matrix::set_cols(int newColCount)
{
  if (newColCount <= 0)
    throw std::out_of_range("Column count should be a positive integer.");

  if (newColCount != cols_)
  {
    double **newMatrix = new double *[rows_];
    for (int i = 0; i < rows_; i++)
    {
      newMatrix[i] = new double[newColCount];
      for (int j = 0; j < newColCount; j++)
      {
        if (j < cols_)
        {
          newMatrix[i][j] = matrix_[i][j];
        }
        else
        {
          newMatrix[i][j] = 0;
        }
      }
    }
    delete_matrix(*this);
    cols_ = newColCount;
    matrix_ = newMatrix;
  }
}

/**
 * @brief Освобождает память, выделенную под матрицу.
 *
 * @param matrixToDeallocate Матрица, для которой необходимо освободить память.
 */
void S21Matrix::delete_matrix(S21Matrix &matrixToDeallocate)
{
  for (int rowIndex = 0; rowIndex < matrixToDeallocate.rows_; rowIndex++)
  {
    delete[] matrixToDeallocate.matrix_[rowIndex];
  }
  delete[] matrixToDeallocate.matrix_;
}

/**
 * @brief Выделяет память под матрицу.
 *
 * @param matrixToAllocate Матрица, для которой необходимо выделить память.
 */
void S21Matrix::allocateMemoryForMatrix(S21Matrix &matrixToAllocate)
{
  matrixToAllocate.matrix_ = new double *[matrixToAllocate.rows_]();
  for (int rowIndex = 0; rowIndex < matrixToAllocate.rows_; rowIndex++)
  {
    matrixToAllocate.matrix_[rowIndex] = new double[matrixToAllocate.cols_]();
  }
}

/**
 * @brief Копирует содержимое одной матрицы в другую.
 *
 * @param sourceMatrix Матрица, содержимое которой нужно скопировать.
 */
void S21Matrix::copyMatrixContent(const S21Matrix &sourceMatrix)
{
  allocateMemoryForMatrix(*this);
  cols_ = sourceMatrix.cols_;
  rows_ = sourceMatrix.rows_;
  for (int rowIndex = 0; rowIndex < rows_; rowIndex++)
  {
    for (int columnIndex = 0; columnIndex < cols_; columnIndex++)
    {
      matrix_[rowIndex][columnIndex] = sourceMatrix.matrix_[rowIndex][columnIndex];
    }
  }
}

/**
 * @brief Проверяет матрицы на равенство.
 *
 * @param other Матрица для сравнения.
 * @return bool true, если матрицы равны, иначе false.
 */
bool S21Matrix::EqMatrix(const S21Matrix &other)
{
  if (other.cols_ != cols_ || other.rows_ != rows_)
    throw std::out_of_range("The matrices cannot be compared. They are not of the same size.");

  int flag = 0;
  for (int i = 0; i < rows_; i++)
  {
    for (int j = 0; j < cols_; j++)
    {
      if (matrix_[i][j] == other.matrix_[i][j])
        flag = 1;
    }
  }
  return flag;
}

/**
 * @brief Суммирует текущую матрицу с другой матрицей.
 *
 * @param other Матрица для суммирования.
 */
void S21Matrix::SumMatrix(const S21Matrix &other)
{
  if (other.cols_ != cols_ || other.rows_ != rows_)
    throw std::out_of_range("Matrices are not the same size.");

  for (int rowIndex = 0; rowIndex < rows_; rowIndex++)
  {
    for (int colIndex = 0; colIndex < cols_; colIndex++)
    {
      matrix_[rowIndex][colIndex] += other.matrix_[rowIndex][colIndex];
    }
  }
}

/**
 * @brief Вычитает из текущей матрицы другую матрицу.
 *
 * @param other Матрица для вычитания.
 */
void S21Matrix::SubMatrix(const S21Matrix &other)
{
  if (other.cols_ != cols_ || other.rows_ != rows_)
    throw std::out_of_range("Matrices are not the same size.");

  for (int rowIndex = 0; rowIndex < rows_; rowIndex++)
  {
    for (int colIndex = 0; colIndex < cols_; colIndex++)
    {
      matrix_[rowIndex][colIndex] -= other.matrix_[rowIndex][colIndex];
    }
  }
}

/**
 * @brief Умножает текущую матрицу на число.
 *
 * @param num Число для умножения.
 */
void S21Matrix::MulNumber(const double num)
{
  for (int rowIndex = 0; rowIndex < rows_; rowIndex++)
  {
    for (int colIndex = 0; colIndex < cols_; colIndex++)
    {
      matrix_[rowIndex][colIndex] *= num;
    }
  }
}

/**
 * @brief Умножает текущую матрицу на другую матрицу.
 *
 * @param other Матрица для умножения.
 */
void S21Matrix::MulMatrix(const S21Matrix &other)
{
  if (other.rows_ != cols_)
    throw std::out_of_range("Number of columns of the first matrix should be equal to the number of rows of the second matrix.");
  S21Matrix res(rows_, other.cols_);
  for (int rowIndex = 0; rowIndex < rows_; rowIndex++)
  {
    for (int colIndex = 0; colIndex < other.cols_; colIndex++)
    {
      res.matrix_[rowIndex][colIndex] = 0;
      for (int k = 0; k < cols_; k++)
      {
        res.matrix_[rowIndex][colIndex] += matrix_[rowIndex][k] * other.matrix_[k][colIndex];
      }
    }
  }
  *this = res;
}

/**
 * @brief Транспонирует текущую матрицу.
 *
 * @return S21Matrix Транспонированная матрица.
 */
S21Matrix S21Matrix::Transpose()
{
  S21Matrix transposedMatrix(cols_, rows_);
  for (int i = 0; i < rows_; i++)
  {
    for (int j = 0; j < cols_; j++)
    {
      transposedMatrix.matrix_[j][i] = matrix_[i][j];
    }
  }
  return transposedMatrix;
}

/**
 * @brief Вычисляет минор матрицы.
 *
 * @param excluded_row Индекс строки, которую нужно исключить.
 * @param excluded_col Индекс столбца, который нужно исключить.
 * @return S21Matrix Минор матрицы.
 */
S21Matrix S21Matrix::calc_minor(int excluded_row, int excluded_col)
{
  int row_offset = 0, col_offset = 0;
  S21Matrix result(rows_ - 1, cols_ - 1);

  for (int i = 0; i < result.rows_; i++)
  {
    if (i == excluded_row)
      row_offset = 1;

    for (int j = 0; j < result.cols_; j++)
    {
      if (j == excluded_col)
      {
        col_offset = 1;
      }
      result.matrix_[i][j] = matrix_[i + row_offset][j + col_offset];
    }

    col_offset = 0;
  }

  return result;
}

/**
 * @brief Вычисляет детерминант матрицы.
 *
 * @return double Детерминант матрицы.
 */
double S21Matrix::Determinant()
{
  if (rows_ != cols_)
    throw std::out_of_range("Matrix is not square, cannot compute determinant.");

  if (rows_ == 1)
  {
    return matrix_[0][0];
  }

  double result = 0;
  for (int i = 0; i < rows_; i++)
  {
    result += pow(-1, i) * matrix_[0][i] * calc_minor(0, i).Determinant();
  }

  return result;
}

/**
 * @brief Вычисляет матрицу сопряженных.
 *
 * @return S21Matrix Матрица сопряженных.
 */
S21Matrix S21Matrix::CalcComplements()
{
  if (rows_ != cols_)
    throw std::out_of_range("Matrix is not square, cannot compute cofactor matrix.");

  S21Matrix cofactorMatrix(rows_, cols_);
  for (int i = 0; i < rows_; i++)
  {
    for (int j = 0; j < cols_; j++)
    {
      cofactorMatrix.matrix_[i][j] = pow(-1, i + j) * calc_minor(i, j).Determinant();
    }
  }

  return cofactorMatrix;
}

/**
 * @brief Оператор присваивания для матрицы.
 *
 * @param otherMatrix Матрица, которую нужно присвоить.
 * @return S21Matrix& Ссылка на текущую матрицу после присваивания.
 */
S21Matrix &S21Matrix::operator=(S21Matrix otherMatrix)
{
  delete_matrix(*this);
  copyMatrixContent(otherMatrix);
  return *this;
}

/**
 * @brief Вычисляет обратную матрицу.
 *
 * @return S21Matrix Обратная матрица.
 */
S21Matrix S21Matrix::InverseMatrix()
{
  double determinant = Determinant();
  if (std::abs(determinant) < 1e-7)
  {
    throw std::logic_error("Cannot invert matrix with zero determinant");
  }

  S21Matrix inverseMatrix = this->Transpose().CalcComplements();
  for (int i = 0; i < inverseMatrix.rows_; i++)
  {
    for (int j = 0; j < inverseMatrix.cols_; j++)
    {
      inverseMatrix.matrix_[i][j] /= determinant;
    }
  }

  return inverseMatrix;
}

/**
 * @brief Оператор "()" для доступа к элементам матрицы.
 *
 * @param i Индекс строки.
 * @param j Индекс столбца.
 * @return double& Ссылка на элемент матрицы.
 */
double &S21Matrix::operator()(int i, int j)
{
  return matrix_[i][j];
}

/**
 * @brief Оператор сложения двух матриц.
 *
 * @param otherMatrix Матрица для сложения.
 * @return S21Matrix Результат сложения двух матриц.
 */
S21Matrix S21Matrix::operator+(S21Matrix &otherMatrix)
{
  S21Matrix resultMatrix(*this);
  resultMatrix.SumMatrix(otherMatrix);
  return resultMatrix;
}

/**
 * @brief Оператор вычитания двух матриц.
 *
 * @param otherMatrix Матрица для вычитания.
 * @return S21Matrix Результат вычитания двух матриц.
 */
S21Matrix S21Matrix::operator-(S21Matrix &otherMatrix)
{
  S21Matrix resultMatrix(*this);
  resultMatrix.SubMatrix(otherMatrix);
  return resultMatrix;
}

/**
 * @brief Оператор умножения двух матриц.
 *
 * @param otherMatrix Матрица для умножения.
 * @return S21Matrix Результат умножения двух матриц.
 */
S21Matrix S21Matrix::operator*(S21Matrix &otherMatrix)
{
  S21Matrix resultMatrix(*this);
  resultMatrix.MulMatrix(otherMatrix);
  return resultMatrix;
}

/**
 * @brief Оператор умножения матрицы на число.
 *
 * @param scalar Число, на которое умножается матрица.
 * @return S21Matrix Результат умножения матрицы на число.
 */
S21Matrix S21Matrix::operator*(double scalar)
{
  S21Matrix resultMatrix(*this);
  resultMatrix.MulNumber(scalar);
  return resultMatrix;
}

/**
 * @brief Оператор сравнения двух матриц.
 *
 * @param otherMatrix Матрица для сравнения.
 * @return bool Результат сравнения двух матриц (true, если равны).
 */
bool S21Matrix::operator==(S21Matrix &otherMatrix)
{
  return EqMatrix(otherMatrix);
}

/**
 * @brief Оператор "+=" для двух матриц.
 *
 * @param otherMatrix Матрица для сложения с текущей матрицей.
 * @return S21Matrix Обновленная текущая матрица после сложения.
 */
S21Matrix S21Matrix::operator+=(S21Matrix &otherMatrix)
{
  SumMatrix(otherMatrix);
  return *this;
}

/**
 * @brief Оператор "-=" для двух матриц.
 *
 * @param otherMatrix Матрица для вычитания из текущей матрицы.
 * @return S21Matrix Обновленная текущая матрица после вычитания.
 */
S21Matrix S21Matrix::operator-=(S21Matrix &otherMatrix)
{
  SubMatrix(otherMatrix);
  return *this;
}

/**
 * @brief Оператор "*=" для двух матриц.
 *
 * @param otherMatrix Матрица для умножения на текущую матрицу.
 * @return S21Matrix Обновленная текущая матрица после умножения.
 */
S21Matrix S21Matrix::operator*=(S21Matrix &otherMatrix)
{
  MulMatrix(otherMatrix);
  return *this;
}

/**
 * @brief Оператор "*=" для умножения матрицы на число.
 *
 * @param scalar Число, на которое умножается текущая матрица.
 * @return S21Matrix Обновленная текущая матрица после умножения.
 */
S21Matrix S21Matrix::operator*=(double scalar)
{
  MulNumber(scalar);
  return *this;
}
