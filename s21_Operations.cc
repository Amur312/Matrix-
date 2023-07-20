#include "s21_matrix_oop.h"

bool S21Matrix::EqMatrix(const S21Matrix &other) const
{
    if (rows_ != other.rows_ || cols_ != other.cols_)
    {
        return false;
    }
    for (int i = 0; i < rows_; i++)
    {
        for (int j = 0; j < cols_; j++)
        {
            if (matrix_[i][j] != other.matrix_[i][j])
            {
                return false; // Матрицы не равны, если их соответствующие элементы не одинаковы
            }
        }
    }
    return true;
}

void S21Matrix::SumMatrix(const S21Matrix &other)
{
    if (rows_ != other.rows_ || cols_ != other.cols_)
    {
        throw std::invalid_argument("Невозможно сложить матрицы разной размерности!");
    }
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
    {
        throw std::invalid_argument("Невозможно вычесть матрицы разной размерности!");
    }
    for (int i = 0; i < rows_; i++)
    {
        for (int j = 0; j < cols_; j++)
        {
            matrix_[i][j] -= other.matrix_[i][j];
        }
    }
}

void S21Matrix::MulNumber(const double num)
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
    {
        throw std::invalid_argument("Количество столбцов в первой матрице должно быть равно количеству строк во второй матрице!");
    }

    // Создание новой матрицы для хранения результата
    double **result = new double *[rows_];
    for (int i = 0; i < rows_; i++)
    {
        result[i] = new double[other.cols_];
    }
    // Вычисление значений для новой матрицы
    /*
    Каждый элемент новой матрицы вычисляется как сумма произведений соответствующих элементов из строки текущей матрицы и столбца другой матрицы.
    */
    for (int i = 0; i < rows_; i++)
    {
        for (int j = 0; j < other.cols_; j++)
        {
            result[i][j] = 0;
            for (int k = 0; k < cols_; k++)
            {
                result[i][j] += matrix_[i][k] * other.matrix_[j][k];
            }
            /*
            matrix_[i][k] * other.matrix_[k][j] это произведение элемента из i-й строки первой матрицы и элемента из j-го столбца второй матрицы.
            */
        }
    }
    // Удаление старой матрицы и замена ее новой
    for (int i = 0; i < rows_; i++)
    {
        delete[] matrix_[i];
    }
    delete[] matrix_;

    matrix_ = result;
    cols_ = other.cols_;
}

S21Matrix S21Matrix::Transpose()
{
    S21Matrix transposed(cols_, rows_); // Создаем новую матрицу с обратными размерами  с помощью конструктора класса
    /*
        При транспонировании строки становятся столбцами, а столбцы становятся строками.
        https://wiki.loginom.ru/articles/transpose.html
    */
    // Копирование элементов в новую матрицу
    for (int i = 0; i < rows_; ++i)
    {
        for (int j = 0; j < cols_; ++j)
        {
            transposed.matrix_[j][i] = matrix_[i][j]; // Переносим элементы из исходной матрицы в транспонированную
        }
    }

    return transposed;
}

/*
    Определители
    http://publish.sutd.ru/e_books/lin_alg_2013/html/matrix_32.html
*/
S21Matrix S21Matrix::CalcComplements()
{
    if (rows_ != cols_)
    {
        throw std::invalid_argument("Матрица не является квадратной!");
    }
    S21Matrix complements(rows_, cols_);

    for (int i = 0; i < rows_; ++i)
    {
        for (int j = 0; j < cols_; ++j)
        {
            // Создаем подматрицу путем удаления i-й строки и j-го столбца
            S21Matrix submatrix(rows_ - 1, cols_ - 1);

            for (int row = 0; row < rows_; ++row)
            {
                for (int col = 0; col < cols_; ++col)
                {
                    if (row != i && col != j)
                    {
                        int new_row = (row < i) ? row : row - 1;
                        int new_col = (col < j) ? col : col - 1;
                        submatrix.matrix_[new_row][new_col] = matrix_[row][col];
                    }
                }
            }

            // Вычисляем определитель подматрицы и умножаем на (-1)^(i+j)
            complements.matrix_[i][j] = std::pow(-1, i + j) * submatrix.Determinant();
        }
    }

    return complements;
}

double S21Matrix::Determinant()
{
    if (rows_ != cols_)
    {
        throw std::invalid_argument("Матрица не является квадратной!");
    }

    if (rows_ == 1)
    {
        // Определитель одномерной матрицы просто равен ее единственному элементу
        return matrix_[0][0];
    }
    else if (rows_ == 2)
    {
        // Определитель двумерной матрицы можно вычислить напрямую
        return matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
    }
    else
    {
        // Для матриц размером 3x3 и больше используем рекурсивное вычисление через миноры
        /*
            https://prog-cpp.ru/matrix-determinant/
        */
        double det = 0;
        for (int j = 0; j < cols_; ++j)
        {
            // Создаем подматрицу путем удаления первой строки и j-го столбца
            S21Matrix submatrix(rows_ - 1, cols_ - 1);

            for (int row = 1; row < rows_; ++row)
            {
                for (int col = 0; col < cols_; ++col)
                {
                    if (col != j)
                    {
                        int new_col = (col < j) ? col : col - 1;
                        submatrix.matrix_[row - 1][new_col] = matrix_[row][col];
                    }
                }
            }

            // Рекурсивно вычисляем определитель подматрицы и умножаем на соответствующий элемент и знак
            det += matrix_[0][j] * std::pow(-1, j) * submatrix.Determinant();
        }

        return det;
    }
}

/*
    Нахождение обратной матрицы
    https://function-x.ru/return_matrix.html
*/
S21Matrix S21Matrix::InverseMatrix()
{
    if (rows_ != cols_)
    {
        throw std::invalid_argument("Матрица не является квадратной!");
    }
    // Если определитель матрицы близок к нулю кидаем исключение
    double det = this->Determinant();
    if (std::abs(det) < 1e-9)
    {
        throw std::invalid_argument("Определитель матрицы равен нулю, обратной матрицы не существует!");
    }

    S21Matrix cofactors = this->CalcComplements();
    S21Matrix cofactors_transposed = cofactors.Transpose();

    // Делим каждый элемент на определитель исходной матрицы
    for (int i = 0; i < rows_; ++i)
    {
        for (int j = 0; j < cols_; ++j)
        {
            cofactors_transposed.matrix_[i][j] /= det;
        }
    }

    return cofactors_transposed;
}

// Перегрузка оператора сложения
S21Matrix S21Matrix::operator+(const S21Matrix &other) const
{
    S21Matrix result(*this); // Создаем копию текущей матрицы
    result.SumMatrix(other); // Прибавляем другую матрицу
    return result;           // Возвращаем результат
}

// Перегрузка оператора вычитания
S21Matrix S21Matrix::operator-(const S21Matrix &other) const
{
    S21Matrix result(*this); // Создаем копию текущей матрицы
    result.SubMatrix(other); // Вычитаем другую матрицу
    return result;           // Возвращаем результат
}

// Перегрузка оператора умножения
S21Matrix S21Matrix::operator*(const S21Matrix &other) const
{
    S21Matrix result(*this); // Создаем копию текущей матрицы
    result.MulMatrix(other); // Умножаем на другую матрицу
    return result;           // Возвращаем результат
}

// Перегрузка оператора умножения на число
S21Matrix S21Matrix::operator*(double num) const
{
    S21Matrix result(*this); // Создаем копию текущей матрицы
    result.MulNumber(num);   // Умножаем на число
    return result;           // Возвращаем результат
}

// Перегрузка оператора сравнения
bool S21Matrix::operator==(const S21Matrix &other) const
{
    return this->EqMatrix(other);
}

// Перегрузка оператора присваивания
S21Matrix &S21Matrix::operator=(const S21Matrix &other)
{
    if (this != &other)
    { // Проверка на самоприсваивание
        // Очистка текущей матрицы
        for (int i = 0; i < rows_; ++i)
        {
            delete[] matrix_[i];
        }
        delete[] matrix_;

        // Копирование другой матрицы
        rows_ = other.rows_;
        cols_ = other.cols_;
        matrix_ = new double *[rows_];
        for (int i = 0; i < rows_; ++i)
        {
            matrix_[i] = new double[cols_];
            for (int j = 0; j < cols_; ++j)
            {
                matrix_[i][j] = other.matrix_[i][j];
            }
        }
    }
    return *this;
}

// Перегрузка оператора присваивания сложения
S21Matrix &S21Matrix::operator+=(const S21Matrix &other)
{
    this->SumMatrix(other);
    return *this;
}

// Перегрузка оператора присваивания вычитания
S21Matrix &S21Matrix::operator-=(const S21Matrix &other)
{
    this->SubMatrix(other);
    return *this;
}

// Перегрузка оператора присваивания умножения
S21Matrix &S21Matrix::operator*=(const S21Matrix &other)
{
    this->MulMatrix(other);
    return *this;
}

// Перегрузка оператора индексации
double &S21Matrix::operator()(int i, int j)
{
    if (i < 0 || i >= rows_ || j < 0 || j >= cols_)
    {
        throw std::out_of_range("Index out of range");
    }
    return matrix_[i][j];
}
