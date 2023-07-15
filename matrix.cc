    #include "s21_matrix_oop.h"
    // Базовый конструктор
    S21Matrix::S21Matrix() {
        rows_ = 0;
        cols_ = 0;
        matrix_ = nullptr;
    }
    // Параметризированный конструктор
    S21Matrix::S21Matrix(int rows, int cols): rows_(rows), cols_(cols){   
        matrix_ = new double*[rows_]; // выделения памяти под массив указателей на double 
        // количество указателей равно количеству строк rows_
        for(int i = 0; i < rows_; ++i){ // проходимся по каждому указателю 
            matrix_[i] = new double[cols_]; //  Для каждого указателя выделяется память под массив
            // чисел, который равняется кол-ву столбцов cols_
        }
    }
    // Конструктор копирования
   S21Matrix::S21Matrix(const S21Matrix& other): rows_(other.rows_), cols_(other.cols_){
    matrix_ = new double*[rows_];
    for(int i = 0; i < rows_; ++i){
        matrix_[i] = new double[cols_];
        for(int j = 0; j < cols_; ++j){
            matrix_[i][j] = other.matrix_[i][j]; // Копирование элементов
        }
    }
}

    // Конструктор переноса
    // UPD - && в данном контексте обозначают rvalue ссылку
    S21Matrix::S21Matrix(S21Matrix&& other): rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_){
        other.rows_ = 0;
        other.cols_ = 0;
        other.matrix_ = nullptr;  
    }

    // Деструктор
    /*
        Деструктор в C++ это специальный метод, который вызывается автоматически 
        при уничтожении объекта. Деструкторы в основном используются для освобождения памяти и 
        других ресурсов, которые были выделены для объекта в течение его существования
    */
    S21Matrix::~S21Matrix(){
        for(int i = 0; i < rows_; ++i){
            // В цикле, освобождаем память каждого подмассива matrix_[i]
            delete[] matrix_[i]; 
        }
        // Затем освобождаем память основного массива matrix_, который содержит указатели на double
        delete[] matrix_; 
    }


 S21Matrix S21Matrix::Transpose() {
    S21Matrix transposed(cols_, rows_);  // Создаем новую матрицу с обратными размерами  с помощью конструктора класса
    
    for(int i = 0; i < rows_; ++i) {
        for(int j = 0; j < cols_; ++j) {
            transposed.matrix_[j][i] = matrix_[i][j];  // Переносим элементы из исходной матрицы в транспонированную
        }
    }
    
    return transposed;  // Возвращаем транспонированную матрицу
}
   