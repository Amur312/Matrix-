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
    for(int i =0; i < rows_; ++i){
        matrix_ = new double*[cols_];
        for(int j = 0;j<cols_; ++j){
            matrix_[i][j] = other.matrix_[i][j]; // Копируем
        }
    }
}
// Конструктор переноса
// UPD - && в данном контексте обозначают rvalue ссылку
S21Matrix::S21Matrix(S21Matrix&& other): rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_){
    other.rows_ = 0;
    other.cols_ = 0;
    other.matrix_ = nullptr;
    // После того , как мы передали все ресурсы 
    // мы обнуляем указатель other.matrix_ и задаем размеры other как 0, 
    // чтобы указать, что эти ресурсы теперь принадлежат новому объекту, а не other. 
    // это предотвращает удаление данных при вызове деструктора для other.  
}

// Деструктор
/*
    Деструктор в C++ это специальный метод, который вызывается автоматически 
    при уничтожении объекта. Деструкторы в основном используются для освобождения памяти и 
    других ресурсов, которые были выделены для объекта в течение его существования
*/
S21Matrix::~S21Matrix(){
    for(int i = 0; i < rows_; ++i){
        |// В цикле, освобождаем память каждого подмассива matrix_[i]
        delete[] matrix_[i]; 
    }
    // Затем освобождаем память основного массива matrix_, который содержит указатели на double
     delete[] matrix_; 
}