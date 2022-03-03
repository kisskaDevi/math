#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
#include <cmath>
#include <iomanip>

//

template<typename type, int n>
void print(type (&a)[n]);                                           //вывести массив a[n]
template<typename type>
void print(type *a, int n);                                         //вывести массив a

//

template<typename type, int n>
void constprod(type (&a)[n], type constant);                        //умножение каждого элемента массива a[n] на число constant
template<typename type>
void constprod(type *a, type constant, int n);                      //умножение каждого элемента массива a на число constant
template<typename type, int n>
void replace(type (&a)[n],int number1, int number2);                //поменять элементы number1 и number2 местами в массиве a[n]
template<typename type>
void replace(type *a,int number1, int number2, int n);              //поменять элементы number1 и number2 местами в массиве a
template<typename type, int n>
void sum(type (&a)[n],type (&b)[n],type (&res)[n]);                 //поэлементное суммирование массивов a[n] и b[n] в массив res[n]
template<typename type>
void sum(type *a,type *b, type *res, int n);                        //поэлементное суммирование массивов a и b в массив res
template<typename type, int n>
type product(type (&a)[n],type (&b)[n]);                            //скалярное произведение массивов a[n] и b[n]
template<typename type>
type product(type *a,type *b, int n);                               //скалярное произведение массивов a и b

//

template<typename type>
void print(type **a, int n, int m);                                 //вывести массив a
template<typename type, int n, int m>
void print(type (&a)[n][m]);                                        //вывести массив a[n][m]

//

template<typename type, int n, int m>
void string_replace(type (&a)[n][m],int number1, int number2);          //поменять строки number1 и number2 в массиве a[n][m]
template<typename type, int n, int m>
void column_replace(type (&a)[n][m],int number1, int number2);          //поменять столбцы number1 и number2 в массиве a[n][m]
template<typename type, int n, int m>
void inversive(type (&a)[n][m],type (&res)[n][m]);                      //обратная матрица res[n][m] матрицы a[n][m]
template<typename type, int n, int m>
type det(type (&a)[n][m]);                                              //определитель матрицы a[n][m]
template<typename type, int n, int m>
type tr(type (&a)[n][m]);                                               //след матрицы a[n][m]
template<typename type>
void string_replace(type **a,int number1, int number2, int n, int m);   //поменять строки number1 и number2 в массиве a
template<typename type>
void column_replace(type **a,int number1, int number2, int n, int m);   //поменять столбцы number1 и number2 в массиве a
template<typename type>
void inversive(type **a,type **res, int n, int m);                      //обратная матрица res матрицы a
template<typename type>
type det(type **a, int n, int m);                                       //определитель матрицы a
template<typename type>
type tr(type **a, int n, int m);                                        //след матрицы a

//

template<typename type, int n, int m>
void constprod(type (&a)[n][m], type constant);                                 //умножение каждого элемента массива a[n][m] на константу constant
template<typename type, int n1, int m1>
void sum(type (&a)[n1][m1], type (&b)[n1][m1], type (&res)[n1][m1]);            //поэлементное суммирование массивов a[n][m], b[n][m] в массив res[n][m]
template<typename type, int n1, int m1, int n2, int m2>
void product(type (&a)[n1][m1], type (&b)[n2][m2], type (&res)[n1][m2]);        //матричное произведение массивов a[n][m], b[n][m] в массив res[n][m]
template<typename type>
void constprod(type **a, type constant, int n, int m);                          //умножение каждого элемента массива a на константу constant
template<typename type>
void sum(type **a, type **b, type **res, int n1, int m1);                       //поэлементное суммирование массивов a, b в массив res
template<typename type>
void product(type **a, type **b, type **res, int n1, int m1, int n2, int m2);   //матричное произведение массивов a, b в массив res

//

template<typename type>
void Gauss(type **a, int n1, int m1);                           //Решение системы методом Гаусса
template<typename type, int n1, int m1>
void Gauss(type (&a)[n1][m1]);                                  //Решение системы методом Гаусса
template<typename type>
void Simple_Gauss(type **a, int n1, int m1);                    //Решение системы упрощённым методом Гаусса
template<typename type, int n1, int m1>
void Simple_Gauss(type (&a)[n1][m1]);                           //Решение системы упрощённым методом Гаусса
template<typename type>
void Tridiagonal(type **a, int n1, int m1);                     //Решение системы методом прогонки
template<typename type, int n1, int m1>
void Tridiagonal(type (&a)[n1][m1]);                            //Решение системы методом прогонки

//

template<typename type, int n>
void copy(type (&from)[n], type (&to)[n]);                      //Копирование матрицы from в матрицу to
template<typename type>
void copy(type *from, type *to, int n);                         //Копирование матрицы from в матрицу to
template<typename type, int n, int m>
void copy(type (&from)[n][m], type (&to)[n][m]);                //Копирование матрицы from в матрицу to
template<typename type>
void copy(type **from, type **to, int n, int m);                //Копирование матрицы from в матрицу to

//

int matrix_out_setprecision=1;
int matrix_out_setw=5;

template <class type>
class matrix
{
    private:
        type **mtx;                                         //динамический массив, в котором хранится введённая информаци
        int n;                                              //количество строк в матрице
        int m;                                              //количество столбцов в матрице
    public:
        matrix();
        matrix(int n, int m);                                           //создаётся объект - матрица n x m
        matrix(int n);                                                  //создаётся вектор - столбец n x 1 (для создания вектора - строки надо использовать конструктор выше, определив матрицу 1 x m)
        matrix(const matrix<type> & other);                             //конструктор копирования
        matrix(type **MTX, int n, int m);
        template <int n, int m> matrix(const type (&MTX)[n][m])
        {
            this->n=n; this->m=m; this->mtx= new type *[this->n];
            for(int i=0;i<this->n;i++){
                this->mtx[i]=new type [this->m];
                for(int j=0;j<this->m;j++){
                    this->mtx[i][j]=MTX[i][j];}}
        }
        ~matrix();                                                      //деструктор

        matrix &    operator =      (const matrix<type> & other);           //   оператор    приравнивания матриц
        bool        operator ==     (const matrix<type> & other);           //   оператор    равенства возвращает значение true если элементы мариц равны
        bool        operator !=     (const matrix<type> & other);           //   оператор    неравенства - обратное оператора равенства
        type &      operator ()     (const int i,const int j);              //   оператор    позволяет задать элемент (i,j) матрицы
        matrix      operator +      (const matrix<type> & other);           //   оператор    матричной сумма
        matrix      operator -      (const matrix<type> & other);           //   оператор    матричной разность
        matrix      operator +      ();                                     //   оператор    унитарный плюс
        matrix      operator -      ();                                     //   оператор    унитарный минус
        matrix      operator *      (const matrix<type> & other);           //   оператор    произведения двух матриц
        matrix      operator *      (const type & constant);                //   оператор    произведения матрицы на константу
        template<class T1, class T2> friend matrix<T2> operator*(const T1 & constant, const matrix<T2> & other);

        matrix & inversive();                                          //обратная матрица
        matrix & T();                                                  //транспонирование матрицы
        matrix & identity();                                           //единичная матрица
        matrix & replaceString(int number1, int number2);              //поменять строки
        matrix & replaceColumn(int number1, int number2);              //поменять столбцы

        void addString(int number);                                    //   добавить    строки в конец
        void addColumn(int number);                                    //   добавить    столбцы в конец
        void addString(int pos, int number);                           //   добавить    number строк в позицию pos
        void addColumn(int pos, int number);                           //   добавить    number столбцов в позицию pos
        void deleteString(int number);                                 //   удалить     строки в конец
        void deleteColumn(int number);                                 //   удалить     столбцы в конец

        int getn() const;                                              //получить количество строк в матрице
        int getm() const;                                              //получить количество столбцов в матрице

        template<class Type> friend std::ostream &operator<< (std::ostream &out, const matrix<Type> & other);               //оператор вывода
        template<class Type> friend std::istream &operator>> (std::istream &in, matrix<Type> &point);                       //оператор ввода

        template<class Type> friend Type tr(const matrix<Type> & other);                                                    //след матрицы
        template<class Type> friend Type det(const matrix<Type> & other);                                                   //определитель матрицы
        template<class Type> friend matrix<Type> inversive(const matrix<Type> & other);                                     //обратная матрица
        template<class Type> friend matrix<Type> T(const matrix<Type> & other);                                             //транспонирование матрицы
        template<class Type> friend matrix<Type> identity(const matrix<Type> & other);                                      //единичная матрицы
        template<class Type> friend matrix<Type> replaceString(const matrix<Type> & other, int number1, int number2);       //поменять строки
        template<class Type> friend matrix<Type> replaceColumn(const matrix<Type> & other, int number1, int number2);       //поменять столбцы
        template<class Type> friend void copy(const matrix<Type> & from, Type **to);                                        //копирование матрицы в массив

        void gauss();           //решение системы линейных уравнений, записанных в виде матрицы
        void tridiagonal();     //решение системы уравнений, в виде тридиагональной матрцы
};

template <>
class matrix<char>{};
template <>
class matrix<char *>{};

//============================================================================================================================================//
//============================================================================================================================================//
//============================================================================================================================================//

template<typename type, int n>
void print(type (&a)[n])
{
    for (int i=0; i<n; i++)
    {
        std::cout<<std::setw(10)<<a[i];
    }
}

template<typename type>
void print(type *a, int n)
{
    for (int i=0; i<n; i++)
    {
        std::cout<<std::setw(10)<<a[i];
    }
}

template<typename type, int n>
void constprod(type (&a)[n], type constant)
{
    for(int i=0;i<n;i++)
    {
        a[i]=constant*a[i];
    }
}

template<typename type>
void constprod(type *a, type constant, int n)
{
    for(int i=0;i<n;i++)
    {
        a[i]=constant*a[i];
    }
}


template<typename type, int n>
void replace(type (&a)[n],int number1, int number2)
{
    if(number1<n&&number2<n)
    {
        type b=0;
        b=a[number1];
        a[number1]=a[number2];
        a[number2]=b;
    }
}

template<typename type>
void replace(type *a,int number1, int number2, int n)
{
    if(number1<n&&number2<n)
    {
        type b=0;
        b=a[number1];
        a[number1]=a[number2];
        a[number2]=b;
    }
}

template<typename type, int n>
void sum(type (&a)[n], type (&b)[n], type (&res)[n])
{
    for(int i=0;i<n;i++)
    {
        res[i]=a[i]+b[i];
    }
}

template<typename type>
void sum(type *a,type *b,type *res, int n)
{
    for(int i=0;i<n;i++)
    {
        res[i]=a[i]+b[i];
    }
}

template<typename type, int n>
type product(type (&a)[n],type (&b)[n])
{
    type res=0;
    for(int i=0;i<n;i++)
    {
        res=res+a[i]*b[i];
    }
    return res;
}

template<typename type>
type product(type *a,type *b, int n)
{
    type res=0;
    for(int i=0;i<n;i++)
    {
        res=res+a[i]*b[i];
    }
    return res;
}

template<typename type>
void print(type **a, int n, int m)
{
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
            std::cout<<std::setw(10)<<a[i][j];
        }
        std::cout<<std::endl;
    }
}

template<typename type, int n, int m>
void print(type (&a)[n][m])
{
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
            std::cout<<std::setw(10)<<a[i][j];
        }
        std::cout<<std::endl;
    }
}

template<typename type, int n, int m>
void string_replace(type (&a)[n][m],int number1, int number2)
{
    type b[m]={0};
    for(int j=0;j<m;j++)
    {
        b[j]=a[number1][j];
        a[number1][j]=a[number2][j];
        a[number2][j]=b[j];
    }
}

template<typename type, int n, int m>
void column_replace(type (&a)[n][m],int number1, int number2)
{
    type b[n]={0};
    for(int i=0;i<n;i++)
    {
        b[i]=a[i][number1];
        a[i][number1]=a[i][number2];
        a[i][number2]=b[i];
    }
}

template<typename type, int n, int m>
void inversive(type (&a)[n][m],type (&res)[n][m])
{
    type d=det(a);
    if(d!=0)
    {
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<m;j++)
            {
                type b[n-1][m-1]={};
                for(int k=0;k<n-1;k++)
                {
                    if(k<i)
                    {
                        for(int s=0;s<n-1;s++)
                        {
                            if(s<j)
                            {
                                b[k][s]=a[k][s];
                            }
                            else
                            {
                                b[k][s]=a[k][s+1];
                            }
                        }
                    }
                    else
                    {
                        for(int s=0;s<n-1;s++)
                        {
                            if(s<j)
                            {
                                b[k][s]=a[k+1][s];
                            }
                            else
                            {
                                b[k][s]=a[k+1][s+1];
                            }
                        }
                    }
                }
                res[j][i]=pow(-1,i+1+j+1)*det(b)/d;
            }
        }
    }
    else
    {
        std::cout<<"inversive matrix is not exist"<<std::endl;
    }
}

template<typename type, int n, int m>
type det(type (&A)[n][m])
{
    type a[n][m];
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            a[i][j]=A[i][j];
        }
    }
    if(n==m)
    {
        type s=exp(0);
        for (int i=0; i<n; i++)
        {
            if(a[i][i]==0)
            {
                for (int k=i+1; k<n; k++)
                {
                    if(a[k][i]!=0)
                    {
                        type b[m]={0};
                        for (int j=0; j<m; j++)
                        {
                            b[j]=a[k][j];
                            a[k][j]=a[i][j];
                            a[i][j]=b[j];
                        }
                        k=n;
                        s=s*(-exp(0));
                    }
                }
            }
            if(a[i][i]!=0)
            {
                for (int k=n-1; k>i; k--)
                {
                    for (int j=m-1; j>=0; j--)
                    {
                        a[k][j]=a[k][j]-a[k][i]*a[i][j]/a[i][i];
                    }
                }
                s=s*a[i][i];
            }
            else
            {
                s=0;
                i=n;
            }
        }
        return s;
    }
    else
    {
        std::cout<<"dementions must be equal"<<std::endl;
        return 0;
    }
}

template<typename type, int n, int m>
type tr(type (&a)[n][m])
{
    if(n==m)
    {
        double s=0;
        for(int i=0;i<n;i++)
        {
            s=s+a[i][i];
        }
        return s;
    }
    else
    {
        std::cout<<"dementions must be equal"<<std::endl;
        return 0;
    }
}

template<typename type>
void string_replace(type **a,int number1, int number2,int m)
{
    type *b=new type[m];
    for(int j=0;j<m;j++)
    {
        b[j]=a[number1][j];
        a[number1][j]=a[number2][j];
        a[number2][j]=b[j];
    }
    delete [] b;
}

template<typename type>
void column_replace(type **a,int number1, int number2, int n)
{
    type *b=new type[n];
    for(int i=0;i<n;i++)
    {
        b[i]=a[i][number1];
        a[i][number1]=a[i][number2];
        a[i][number2]=b[i];
    }
    delete [] b;
}

template<typename type>
void inversive(type **a,type **res, int n, int m)
{
    type d=det(a,n,m);
    if(d!=0)
    {
        type **b=new type*[n-1];
        for(int i=0;i<n-1;i++)
        {
            b[i]=new type[m-1];
        }
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<m;j++)
            {
                for(int k=0;k<n-1;k++)
                {
                    if(k<i)
                    {
                        for(int s=0;s<n-1;s++)
                        {
                            if(s<j)
                            {
                                b[k][s]=a[k][s];
                            }
                            else
                            {
                                b[k][s]=a[k][s+1];
                            }
                        }
                    }
                    else
                    {
                        for(int s=0;s<n-1;s++)
                        {
                            if(s<j)
                            {
                                b[k][s]=a[k+1][s];
                            }
                            else
                            {
                                b[k][s]=a[k+1][s+1];
                            }
                        }
                    }
                }
                res[j][i]=pow(-1,i+1+j+1)*det(b,n-1,m-1)/d;
            }
        }
        for(int i=0;i<n-1;i++)
        {
            delete [] b[i];
        }
        delete [] b;
    }
    else
    {
        std::cout<<"inversive matrix is not exist"<<std::endl;
    }
}

template<typename type>
type det(type **A, int n, int m)
{
    type **a=new type*[n];
    type *b=new type[m];
    for(int i=0;i<n;i++)
    {
        a[i]=new type[m];
    }
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            a[i][j]=A[i][j];
        }
    }
    if(n==m)
    {
        type s=exp(0);
        for (int i=0; i<n; i++)
        {
            if(a[i][i]==0)
            {
                for (int k=i+1; k<n; k++)
                {
                    if(a[k][i]!=0)
                    {
                        for (int j=0; j<m; j++)
                        {
                            b[j]=a[k][j];
                            a[k][j]=a[i][j];
                            a[i][j]=b[j];
                        }
                        k=n;
                        s=s*(-exp(0));
                    }
                }
            }
            if(a[i][i]!=0)
            {
                for (int k=n-1; k>i; k--)
                {
                    for (int j=m-1; j>=0; j--)
                    {
                        a[k][j]=a[k][j]-a[k][i]*a[i][j]/a[i][i];
                    }
                }
                s=s*a[i][i];
            }
            else
            {
                s=0;
                i=n;
            }
        }
        for(int i=0;i<n;i++)
        {
            delete [] a[i];
        }
        delete [] a;
        delete [] b;
        return s;
    }
    else
    {
        for(int i=0;i<n;i++)
        {
            delete [] a[i];
        }
        delete [] a;
        delete [] b;
        std::cout<<"dementions must be equal"<<std::endl;
        return 0;
    }
}

template<typename type>
type tr(type **a, int n, int m)
{
    if(n==m)
    {
        double s=0;
        for(int i=0;i<n;i++)
        {
            s=s+a[i][i];
        }
        return s;
    }
    else
    {
        std::cout<<"dementions must be equal"<<std::endl;
        return 0;
    }
}

template<typename type, int n, int m>
void constprod(type (&a)[n][m], type constant)
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            a[i][j]=constant*a[i][j];
        }
    }
}

template<typename type, int n1, int m1>
void sum(type (&a)[n1][m1], type (&b)[n1][m1], type (&res)[n1][m1])
{
    for(int i=0;i<n1;i++)
    {
        for(int j=0;j<m1;j++)
        {
            res[i][j]=a[i][j]+b[i][j];
        }
    }
}

template<typename type, int n1, int m1, int n2, int m2>
void product(type (&a)[n1][m1], type (&b)[n2][m2], type (&res)[n1][m2])
{
    if(m1==n2)
    {
        for(int i=0;i<n1;i++)
        {
            for(int j=0;j<m2;j++)
            {
                res[i][j]=0;
                for(int k=0;k<m1;k++)
                {
                    res[i][j]=res[i][j]+a[i][k]*b[k][j];
                }
            }
        }
    }
    else
    {
        std::cout<<"dementions are not equal"<<std::endl;
    }
}

template<typename type>
void constprod(type **a, type constant, int n, int m)
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            a[i][j]=constant*a[i][j];
        }
    }
}

template<typename type>
void sum(type **a, type **b, type **res, int n1, int m1)
{
    for(int i=0;i<n1;i++)
    {
        for(int j=0;j<m1;j++)
        {
            res[i][j]=a[i][j]+b[i][j];
        }
    }
}

template<typename type>
void product(type **a, type **b, type **res, int n1, int m1, int n2, int m2)
{
    if(m1==n2)
    {
        for(int i=0;i<n1;i++)
        {
            for(int j=0;j<m2;j++)
            {
                res[i][j]=0;
                for(int k=0;k<m1;k++)
                {
                    res[i][j]=res[i][j]+a[i][k]*b[k][j];
                }
            }
        }
    }
    else
    {
        std::cout<<"dementions are not equal"<<std::endl;
    }
}

template<typename type>
void Gauss(type **a, int n1, int m1)
{
    int n=n1;
    int m=m1;
    type *b=new type[m];
    for (int i=0; i<n; i++)
    {
        if(a[i][i]==0)
        {
            for (int k=i+1; k<n; k++)
            {
                if(a[k][i]!=0)
                {
                    for (int j=0; j<m; j++)
                    {
                        b[j]=a[k][j];
                        a[k][j]=a[i][j];
                        a[i][j]=b[j];
                    }
                    k=n;
                }
            }
        }
        for (int k=n-1; k>i; k--)
        {
            int s=0;
            for (int j=m-1; j>=0; j--)
            {
                a[k][j]=a[k][j]-a[k][i]*a[i][j]/a[i][i];
                if(a[k][j]==0&&j!=m-1)
                {
                    s++;
                }
            }
            if(s==m-1 && a[k][m-1]!=0)
            {
                std::cout<<"solution cunnot be founded\n"<<std::endl;
            }
            if(s==m-1 && a[k][m-1]==0)
            {
                for (int j=0; j<m; j++)
                {
                    a[k][j]=a[n-1][j];
                    a[n-1][j]=b[j];
                }
                n--;
            }
        }
    }
    for (int i=n-1; i>=0; i--)
    {
        if(a[i][i]!=0)
        {
            for(int j=m-1; j>=i; j--)
            {
                a[i][j]= a[i][j]/a[i][i];
                if(i>0)
                {
                    for(int k=i-1; k>=0; k--)
                    {
                        a[k][j]=a[k][j]-a[i][j]*a[k][i];
                    }
                }
            }
        }
    }
    delete [] b;
}

template<typename type, int n1, int m1>
void Gauss(type (&a)[n1][m1])
{
    int n=n1;
    int m=m1;
    for (int i=0; i<n; i++)
    {
        if(a[i][i]==0)
        {
            for (int k=i+1; k<n; k++)
            {
                if(a[k][i]!=0)
                {
                    type b[m]={0};
                    for (int j=0; j<m; j++)
                    {
                        b[j]=a[k][j];
                        a[k][j]=a[i][j];
                        a[i][j]=b[j];
                    }
                    k=n;
                }
            }
        }
        for (int k=n-1; k>i; k--)
        {
            int s=0;
            for (int j=m-1; j>=0; j--)
            {
                a[k][j]=a[k][j]-a[k][i]*a[i][j]/a[i][i];
                if(a[k][j]==0&&j!=m-1)
                {
                    s++;
                }
            }
            if(s==m-1 && a[k][m-1]!=0)
            {
                std::cout<<"solution cunnot be founded\n"<<std::endl;
            }
            if(s==m-1 && a[k][m-1]==0)
            {
                type b[m]={0};
                for (int j=0; j<m; j++)
                {
                    a[k][j]=a[n-1][j];
                    a[n-1][j]=b[j];
                }
                n--;
            }
        }
    }
    for (int i=n-1; i>=0; i--)
    {
        if(a[i][i]!=0)
        {
            for(int j=m-1; j>=i; j--)
            {
                a[i][j]= a[i][j]/a[i][i];
                if(i>0)
                {
                    for(int k=i-1; k>=0; k--)
                    {
                        a[k][j]=a[k][j]-a[i][j]*a[k][i];
                    }
                }
            }
        }
    }
}

template<typename type>
void Simple_Gauss(type **a, int n1, int m1)
{
    int n=n1;
    int m=m1;
    type *b=new type[m];
    for (int i=0; i<n; i++)
    {
        for (int k=n-1; k>i; k--)
        {
            if(a[k][i]!=0)
            {
                for (int j=m-1; j>=0; j--)
                {
                    if(a[i][j]!=0)
                    {
                        a[k][j]=a[k][j]-a[k][i]*a[i][j]/a[i][i];
                    }
                }
            }
        }
    }
    for (int i=n-1; i>=0; i--)
    {
        for(int j=m-1; j>=i; j--)
        {
            a[i][j]= a[i][j]/a[i][i];
            if(a[i][j]!=0)
            {
                for(int k=i-1; k>=0; k--)
                {
                    if(a[k][i]!=0)
                    {
                        a[k][j]=a[k][j]-a[i][j]*a[k][i];
                    }
                }
            }
        }
    }
    int i=0;
    for(int j=m-1; j>=i; j--)
    {
        a[i][j]= a[i][j]/a[i][i];
        if(a[i][j]!=0)
        {
            for(int k=i-1; k>=0; k--)
            {
                if(a[k][i]!=0)
                {
                    a[k][j]=a[k][j]-a[i][j]*a[k][i];
                }
            }
        }
    }
    delete [] b;
}

template<typename type, int n1, int m1>
void Simple_Gauss(type (&a)[n1][m1])
{
    int n=n1;
    int m=m1;
    for (int i=0; i<n; i++)
    {
        for (int k=n-1; k>i; k--)
        {
            if(a[k][i]!=0)
            {
                for (int j=m-1; j>=0; j--)
                {
                    if(a[i][j]!=0)
                    {
                        a[k][j]=a[k][j]-a[k][i]*a[i][j]/a[i][i];
                    }
                }
            }
        }
    }
    for (int i=n-1; i>=0; i--)
    {
        for(int j=m-1; j>=i; j--)
        {
            a[i][j]= a[i][j]/a[i][i];
            if(a[i][j]!=0)
            {
                for(int k=i-1; k>=0; k--)
                {
                    if(a[k][i]!=0)
                    {
                        a[k][j]=a[k][j]-a[i][j]*a[k][i];
                    }
                }
            }
        }
    }
    int i=0;
    for(int j=m-1; j>=i; j--)
    {
        a[i][j]= a[i][j]/a[i][i];
        if(a[i][j]!=0)
        {
            for(int k=i-1; k>=0; k--)
            {
                if(a[k][i]!=0)
                {
                    a[k][j]=a[k][j]-a[i][j]*a[k][i];
                }
            }
        }
    }
}

template<typename type>
void Tridiagonal(type **a, int n, int m)
{
    type *A= new type[n];
    type *B= new type[n];
    type *x= new type[n];
    A[0]=-a[0][1]/a[0][0];
    B[0]=a[0][m-1]/a[0][0];
    for(int i=1;i<n-1;i++)
    {
        A[i]=-a[i][i+1]/(a[i][i-1]*A[i-1]+a[i][i]);
        B[i]=(a[i][m-1]-a[i][i-1]*B[i-1])/(a[i][i-1]*A[i-1]+a[i][i]);
    }
    x[n-1]=(a[n-1][m-1]-a[n-1][m-3]*B[n-2])/(a[n-1][m-3]*A[n-2]+a[n-1][m-2]);
    for(int i=n-2;i>=0;i--)
    {
        x[i]=A[i]*x[i+1]+B[i];
    }
    a[0][m-1]=x[0];
    a[0][0]=1;
    a[0][1]=0;
    for(int i=1;i<n-1;i++)
    {
        a[i][m-1]=x[i];
        a[i][i]=1;
        a[i][i-1]=0;
        a[i][i+1]=0;
    }
    a[n-1][m-1]=x[n-1];
    a[n-1][n-1]=1;
    a[n-1][n-2]=0;
    delete [] A;
    delete [] B;
    delete [] x;
}

template<typename type, int n, int m>
void Tridiagonal(type (&a)[n][m])
{
    type A[n]={0};
    type B[n]={0};
    type x[n]={0};
    A[0]=-a[0][1]/a[0][0];
    B[0]=a[0][m-1]/a[0][0];
    for(int i=1;i<n-1;i++)
    {
        A[i]=-a[i][i+1]/(a[i][i-1]*A[i-1]+a[i][i]);
        B[i]=(a[i][m-1]-a[i][i-1]*B[i-1])/(a[i][i-1]*A[i-1]+a[i][i]);
    }
    x[n-1]=(a[n-1][m-1]-a[n-1][m-3]*B[n-2])/(a[n-1][m-3]*A[n-2]+a[n-1][m-2]);
    for(int i=n-2;i>=0;i--)
    {
        x[i]=A[i]*x[i+1]+B[i];
    }
    a[0][m-1]=x[0];
    a[0][0]=1;
    a[0][1]=0;
    for(int i=1;i<n-1;i++)
    {
        a[i][m-1]=x[i];
        a[i][i]=1;
        a[i][i-1]=0;
        a[i][i+1]=0;
    }
    a[n-1][m-1]=x[n-1];
    a[n-1][n-1]=1;
    a[n-1][n-2]=0;
}


template<typename type, int n>
void copy(type (&from)[n], type (&to)[n])
{
    for(int i=0;i<n;i++)
    {
        to[i]=from[i];
    }
}

template<typename type>
void copy(type *from, type *to, int n)
{
    for(int i=0;i<n;i++)
    {
        to[i]=from[i];
    }
}

template<typename type, int n, int m>
void copy(type (&from)[n][m], type (&to)[n][m])
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            to[i][j]=from[i][j];
        }
    }
}

template<typename type>
void copy(type **from, type **to, int n, int m)
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            to[i][j]=from[i][j];
        }
    }
}

//============================================================================================================================================//
//============================================================================================================================================//
//============================================================================================================================================//

template<typename type>
matrix<type>::matrix()
{
    this->mtx = nullptr;
}

template <typename type>
matrix<type>::matrix(int n, int m)
{
    this->n=n;
    this->m=m;
    this->mtx=new type *[n];
    for(int i=0; i<n;i++)
    {
        this->mtx[i]=new type [m];
        for(int j=0; j<m;j++)
        {
            this->mtx[i][j]=0;
        }
    }
}

template <typename type>
matrix<type>::matrix(int n)
{
    this->n=n;
    this->m=1;
    this->mtx=new type *[n];
    for(int i=0; i<n;i++)
    {
        this->mtx[i]=new type [m];
        for(int j=0; j<m;j++)
        {
            this->mtx[i][j]=0;
        }
    }
}

template<class type>
matrix<type>::matrix(const matrix<type> &other)
{
    this->n=other.n;
    this->m=other.m;

    this->mtx= new type *[this->n];
    for(int i=0; i<n;i++)
    {
        this->mtx[i]=new type [this->m];
    }
    for(int i=0;i<this->n;i++)
    {
        for(int j=0;j<this->m;j++)
        {
            this->mtx[i][j]=other.mtx[i][j];
        }
    }
}

template<class type>
matrix<type>::matrix(type **MTX, int n, int m)
{
    this->n=n;
    this->m=m;

    this->mtx= new type *[this->n];
    for(int i=0; i<n;i++)
    {
        this->mtx[i]=new type [this->m];
    }
    for(int i=0;i<this->n;i++)
    {
        for(int j=0;j<this->m;j++)
        {
            this->mtx[i][j]=MTX[i][j];
        }
    }
}

template <typename type>
matrix<type>::~matrix()
{
    for(int i=0; i<n;i++)
    {
        delete[] mtx[i];
    }
    delete[] mtx;
}

template <typename type>
matrix<type> & matrix<type>:: operator=(const matrix<type> & other)
{
    if(this->mtx!=nullptr)
    {
        for(int i=0; i<this->n;i++)
        {
            delete[] this->mtx[i];
        }
        delete[] this->mtx;
    }
    this->n=other.n;
    this->m=other.m;
    this->mtx= new type *[other.n];
    for(int i=0; i<other.n;i++)
    {
        this->mtx[i]=new type [other.m];
    }
    for(int i=0;i<other.n;i++)
    {
        for(int j=0;j<other.m;j++)
        {
            this->mtx[i][j]=other.mtx[i][j];
        }
    }
    return *this;
}

template <typename type>
bool matrix<type>::operator == (const matrix<type> & other)
{
    for(int i=0; i<other.n;i++)
    {
        for(int j=0; j<other.m;j++)
        {
            if(this->mtx[i][j]==other.mtx[i][j]){}
            else{return false;}
        }
    }
    return true;
}

template <typename type>
bool matrix<type>::operator != (const matrix<type> & other)
{
    for(int i=0; i<other.n;i++)
    {
        for(int j=0; j<other.m;j++)
        {
            if(this->mtx[i][j]==other.mtx[i][j]){}
            else{return true;}
        }
    }
    return false;
}

template <typename type>
type & matrix<type>:: operator()(const int i,const int j)
{
    return mtx[i][j];
}

template <typename type>
matrix<type> matrix<type>:: operator+(const matrix<type> & other)
{
    matrix<type> res(*this);
    if(this->n==other.n&&this->m==other.m)
    {
        for(int i=0;i<other.n;i++)
        {
            for(int j=0;j<other.m;j++)
            {
                res.mtx[i][j]=this->mtx[i][j]+other.mtx[i][j];
            }
        }
    }
    return res;
}

template <typename type>
matrix<type> matrix<type>:: operator-(const matrix<type> & other)
{
    matrix<type> res(*this);
    if(this->n==other.n&&this->m==other.m)
    {
        for(int i=0;i<other.n;i++)
        {
            for(int j=0;j<other.m;j++)
            {
                res.mtx[i][j]=this->mtx[i][j]-other.mtx[i][j];
            }
        }
    }
    return res;
}

template<class type>
matrix<type> matrix<type>::operator+()
{
    matrix<type> res(*this);
    return res;
}

template<class type>
matrix<type> matrix<type>::operator-()
{
    matrix<type> res(*this);
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            res.mtx[i][j]=-mtx[i][j];
        }
    }
    return res;
}

template <typename type>
matrix<type> matrix<type>:: operator*(const matrix<type> & other)
{
    matrix<type> res(this->n,other.m);
    if(this->m==other.n)
    {
        for(int i=0;i<this->n;i++)
        {
            for(int j=0;j<other.m;j++)
            {
                res.mtx[i][j]=this->mtx[0][0]*0.0;
                for(int k=0;k<other.n;k++)
                {
                    res.mtx[i][j]=res.mtx[i][j]+this->mtx[i][k]*other.mtx[k][j];
                }
            }
        }
    }
    return res;
}

template <typename type>
matrix<type> matrix<type>:: operator*(const type & constant)
{
    matrix<type> res(*this);
    for(int i=0;i<this->n;i++)
    {
        for(int j=0;j<this->m;j++)
        {
            res.mtx[i][j]=constant*this->mtx[i][j];
        }
    }
    return res;
}

template<class T1, class T2>
matrix<T2> operator*(const T1 & constant, const matrix<T2> & other)
{

    matrix<T2> res(other);
    for(int i=0;i<other.n;i++)
    {
        for(int j=0;j<other.m;j++)
        {
            res.mtx[i][j]*=constant;
        }
    }
    return res;
}

template<class type>
void matrix<type>::addString(int number)
{
    type **newmtx = new type*[n+number];
    for(int i=0;i<n;i++)
    {
        newmtx[i] = new type[m];
        for(int j=0;j<m;j++)
        {
            newmtx[i][j] = mtx[i][j];
        }
    }
    for(int i=n;i<n+number;i++)
    {
        newmtx[i] = new type[m];
        for(int j=0;j<m;j++)
        {
            newmtx[i][j] = 0;
        }
    }
    for(int i=0;i<n;i++)
    {
        delete [] mtx[i];
    }
    delete [] mtx;

    mtx = newmtx;

    n=n+number;
}

template<class type>
void matrix<type>::addColumn(int number)
{
    type **newmtx = new type*[n];
    for(int i=0;i<n;i++)
    {
        newmtx[i] = new type[m+number];
        for(int j=0;j<m;j++)
        {
            newmtx[i][j] = mtx[i][j];
        }
        for(int j=m;j<m+number;j++)
        {
            newmtx[i][j] = 0;
        }
    }
    for(int i=0;i<n;i++)
    {
        delete [] mtx[i];
    }
    delete [] mtx;

    mtx = newmtx;

    m=m+number;
}

template<class type>
void matrix<type>::addString(int pos, int number)
{
    type **newmtx = new type*[n+number];
    for(int i=0;i<pos;i++)
    {
        newmtx[i] = new type[m];
        for(int j=0;j<m;j++)
        {
            newmtx[i][j] = mtx[i][j];
        }
    }
    for(int i=pos;i<pos+number;i++)
    {
        newmtx[i] = new type[m];
        for(int j=0;j<m;j++)
        {
            newmtx[i][j] = 0;
        }
    }
    for(int i=pos+number;i<n+number;i++)
    {
        newmtx[i] = new type[m];
        for(int j=0;j<m;j++)
        {
            newmtx[i][j] = mtx[i-number][j];
        }
    }
    for(int i=0;i<n;i++)
    {
        delete [] mtx[i];
    }
    delete [] mtx;

    mtx = newmtx;

    n=n+number;
}

template<class type>
void matrix<type>::addColumn(int pos, int number)
{
    type **newmtx = new type*[n];
    for(int i=0;i<n;i++)
    {
        newmtx[i] = new type[m+number];
        for(int j=0;j<pos;j++)
        {
            newmtx[i][j] = mtx[i][j];
        }
        for(int j=pos;j<pos+number;j++)
        {
            newmtx[i][j] = 0;
        }
        for(int j=pos+number;j<m+number;j++)
        {
            newmtx[i][j] = mtx[i][j-number];
        }
    }
    for(int i=0;i<n;i++)
    {
        delete [] mtx[i];
    }
    delete [] mtx;

    mtx = newmtx;

    m=m+number;
}

template<class type>
void matrix<type>::deleteString(int number)
{
    type **newmtx = new type*[n-1];
    for(int i=0;i<number;i++)
    {
        newmtx[i] = new type[m];
        for(int j=0;j<m;j++)
        {
            newmtx[i][j] = mtx[i][j];
        }
    }
    for(int i=number;i<n-1;i++)
    {
        newmtx[i] = new type[m];
        for(int j=0;j<m;j++)
        {
            newmtx[i][j] = mtx[i+1][j];
        }
    }
    for(int i=0;i<n;i++)
    {
        delete [] mtx[i];
    }
    delete [] mtx;

    mtx = newmtx;

    n=n-1;
}

template<class type>
void matrix<type>::deleteColumn(int number)
{
    type **newmtx = new type*[n];
    for(int i=0;i<n;i++)
    {
        newmtx[i] = new type[m-1];
        for(int j=0;j<number;j++)
        {
            newmtx[i][j] = mtx[i][j];
        }
        for(int j=number;j<m-1;j++)
        {
            newmtx[i][j] = mtx[i][j+1];
        }
    }
    for(int i=0;i<n;i++)
    {
        delete [] mtx[i];
    }
    delete [] mtx;

    mtx = newmtx;

    m=m-1;
}

template<class type>
int matrix<type>::getn() const
{
    return n;
}

template<class type>
int matrix<type>::getm() const
{
    return m;
}

template <typename type>
matrix<type> & matrix<type>::inversive()
{
    matrix<type> res(*this);
    type d=det(*this);
    if(d!=0)
    {
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<m;j++)
            {
                type **b;
                b=new type *[n-1];
                for(int i=0; i<n-1;i++)
                {
                    b[i]=new type [m-1];
                }
                for(int k=0;k<n-1;k++)
                {
                    if(k<i)
                    {
                        for(int s=0;s<n-1;s++)
                        {
                            if(s<j)
                            {
                                b[k][s]=mtx[k][s];
                            }
                            else
                            {
                                b[k][s]=mtx[k][s+1];
                            }
                        }
                    }
                    else
                    {
                        for(int s=0;s<n-1;s++)
                        {
                            if(s<j)
                            {
                                b[k][s]=mtx[k+1][s];
                            }
                            else
                            {
                                b[k][s]=mtx[k+1][s+1];
                            }
                        }
                    }
                }
                matrix<type> b1(b,n-1,m-1);
                res.mtx[j][i]=det(b1)/d*pow(-1,i+1+j+1);
                for(int i=0; i<n-1;i++)
                {
                    delete[] b[i];
                }
                delete[] b;
            }
        }
    }
    else
    {
        std::cout<<"inversive matrix is not exist"<<std::endl;
    }
    *this=res;
    return *this;
}

template <typename type>
matrix<type> & matrix<type>::T()
{
    matrix<type> res(m,n);
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n;j++)
        {
            res.mtx[i][j]=mtx[j][i];
        }
    }
    *this=res;
    return *this;
}

template<class type>
matrix<type> & matrix<type>::identity()
{
    if(n==m)
    {
        for(int i=0;i<n;i++)
        {
            this->mtx[i][i]=1.0;
            for(int j=i+1;j<n;j++)
            {
                this->mtx[i][j]=0.0;
                this->mtx[j][i]=0.0;
            }
        }
    }
    return *this;
}

template <typename type>
matrix<type> & matrix<type>::replaceString(int number1, int number2)
{
    type save;
    for(int j=0;j<m;j++)
    {
        save=mtx[number1][j];
        mtx[number1][j]=mtx[number2][j];
        mtx[number2][j]=save;
    }
    return *this;
}

template <typename type>
matrix<type> & matrix<type>::replaceColumn(int number1, int number2)
{
    type save;
    for(int j=0;j<n;j++)
    {
        save=mtx[j][number1];
        mtx[j][number1]=mtx[j][number2];
        mtx[j][number2]=save;
    }
    return *this;
}

template<class T>
std::ostream & operator<<(std::ostream &out, const matrix<T> &other)
{
    for(int i=0;i<other.n;i++)
    {
        for(int j=0;j<other.m;j++)
        {
            out<<std::setw(matrix_out_setw)<<std::fixed<<std::setprecision(matrix_out_setprecision)<<other.mtx[i][j];
        }
        out<<std::endl;
    }
    return out;
}

template<class T>
std::istream & operator>>(std::istream &in, matrix<T> &point)
{
    for(int i=0;i<point.n;i++)
    {
        for(int j=0;j<point.m;j++)
        {
            std::cout<<"mtx["<<i<<"]["<<j<<"]=";
            in >>point.mtx[i][j];
        }
    }
    return in;
}

template<class Type>
matrix<Type> T(const matrix<Type> & other)
{
    matrix<Type> res(other.m,other.n);
    for(int i=0;i<other.m;i++)
    {
        for(int j=0;j<other.n;j++)
        {
            res.mtx[i][j]=other.mtx[j][i];
        }
    }
    return res;
}

template<class Type>
matrix<Type> replaceString(const matrix<Type> & other, int number1, int number2)
{
    Type save;
    matrix<Type> res(other);
    res=other;
    for(int j=0;j<other.m;j++)
    {
        save=res.mtx[number1][j];
        res.mtx[number1][j]=res.mtx[number2][j];
        res.mtx[number2][j]=save;
    }
    return res;
}

template<class Type>
matrix<Type> replaceColumn(const matrix<Type> & other, int number1, int number2)
{
    Type save;
    matrix<Type> res(other);
    res=other;
    for(int j=0;j<other.n;j++)
    {
        save=res.mtx[j][number1];
        res.mtx[j][number1]=res.mtx[j][number2];
        res.mtx[j][number2]=save;
    }
    return res;
}

template<class Type>
Type tr(const matrix<Type> & other)
{
    Type res=other.mtx[0][0]*0.0;
    if(other.n==other.m)
    {
        for(int j=0;j<other.n;j++)
        {
            res=res+other.mtx[j][j];
        }
    }
    return res;
}

template<class Type>
Type det(const matrix<Type> & other)
{
    if(other.n==other.m)
    {
        Type *b;
        b=new Type [other.m];
        Type res=other.mtx[0][0]/other.mtx[0][0];
        for (int i=0; i<other.n; i++)
        {
            if(other.mtx[i][i]==0)
            {
                for (int k=i+1; k<other.n; k++)
                {
                    if(other.mtx[k][i]!=0)
                    {
                        for (int j=0; j<other.m; j++)
                        {
                            b[j]=other.mtx[k][j];
                            other.mtx[k][j]=other.mtx[i][j];
                            other.mtx[i][j]=b[j];
                        }
                        k=other.n;
                        res=res*(-exp(0));
                    }
                }
            }
            if(other.mtx[i][i]!=0)
            {
                for (int k=other.n-1; k>i; k--)
                {
                    for (int j=other.m-1; j>=0; j--)
                    {
                        other.mtx[k][j]=other.mtx[k][j]-other.mtx[k][i]*other.mtx[i][j]/other.mtx[i][i];
                    }
                }
                res=res*other.mtx[i][i];
            }
            else
            {
                res=0;
                i=other.n;
            }
        }
        delete [] b;
        return res;
    }
    else
    {
        std::cout<<"dementions must be equal"<<std::endl;
        return 0;
    }
}

template<class Type>
matrix<Type> identity(const matrix<Type> & other)
{
    matrix<Type> res(other);
    if(res.n==res.m)
    {
        for(int i=0;i<res.n;i++)
        {
            res.mtx[i][i]=1.0;
            for(int j=i+1;j<res.n;j++)
            {
                res.mtx[i][j]=0.0;
                res.mtx[j][i]=0.0;
            }
        }
    }
    return res;
}

template<class Type>
matrix<Type> inversive(const matrix<Type> & other)
{
    matrix<Type> res(other);
    Type d=det(other);
    if(d!=0)
    {
        for(int i=0;i<other.n;i++)
        {
            for(int j=0;j<other.m;j++)
            {
                Type **b;
                b=new Type *[other.n-1];
                for(int i=0; i<other.n-1;i++)
                {
                    b[i]=new Type [other.m-1];
                }
                for(int k=0;k<other.n-1;k++)
                {
                    if(k<i)
                    {
                        for(int s=0;s<other.n-1;s++)
                        {
                            if(s<j)
                            {
                                b[k][s]=other.mtx[k][s];
                            }
                            else
                            {
                                b[k][s]=other.mtx[k][s+1];
                            }
                        }
                    }
                    else
                    {
                        for(int s=0;s<other.n-1;s++)
                        {
                            if(s<j)
                            {
                                b[k][s]=other.mtx[k+1][s];
                            }
                            else
                            {
                                b[k][s]=other.mtx[k+1][s+1];
                            }
                        }
                    }
                }
                matrix<Type> b1(b,other.n-1,other.m-1);
                res.mtx[j][i]=det(b1)/d*pow(-1,i+1+j+1);
                for(int i=0; i<other.n-1;i++)
                {
                    delete[] b[i];
                }
                delete[] b;
            }
        }
    }
    else
    {
        std::cout<<"inversive matrix is not exist"<<std::endl;
    }
    return res;
}

template<class Type>
void copy(const matrix<Type> & from, Type **to)
{
    for(int i=0;i<from.n;i++)
    {
        for(int j=0;j<from.m;j++)
        {
            to[i][j]=from.mtx[i][j];
        }
    }
}

template<class type>
void matrix<type>::gauss()
{
    Gauss(mtx,n,m);
}

template<class type>
void matrix<type>::tridiagonal()
{
    Tridiagonal(mtx,n,m);
}


#endif
