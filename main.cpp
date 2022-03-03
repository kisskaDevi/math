#include <iostream>

#include "kisskaMath/complex.h"
#include "kisskaMath/matrix.h"
#include "kisskaMath/function.h"
#include "timer.h"
#include "saveData.h"

void example1();
void example2();
void example3();

int main()
{
    timer<std::chrono::milliseconds> Timer;

std::cout<<"example1 complex================================"<<std::endl;
    example1();
std::cout<<"example2 matrix================================="<<std::endl;
    example2();
std::cout<<"example3 function==============================="<<std::endl;
    example3();

    std::cout<<"time "<<Timer.getTime()<<" milliseconds"<<endl;
    return 0;
}

//примеры

void example1()
{
    cn::complex<float> z(1.0f,2.0f);                            //определение комплексного числа z = 1 + 2i
    std::cout<<"z = "<<z<<std::endl;                            //выведем его
std::cout<<std::endl;
    cn::complex<float> a = z;                                   //определим аналогично второе число a и присвоим ему значение z
    if(a == z){std::cout<<"a = z"<<std::endl;}                  //если действительные и комплексные части равны, то оператор == возвращает true
std::cout<<std::endl;
    a = 2*z;                                                    //пусть теперь a = 2z
    if(a != z){std::cout<<"a = "<<a<<std::endl;}                //выведем результат
std::cout<<std::endl;
    std::cout<<"a+z = "<<a+z<<std::endl;                        //примеры математических операций
    std::cout<<"a-z = "<<a-z<<std::endl;
    std::cout<<"a*z = "<<a*z<<std::endl;
    std::cout<<"a/z = "<<a/z<<std::endl;                        //если нужна проверка деления на 0, можно подключить макрос #zerocheck, базово она отключена для повышения скорости вычислений
std::cout<<std::endl;
    std::cout<<"Re(a)   = "<<a.Real()<<std::endl;               //базовые операции над комплексными числами
    std::cout<<"Im(a)   = "<<a.Imaginary()<<std::endl;
    std::cout<<"abs(a)  = "<<a.abs()<<std::endl;
    std::cout<<"arg(a)  = "<<a.arg()<<std::endl;
    std::cout<<"con(a)  = "<<a.con()<<std::endl;                //комплексное сопряжение вызванное таким образом изменяет комплексное число, чтобы этого избежать можно совершить оту операцию через функцию
std::cout<<std::endl;
    std::cout<<"Re(z)   = "<<cn::Re(z)<<std::endl;              //то же самое через функции
    std::cout<<"Im(z)   = "<<cn::Im(z)<<std::endl;
    std::cout<<"abs(z)  = "<<cn::abs(z)<<std::endl;
    std::cout<<"arg(z)  = "<<cn::arg(z)<<std::endl;
    std::cout<<"con(z)  = "<<cn::con(z)<<std::endl;
std::cout<<std::endl;
    cn_out_setw = 4;
    cn_out_setprecision = 8;
    std::cout<<"z = "<<z<<std::endl;
    std::cout<<"a = "<<a<<std::endl;
std::cout<<std::endl;
    std::cout<<"exp(z) = "<<cn::exp(z)<<std::endl;              //примеры некоторых математических функций
    std::cout<<"cos^2(a) + sin^2(a) = "
    <<cn::cos(a)*cn::cos(a)+cn::sin(a)*cn::sin(a)<<std::endl;   //основное тригонометрическое тождество
    std::cout<<"pow(z,0.5) = "<<cn::pow(z,0.5f,0)<<std::endl;   //при возведении в степень так же необходимо указать кратность корня, поскольку у комплексного числа существует несколько корней, в данном случае k=0
    std::cout<<"tgh(Arctgh(z)) = "
            <<cn::tgh(cn::Arctgh(z,0))<<std::endl;
std::cout<<std::endl;
}

void example2()
{
    matrix<float> M = {{
        {2.5f, 0.2f},
        {5.0f,-5.0f} }};                                        //определим матрицу например так
    std::cout<<"M = \n"<<M<<std::endl;

    matrix<float> A = {{
        {1.5f, 2.2f, 1.0f, -5.0f},
        {9.0f, 0.0f, 3.0f, 11.0f},
        {7.5f, 0.7f, 8.0f,  5.0f} }};                           //или побольше
    std::cout<<"A = \n"<<A<<std::endl;
std::cout<<std::endl;
    M.addColumn(1,1);                                           //добавим в позицию 1 одну сроку в матрице M
    M(0,1) = 1.0; M(1,1) = -1.0;                                //и заполним пустые метса, по умолчания они равны 0
    std::cout<<"M = \n"<<M<<std::endl;
std::cout<<std::endl;
    std::cout<<"M*A = \n"<<M*A<<std::endl;                      //теперь мы их можем перемножить, но только в таком порядке, поскольку умножение матриц операция не коомутативная
std::cout<<std::endl;
    A.deleteColumn(3);                                          //удалим последнюю строку
    std::cout<<"A = \n"<<A<<std::endl;                          //теперь матрица квадратная

    matrix<float> B = replaceString(replaceColumn(A,0,2),0,2);  //добавим новую матрицу, равную матрице A с перевёрнутыми последними строками и столбцами
    std::cout<<"B = \n"<<B<<std::endl;
std::cout<<std::endl;
    std::cout<<"A+B = \n"<<A+B<<std::endl;                      //проверим все математические операции
    std::cout<<"A-B = \n"<<A-B<<std::endl;
    std::cout<<"A*B = \n"<<A*B<<std::endl;
std::cout<<std::endl;
    std::cout<<"tr(A) = "<<tr(A)<<std::endl;                    //поддерживаемые матричные операции
    std::cout<<"det(A) = "<<det(A)<<std::endl<<std::endl;
    std::cout<<"A^-1 = \n"<<T(A)<<std::endl;
    std::cout<<"A^T = \n"<<T(A)<<std::endl;
std::cout<<std::endl;
    cn::complex<float> i(0.0f,1.0f);
    cn_out_setw = 1;
    cn_out_setprecision = 5;

    matrix<cn::complex<float>> N = {{
        {2.0f+1.0f*i, 1.0f+2.0f*i},
        {5.0f+3.0f*i,-5.0f-4.0f*i} }};                          //можно определить матрицу и от комплексного числа
    std::cout<<"N = \n"<<N<<std::endl;

    std::cout<<"det = "<<det(N)<<std::endl;                     //можем вычислить, например, определитель комплексной матрицы
std::cout<<std::endl;
}

void example3()
{
    float pi = 4*std::atan(1.0f);

    X::function<float> x(-5.0f,5.0f,51);                        //задём область определения функции [-5,5], в данном случае сетка равномерная, но можно задать неравномерную через массив
    X::function<float> y = X::exp(-0.25f*x*x)*cos(pi*x);
    y.setout(4,6);

    std::cout<<y<<std::endl;                                    //вывод
std::cout<<std::endl;
    std::cout<<"integral(y) = "<<y.integralS2(0,51)<<std::endl; //интегрирование
std::cout<<std::endl;
    y=X::pow(y,2.0f);                                           //возведение в квадрат
    y.normalize();                                              //нормировка
    std::cout<<y<<std::endl;                                    //теперь у нас есть некоторая плотность вероятности
std::cout<<std::endl;
    std::cout<<"integral(y) = "<<y.integralS2(0,51)<<std::endl; //при повторном интегрировании результат будет равен 1
    std::cout<<"D(y) = "<<y.dispersion()<<std::endl;            //можем найти дисперсию
std::cout<<std::endl;
}
