#ifndef _COMPLEX_H_
#define _COMPLEX_H_

#include <iostream>
#include <cmath>
#include <iomanip>

int cn_out_setprecision=1;
int cn_out_setw=2;

namespace cn
{
    template <typename type>
    class complex;

    template<class T1,class T2> complex<T2> operator+(const T1 & C, const complex<T2> & other);   //такое определение позволяет умножать любой тип T на тип complex
    template<class T1,class T2> complex<T2> operator*(const T1 & C, const complex<T2> & other);   //
    template<class T1,class T2> complex<T2> operator-(const T1 & C, const complex<T2> & other);   //
    template<class T1,class T2> complex<T2> operator/(const T1 & C, const complex<T2> & other);   //

    template<class T> std::ostream & operator<< (std::ostream & out, const complex<T> & point);  //оператор вывода
    template<class T> std::istream & operator>> (std::istream & in, complex<T> & point);         //оператор ввода

    template<class T> T Re(const complex<T> & z);                                    //получить действительную часть в виде переменной типа T
    template<class T> T Im(const complex<T> & z);                                    //получить мнимую часть в виде переменной типа T
    template<class T> T abs(const complex<T> & z);                                   //получить модуль числа в виде переменной типа T
    template<class T> T arg(const complex<T> & z);                                   //получить аругмет числа в виде переменной типа T
    template<class T> complex<T> con(const complex<T> & z);                          //комплексное сопряжение

    template<class T> complex<T> exp(const complex<T> & z);                          //определение математических функций//
    template<class T> complex<T> Ln(const complex<T> & z, int k);                    //
    template<class T> complex<T> pow(const complex<T> & z, T deg, int k);            //
    template<class T> complex<T> sin(const complex<T> & z);                          //
    template<class T> complex<T> cos(const complex<T> & z);                          //
    template<class T> complex<T> tg(const complex<T> & z);                           //
    template<class T> complex<T> ctg(const complex<T> & z);                          //
    template<class T> complex<T> sinh(const complex<T> & z);                         //
    template<class T> complex<T> cosh(const complex<T> & z);                         //
    template<class T> complex<T> tgh(const complex<T> & z);                          //
    template<class T> complex<T> ctgh(const complex<T> & z);                         //
    template<class T> complex<T> Arcsin(const complex<T> & z, int k, int root);      //
    template<class T> complex<T> Arccos(const complex<T> & z, int k, int root);      //
    template<class T> complex<T> Arctg(const complex<T> & z, int k);                 //
    template<class T> complex<T> Arcctg(const complex<T> & z, int k);                //
    template<class T> complex<T> Arcsinh(const complex<T> & z, int k, int root);     //
    template<class T> complex<T> Arccosh(const complex<T> & z, int k, int root);     //
    template<class T> complex<T> Arctgh(const complex<T> & z, int k);                //
    template<class T> complex<T> Arcctgh(const complex<T> & z, int k);               //______определение закончено______//

    template <typename type>
    class complex
    {
        private:
            type Re;                                    //параметр числа - действительная часть
            type Im;                                    //параметр числа - мнимая часть
            int setprecision=cn_out_setprecision;       //параметр вывода - точность числа после запятой
            int setw=cn_out_setw;                       //параметр вывода - выделяемое под одно число пространство в консоли
        public:
            complex();                                  //Re=0; Im=0;
            complex(type Re);                           //конструктор задающий толкьо действительную часть
            complex(type Re, type Im);                  //конструктор
            complex(const complex<type>&other);         //конструктор копирования
            ~complex();                                 //деструктор

            bool        operator == (const complex<type> & other) const;                                         //_____переопределение операторов_____//
            bool        operator != (const complex<type> & other) const;                                         //
            bool        operator != (type other) const;                                                          //
            complex &   operator =  (const complex<type> & other);                                               //
            complex &   operator =  (const type & other);                                                        //
            complex     operator +  (const complex<type> & other) const;                                         //
            complex     operator -  (const complex<type> & other) const;                                         //
            complex     operator -  () const;                                                                    //
            complex     operator *  (const complex<type> & other) const;                                         //
            complex     operator /  (const complex<type> & other) const;                                         //
            complex     operator *  (const type & C) const;                                                      //
            complex     operator /  (const type & C) const;                                                      //
            complex &   operator += (const complex<type> & other);                                               //
            complex &   operator *= (const complex<type> & other);                                               //
            complex &   operator -= (const complex<type> & other);                                               //
            complex &   operator /= (const complex<type> & other);                                               //
            complex &   operator ++ ();                                                                          //
            complex &   operator ++ (int value);                                                                 //
            complex &   operator -- ();                                                                          //
            complex &   operator -- (int value);                                                                 //
            template<class T1,class T2> friend complex<T2> operator+(const T1 & C, const complex<T2> & other);   //такое определение позволяет умножать любой тип T на тип complex
            template<class T1,class T2> friend complex<T2> operator*(const T1 & C, const complex<T2> & other);   //
            template<class T1,class T2> friend complex<T2> operator-(const T1 & C, const complex<T2> & other);   //
            template<class T1,class T2> friend complex<T2> operator/(const T1 & C, const complex<T2> & other);   //

            void set(type Real, type Imm);                            //задать параметры числа Re & Im
            void setout(int setprecision, int setw);                  //задать параметры вывода

            type Real() const;                  //получить действительную часть
            type Imaginary() const;             //получить мнимую часть
            type abs() const;                   //получить модуль числа
            type arg() const;                   //получить аругмет числа
            complex<type> & con();              //комплексное сопряжение

            template<class T> friend std::ostream & operator<< (std::ostream & out, const complex<T> & point);  //оператор вывода
            template<class T> friend std::istream & operator>> (std::istream & in, complex<T> & point);         //оператор ввода

            template<class T> friend T Re(const complex<T> & z);                                    //получить действительную часть в виде переменной типа T
            template<class T> friend T Im(const complex<T> & z);                                    //получить мнимую часть в виде переменной типа T
            template<class T> friend T abs(const complex<T> & z);                                   //получить модуль числа в виде переменной типа T
            template<class T> friend T arg(const complex<T> & z);                                   //получить аругмет числа в виде переменной типа T
            template<class T> friend complex<T> con(const complex<T> & z);                          //комплексное сопряжение

            template<class T> friend complex<T> exp(const complex<T> & z);                          //определение математических функций//
            template<class T> friend complex<T> Ln(const complex<T> & z, int k);                    //
            template<class T> friend complex<T> pow(const complex<T> & z, T deg, int k);            //
            template<class T> friend complex<T> sin(const complex<T> & z);                          //
            template<class T> friend complex<T> cos(const complex<T> & z);                          //
            template<class T> friend complex<T> tg(const complex<T> & z);                           //
            template<class T> friend complex<T> ctg(const complex<T> & z);                          //
            template<class T> friend complex<T> sinh(const complex<T> & z);                         //
            template<class T> friend complex<T> cosh(const complex<T> & z);                         //
            template<class T> friend complex<T> tgh(const complex<T> & z);                          //
            template<class T> friend complex<T> ctgh(const complex<T> & z);                         //
            template<class T> friend complex<T> Arcsin(const complex<T> & z, int k, int root);      //
            template<class T> friend complex<T> Arccos(const complex<T> & z, int k, int root);      //
            template<class T> friend complex<T> Arctg(const complex<T> & z, int k);                 //
            template<class T> friend complex<T> Arcctg(const complex<T> & z, int k);                //
            template<class T> friend complex<T> Arcsinh(const complex<T> & z, int k, int root);     //
            template<class T> friend complex<T> Arccosh(const complex<T> & z, int k, int root);     //
            template<class T> friend complex<T> Arctgh(const complex<T> & z, int k);                //
            template<class T> friend complex<T> Arcctgh(const complex<T> & z, int k);               //

            complex                             & exp();
    };

    //

    template <typename type>
    class function: public complex<type>
    {
        private:
            type **x;                       //действительная часть аргумена
            type **y;                       //мнимая часть аргумента
            int nx;
            int ny;
            cn::complex<type> **z;          //функция от аргумента
            int setprecision=1, setw=2;     //параметны вывода
        public:
            function(type **x, type **y, cn::complex<type> **z, int nx, int ny);        //конструктор задаёт как аргумент, так и функцию
            function(type **x, type **y, int nx, int ny);                               //конструктор задаёт аргумент, функция задаётся автоматически в виде z=x+iy, удобно использовать как аргумент функций
            function(const function<type> & other);                                     //конструктор копирования
            ~function();                                                                //деструктор

            bool operator == (const function<type> & other);                            //операторы
            bool operator != (const function<type> & other);                            //
            function operator + (const function<type> & other);                         //
            function operator - (const function<type> & other);                         //
            function operator - ();                                                     //
            function operator * (const function<type> & other);                         //
            function operator / (const function<type> & other);                         //
            function operator * (const type & C);                                       //
            function operator / (const type & C);                                       //
            function &operator = (const function<type> & other);                        //
            complex<type> & operator () (const int i, const int j);                     //доступ к элементу в точке (i;j)

            void setf(int i, int j, const complex<type> & other);                       //задать элемент в точке (i;j)
            void setx(int i, int j, const type & other);
            void sety(int i, int j, const type & other);
            int getnx() const;
            int getny() const;
            type getx(int i, int j) const;
            type gety(int i, int j) const;
            complex<type> getf(int i, int j) const;

            type Re(int i, int j) const;                                      //действительное значение функции в точке (i;j)
            type Im(int i, int j) const;                                      //мнимое значение функции в точке (i;j)
            type abs(int i, int j) const;                                     //модуль функции в точке (i;j)
            type arg(int i, int j) const;                                     //аргумент функции в точке (i;j)
            void setout(int setprecision, int setw);                          //задать параметры вывода

            template<class T> friend std::ostream & operator<< (std::ostream & out, const function<T> & point);  //вывести функцию

            template<class T> friend function<T> con(const function<T> & other);                                 //комплексное сопряжение
            template<class T> friend function<T> exp(const function<T> & other);                                 //определение математических функций
            template<class T> friend function<T> Ln(const function<T> & other, int k);                           //
            template<class T> friend function<T> pow(const function<T> & other, T deg, int k);                   //
            template<class T> friend function<T> sin(const function<T> & other);                                 //
            template<class T> friend function<T> cos(const function<T> & other);                                 //
            template<class T> friend function<T> tg(const function<T> & other);                                  //
            template<class T> friend function<T> ctg(const function<T> & other);                                 //
            template<class T> friend function<T> sinh(const function<T> & other);                                //
            template<class T> friend function<T> cosh(const function<T> & other);                                //
            template<class T> friend function<T> tgh(const function<T> & other);                                 //
            template<class T> friend function<T> ctgh(const function<T> & other);                                //
            template<class T> friend function<T> Arcsin(const function<T> & other, int k, int root);             //
            template<class T> friend function<T> Arccos(const function<T> & other, int k, int root);             //
            template<class T> friend function<T> Arctg(const function<T> & other, int k);                        //
            template<class T> friend function<T> Arcctg(const function<T> & other, int k);                       //
            template<class T> friend function<T> Arcsinh(const function<T> & other, int k, int root);            //
            template<class T> friend function<T> Arccosh(const function<T> & other, int k, int root);            //
            template<class T> friend function<T> Arctgh(const function<T> & other, int k);                       //
            template<class T> friend function<T> Arcctgh(const function<T> & other, int k);                      //
    };

    template <>
    class complex<char>{};
    template <>
    class function<char>{};
}

//============================================================================================================================================//
//============================================================================================================================================//
//============================================================================================================================================//

namespace cn {

    template <typename type>
    cn::complex<type>::complex()
    {
        this->Re=0;
        this->Im=0;
    }

    template <typename type>
    cn::complex<type>::complex(type Re)
    {
        this->Re=Re;
        this->Im=0;
    }

    template <typename type>
    cn::complex<type>::complex(type Re, type Im)
    {
        this->Re=Re;
        this->Im=Im;
    }

    template <typename type>
    cn::complex<type>::~complex(){}

    template <typename type>
    cn::complex<type>::complex(const complex<type>&other)
    {
        this->Re=other.Re;
        this->Im=other.Im;
        this->setw=other.setw;
        this->setprecision=other.setprecision;
    }

    template <typename type>
    cn::complex<type> & cn::complex<type>::operator = (const complex<type>&other)
    {
        this->Re=other.Re;
        this->Im=other.Im;
        this->setw=other.setw;
        this->setprecision=other.setprecision;
        return *this;
    }

    template <typename type>
    cn::complex<type> & cn::complex<type>::operator = (const type& other)
    {
        this->Re=other;
        this->Im=0;
        return *this;
    }

    template <typename type>
    bool cn::complex<type>::operator == (const complex<type> & other) const
    {
        return this->Re==other.Re && this->Im==other.Im;
    }

    template <typename type>
    bool cn::complex<type>:: operator != (const complex<type> & other) const
    {
        return !(this->Re==other.Re && this->Im==other.Im);
    }

    template <typename type>
    bool cn::complex<type>:: operator != (type other) const
    {
        return !(std::sqrt(Re*Re+Im*Im)==other);
    }

    template <typename type>
    cn::complex<type> cn::complex<type>:: operator + (const complex<type>&other) const
    {
        complex<type> res(
        this->Re + other.Re,
        this->Im + other.Im);
        return res;
    }

    template <typename type>
    cn::complex<type> cn::complex<type>:: operator - (const complex<type>&other) const
    {
        complex<type> res(
        this->Re - other.Re,
        this->Im - other.Im);
        return res;
    }

    template <typename type>
    cn::complex<type> cn::complex<type>:: operator - () const
    {
        complex<type> res(
        this->Re*(-1),
        this->Im*(-1));
        return res;
    }

    template <typename type>
    cn::complex<type> cn::complex<type>:: operator * (const complex<type>&other) const
    {
        complex<type> res(
        this->Re*other.Re - this->Im*other.Im,
        this->Re*other.Im + this->Im*other.Re);
        return res;
    }

    template <typename type>
    cn::complex<type> cn::complex<type>:: operator / (const complex<type>&other) const
    {
        #ifdef zerocheck
            if(other.Re*other.Re+other.Im*other.Im!=0)
            {
                complex<type> res(*this);
                res.Re = (this->Re*other.Re + this->Im*other.Im)/(other.Re*other.Re+other.Im*other.Im);
                res.Im = (this->Im*other.Re - this->Re*other.Im)/(other.Re*other.Re+other.Im*other.Im);
                return res;
            }
            else
            {
                std::cout<<"undefined operation in cn::complex<type> cn::complex<type>:: operator / (const complex<type>&other), returned 0"<<std::endl;
                return 0;
            }
        #else
            complex<type> res(
            (this->Re*other.Re + this->Im*other.Im)/(other.Re*other.Re+other.Im*other.Im),
            (this->Im*other.Re - this->Re*other.Im)/(other.Re*other.Re+other.Im*other.Im));
            return res;
        #endif
    }

    template <typename type>
    cn::complex<type> cn::complex<type>:: operator * (const type & C) const
    {
        complex<type> res(this->Re*C, this->Im*C);
        return res;
    }

    template <typename type>
    cn::complex<type> cn::complex<type>:: operator / (const type & C) const
    {
        #ifdef zerocheck
            if(C!=0)
            {
                complex<type> res(*this);
                res.Re = (this->Re)/C;
                res.Im = (this->Im)/C;
                return res;
            }
            else
            {
                std::cout<<"undefined operation in cn::complex<type> cn::complex<type>:: operator / (const type & C), returned 0"<<std::endl;
                return 0;
            }
        #else
            complex<type> res((this->Re)/C,(this->Im)/C);
            return res;
        #endif
    }

    template<typename type>
    complex<type> & complex<type>::operator +=(const complex<type> &other)
    {
        this->Re += other.Re;
        this->Im += other.Im;
        return *this;
    }


    template<typename type>
    complex<type> & complex<type>::operator *=(const complex<type> &other)
    {
        type re = Re;
        type im = Im;
        this->Re = re*other.Re - im*other.Im;
        this->Im = re*other.Im + im*other.Re;
        return *this;
    }

    template<typename type>
    complex<type> & complex<type>::operator -=(const complex<type> &other)
    {
        this->Re -= other.Re;
        this->Im -= other.Im;
        return *this;
    }

    template<typename type>
    complex<type> & complex<type>::operator /=(const complex<type> &other)
    {
        type re = Re;
        type im = Im;
        this->Re = (re*other.Re + im*other.Im)/(other.Re*other.Re+other.Im*other.Im);
        this->Im = (im*other.Re - re*other.Im)/(other.Re*other.Re+other.Im*other.Im);
        return *this;
    }

    template <typename type>
    cn::complex<type> & cn::complex<type>:: operator ++()
    {
        this->Re++;
        this->Im++;
        return *this;
    }

    template <typename type>
    cn::complex<type> & cn::complex<type>:: operator ++(int value)
    {
        (void)value;
        this->Re++;
        this->Im++;
        return *this;
    }

    template <typename type>
    cn::complex<type> & cn::complex<type>:: operator --()
    {
        this->Re--;
        this->Im--;
        return *this;
    }

    template <typename type>
    cn::complex<type> & cn::complex<type>:: operator --(int value)
    {
        (void)value;
        this->Re--;
        this->Im--;
        return *this;
    }

    template<class T1,class T2>
    complex<T2> operator*(const T1 & C, const complex<T2> & other)
    {
        complex<T2> res(other.Re*C,other.Im*C);
        res.setw = other.setw;
        res.setprecision = other.setprecision;
        return res;
    }

    template<class T1,class T2>
    complex<T2> operator+(const T1 & C, const complex<T2> & other)
    {
        complex<T2> res(C + other.Re,other.Im);
        res.setw = other.setw;
        res.setprecision = other.setprecision;
        return res;
    }

    template<class T1,class T2>
    complex<T2> operator-(const T1 & C, const complex<T2> & other)
    {
        complex<T2> res(C - other.Re,- other.Im);
        res.setw = other.setw;
        res.setprecision = other.setprecision;
        return res;
    }

    template<class T1,class T2>
    complex<T2> operator/(const T1 & C, const complex<T2> & other)
    {
        complex<T2> res(
        (C*other.Re)/(other.Re*other.Re+other.Im*other.Im),
        (-C*other.Im)/(other.Re*other.Re+other.Im*other.Im));
        res.setw = other.setw;
        res.setprecision = other.setprecision;
        return res;
    }

    template<typename type>
    void complex<type>::setout(int setprecision, int setw)
    {
        this->setprecision=setprecision;
        this->setw=setw;
    }

    template <typename type>
    void cn::complex<type>::set(type Real, type Imm)
    {
        Re=Real;
        Im=Imm;
    }

    template <typename type>
    type cn::complex<type>::Real() const
    {
        return Re;
    }

    template <typename type>
    type cn::complex<type>::Imaginary() const
    {
        return Im;
    }

    template <typename type>
    type cn::complex<type>::abs() const
    {
        return sqrt(Re*Re+Im*Im);
    }

    template <typename type>
    type cn::complex<type>::arg() const
    {
        type pi=4*std::atan(1);
        type arg;
        if(Re>0)
        {
            arg=std::atan(Im/Re);
        }else if(Re<0)
        {
            if(Im>=0)
            {
                arg=pi+std::atan(Im/Re);
            }
            else
            {
                arg=-pi+std::atan(Im/Re);
            }
        }else if(Re==0)
        {
            if(Im>0)
            {
                arg=0.5*pi;
            }else if(Im<0)
            {
                arg=-0.5*pi;
            }else if(Im==0)
            {
                std::cout<<"null vector, returned 0"<<std::endl;
                arg=0;
            }
        }
        return arg;
    }

    template <typename type>
    cn::complex<type> & cn::complex<type>::con()
    {
        Im=-Im;
        return *this;
    }

    template<class T>
    std::ostream & operator<<(std::ostream &out, const cn::complex<T> &point)
    {
        cn::complex<T> z0=point;
        if(std::abs(point.Im)<1e-8){z0.Im=0;}
        if(std::abs(point.Re)<1e-8){z0.Re=0;}
        if(z0.Im<0)
            out <<std::setw(point.setw)<<std::fixed<<std::setprecision(point.setprecision)<<z0.Re<<"-"<<std::setw(point.setw)<<std::abs(z0.Im)<<"i";
        else
            out <<std::setw(point.setw)<<std::fixed<<std::setprecision(point.setprecision)<<z0.Re<<"+"<<std::setw(point.setw)<<z0.Im<<"i";
        return out;
    }

    template<class T>
    std::istream & operator>>(std::istream &in, cn::complex<T> &point)
    {
        std::cout<<"Re = ";
        in >> point.Re;
        std::cout<<"Im = ";
        in >> point.Im;
        return in;
    }

    template<class T>
    cn::complex<T> con(const complex<T> & z)
    {
        cn::complex<T> z0(z);
        z0.Re=z.Re;
        z0.Im=-z.Im;
        return z0;
    }

    template<class T>
    T Re(const complex<T> & z)
    {
        return z.Re;
    }

    template<class T>
    T Im(const complex<T> & z)
    {
        return z.Im;
    }

    template<class T>
    T abs(const complex<T> & z)
    {
        return sqrt(z.Re*z.Re+z.Im*z.Im);
    }

    template<class T>
    T arg(const complex<T> & z)
    {
        T pi=4*std::atan(1);
        T arg;
        if(z.Re>0)
        {
            arg=std::atan(z.Im/z.Re);
        }else if(z.Re<0)
        {
            if(z.Im>=0)
            {
                arg=pi+std::atan(z.Im/z.Re);
            }
            else
            {
                arg=-pi+std::atan(z.Im/z.Re);
            }
        }else if(z.Re==0)
        {
            if(z.Im>0)
            {
                arg=0.5*pi;
            }else if(z.Im<0)
            {
                arg=-0.5*pi;
            }else if(z.Im==0)
            {
                std::cout<<"null vector, returned 0"<<std::endl;
                arg=0;
            }
        }
        return arg;
    }

    template<class T>
    complex<T> exp(const complex<T>& z)
    {
        cn::complex<T> z0(z);
        z0.Re=std::exp(z.Re)*std::cos(z.Im);
        z0.Im=std::exp(z.Re)*std::sin(z.Im);
        return z0;
    }

    template<typename type>
    complex<type> &complex<type>::exp()
    {
        type re = Re;
        type im = Im;
        this->Re=std::exp(re)*std::cos(im);
        this->Im=std::exp(re)*std::sin(im);
        return *this;
    }

    template<class T>
    complex<T> Ln(const complex<T> & z, int k)
    {
        T pi=4*atan(1);
        cn::complex<T> z0(z);
        z0.Re=std::log(abs(z));
        z0.Im=arg(z)+2*pi*k;
        return z0;
    }

    template<class T>
    complex<T> pow(const complex<T> & z, T deg, int k)
    {
        T pi=4*atan(1);
        cn::complex<T> z0(z);
        z0.Re=std::pow(abs(z),deg)*std::cos(deg*(arg(z)+2*pi*k));
        z0.Im=std::pow(abs(z),deg)*std::sin(deg*(arg(z)+2*pi*k));
        return z0;
    }

    template<class T>
    complex<T> sin(const complex<T> & z)
    {
        cn::complex<T> z0(z);
        z0.Re=std::sin(z.Re)*std::cosh(z.Im);
        z0.Im=std::cos(z.Re)*std::sinh(z.Im);
        return z0;
    }

    template<class T>
    complex<T> cos(const complex<T> & z)
    {
        cn::complex<T> z0(z);
        z0.Re=std::cos(z.Re)*std::cosh(z.Im);
        z0.Im=-std::sin(z.Re)*std::sinh(z.Im);
        return z0;
    }

    template<class T>
    complex<T> tg(const complex<T> & z)
    {
        cn::complex<T> z0(z);
        z0=cn::sin(z)/cn::cos(z);
        return z0;
    }

    template<class T>
    complex<T> ctg(const complex<T> & z)
    {
        cn::complex<T> z0(z);
        z0=cn::cos(z)/cn::sin(z);
        return z0;
    }

    template<class T>
    complex<T> sinh(const complex<T> & z)
    {
        cn::complex<T> z0(z);
        z0=cn::exp(z)-cn::exp(-z);
        return z0/2;
    }

    template<class T>
    complex<T> cosh(const complex<T> & z)
    {
        cn::complex<T> z0(z);
        z0=cn::exp(z)+cn::exp(-z);
        return z0/2;
    }

    template<class T>
    complex<T> tgh(const complex<T> & z)
    {
        cn::complex<T> z0(z);
        z0=cn::sinh(z)/cn::cosh(z);
        return z0;
    }

    template<class T>
    complex<T> ctgh(const complex<T> & z)
    {
        cn::complex<T> z0(z);
        z0=cn::cosh(z)/cn::sinh(z);
        return z0;
    }

    template<class T>
    complex<T> Arcsin(const complex<T> & z, int k, int root)
    {
        cn::complex<T> z0(z);
        cn::complex<T> i(0,1);
        cn::complex<T> one(1,0);
        T deg=std::sin(2*std::atan(1)/3);
        return z0=-i*cn::Ln(i*z+cn::pow(one-z*z,deg,root),k);
    }

    template<class T>
    complex<T> Arccos(const complex<T> & z, int k, int root)
    {
        cn::complex<T> z0(z);
        cn::complex<T> i(0,1);
        cn::complex<T> one(1,0);
        cn::complex<T> pi(4*atan(1),0);
        T deg=std::sin(2*std::atan(1)/3);
        return z0=pi/2+i*cn::Ln(i*z+cn::pow(one-z*z,deg,root),k);
    }

    template<class T>
    complex<T> Arctg(const complex<T> & z, int k)
    {
        cn::complex<T> z0(z);
        cn::complex<T> i(0,1);
        cn::complex<T> one(1,0);
        return z0=-i*cn::Ln((one+i*z)/(one-i*z),k)/2;
    }

    template<class T>
    complex<T> Arcctg(const complex<T> & z, int k)
    {
        cn::complex<T> z0(z);
        cn::complex<T> i(0,1);
        cn::complex<T> one(1,0);
        return z0=-i*cn::Ln((i*z-one)/(i*z+one),k)/2;
    }

    template<class T>
    complex<T> Arcsinh(const complex<T> & z, int k, int root)
    {
        cn::complex<T> z0(z);
        cn::complex<T> i(0,1);
        cn::complex<T> one(1,0);
        T deg=std::sin(2*std::atan(1)/3);
        return z0=cn::Ln(z+cn::pow(one+z*z,deg,root),k);
    }

    template<class T>
    complex<T> Arccosh(const complex<T> & z, int k, int root)
    {
        cn::complex<T> z0(z);
        cn::complex<T> i(0,1);
        cn::complex<T> one(1,0);
        T deg=std::sin(2*std::atan(1)/3);
        return z0=cn::Ln(z+cn::pow(z*z-one,deg,root),k);
    }

    template<class T>
    complex<T> Arctgh(const complex<T> & z, int k)
    {
        cn::complex<T> z0(z);
        cn::complex<T> i(0,1);
        cn::complex<T> one(1,0);
        return z0=cn::Ln((one+z)/(one-z),k)/2;
    }

    template<class T>
    complex<T> Arcctgh(const complex<T> & z, int k)
    {
        cn::complex<T> z0(z);
        cn::complex<T> i(0,1);
        cn::complex<T> one(1,0);
        return z0=cn::Ln((z+one)/(z-one),k)/2;
    }

    //=============================================================================//
    //==============================COMPLEX=FUNCTION==============================//
    //=============================================================================//

    template<typename type>
    function<type> :: function (type **x, type **y, cn::complex<type> **z, int nx, int ny)
    {
        this->nx = nx;
        this->ny = ny;
        this->x = new type*[nx];
        this->y = new type*[nx];
        this->z = new cn::complex<type>*[nx];
        for(int i=0;i<nx;i++)
        {
            this->x[i]= new type[ny];
            this->y[i]= new type[ny];
            this->z[i]= new cn::complex<type>[ny];
            for(int j=0;j<ny;j++)
            {
                this->x[i][j]=x[i][j];
                this->y[i][j]=y[i][j];
                this->z[i][j]=z[i][j];
            }
        }
    }

    template<typename type>
    function<type> :: function (type **x, type **y, int nx, int ny)
    {
        this->nx = nx;
        this->ny = ny;
        this->x = new type*[nx];
        this->y = new type*[nx];
        this->z = new cn::complex<type>*[nx];
        for(int i=0;i<nx;i++)
        {
            this->x[i]= new type[ny];
            this->y[i]= new type[ny];
            this->z[i]= new cn::complex<type>[ny];
            for(int j=0;j<ny;j++)
            {
                this->x[i][j]=x[i][j];
                this->y[i][j]=y[i][j];
                this->z[i][j]=z[i][j].set(x[i][j],y[i][j]);
            }
        }
    }

    template<typename type>
    cn::function<type>::function(const cn::function<type> & other) : complex<type>(other)
    {
        this->setprecision=other.setprecision;
        this->setw=other.setw;
        this->nx = other.nx;
        this->ny = other.ny;
        this->x = new type*[nx];
        this->y = new type*[nx];
        this->z = new cn::complex<type>*[nx];
        for(int i=0;i<nx;i++)
        {
            this->x[i]= new type[ny];
            this->y[i]= new type[ny];
            this->z[i]= new cn::complex<type>[ny];
            for(int j=0;j<ny;j++)
            {
                this->x[i][j]=other.x[i][j];
                this->y[i][j]=other.y[i][j];
                this->z[i][j]=other.z[i][j];
                this->z[i][j].setout(setprecision,setw);
            }
        }
    }

    template<typename type>
    function<type>::~function()
    {
        for(int i=0;i<nx;i++)
        {
            delete [] this->x[i];
            delete [] this->y[i];
            delete [] this->z[i];
        }
        delete [] this->x;
        delete [] this->y;
        delete [] this->z;
    }

    template<typename type>
    bool function<type>::operator ==(const function<type> &other)
    {
        for(int i=0; i<nx;i++)
        {
            for(int j=0; j<ny;j++)
            {
                if(
                    (this->z[i][j]==other.z[i][j])&&
                    (this->x[i][j]==other.x[i][j])&&
                    (this->y[i][j]==other.y[i][j])
                        ){}
                else{return false;}
            }
        }
        return true;
    }

    template<typename type>
    bool function<type>::operator !=(const function<type> &other)
    {
        for(int i=0; i<nx;i++)
        {
            for(int j=0; j<ny;j++)
            {
                if(
                    (this->z[i][j]==other.z[i][j])&&
                    (this->x[i][j]==other.x[i][j])&&
                    (this->y[i][j]==other.y[i][j])
                        ){}
                else{return true;}
            }
        }
        return false;
    }

    template<typename type>
    function<type> function<type>::operator +(const function<type> &other)
    {
        function<type> res(*this);
        for(int i=0; i<nx;i++)
        {
            for(int j=0; j<ny;j++)
            {
                res.z[i][j]=this->z[i][j]+other.z[i][j];
            }
        }
        return res;
    }

    template<typename type>
    function<type> function<type>::operator -(const function<type> &other)
    {
        function<type> res(*this);
        for(int i=0; i<nx;i++)
        {
            for(int j=0; j<ny;j++)
            {
                res.z[i][j]=this->z[i][j]-other.z[i][j];
            }
        }
        return res;
    }

    template<typename type>
    function<type> function<type>::operator -()
    {
        function<type> res(*this);
        for(int i=0; i<nx;i++)
        {
            for(int j=0; j<ny;j++)
            {
                res.z[i][j]=this->z[i][j]*(-1);
            }
        }
        return res;
    }

    template<typename type>
    function<type> function<type>::operator * (const function<type> &other)
    {
        function<type> res(*this);
        for(int i=0; i<nx;i++)
        {
            for(int j=0; j<ny;j++)
            {
                res.z[i][j]=this->z[i][j]*other.z[i][j];
            }
        }
        return res;
    }

    template<typename type>
    function<type> function<type>::operator / (const function<type> &other)
    {
        function<type> res(*this);
        for(int i=0; i<nx;i++)
        {
            for(int j=0; j<ny;j++)
            {
                res.z[i][j]=this->z[i][j]/other.z[i][j];
            }
        }
        return res;
    }

    template<typename type>
    function<type> function<type>::operator * (const type & C)
    {
        function<type> res(*this);
        for(int i=0; i<nx;i++)
        {
            for(int j=0; j<ny;j++)
            {
                res.z[i][j]=this->z[i][j]*C;
            }
        }
        return res;
    }

    template<typename type>
    function<type> function<type>::operator / (const type & C)
    {
        function<type> res(*this);
        if(C!=0)
        {
            for(int i=0; i<nx;i++)
            {
                for(int j=0; j<ny;j++)
                {
                    res.z[i][j]=this->z[i][j]/C;
                }
            }
        }
        else
        {
            std::cout<<"undefined operation in cn::function<type> cn::function<type>::operator / (const type & C), returned 0"<<std::endl;
            for(int i=0; i<nx;i++)
            {
                for(int j=0; j<ny;j++)
                {
                    res.z[i][j]=0;
                }
            }
        }
        return res;
    }

    template<typename type>
    function<type> & function<type>::operator =(const function<type> &other)
    {
        setprecision=other.setprecision;
        setw=other.setw;
        if(this->x!=nullptr)
        {
            for(int i=0; i<nx;i++)
            {
                delete[] this->x[i];
            }
            delete[] this->x;
        }
        if(this->y!=nullptr)
        {
            for(int i=0; i<nx;i++)
            {
                delete[] this->y[i];
            }
            delete[] this->y;
        }
        if(this->z!=nullptr)
        {
            for(int i=0; i<nx;i++)
            {
                delete[] this->z[i];
            }
            delete[] this->z;
        }
        this->nx = other.nx;
        this->ny = other.ny;
        this->x = new type*[nx];
        this->y = new type*[nx];
        this->z = new cn::complex<type>*[nx];
        for(int i=0;i<nx;i++)
        {
            this->x[i]= new type[ny];
            this->y[i]= new type[ny];
            this->z[i]= new cn::complex<type>[ny];
            for(int j=0;j<ny;j++)
            {
                this->x[i][j]=other.x[i][j];
                this->y[i][j]=other.y[i][j];
                this->z[i][j]=other.z[i][j];
                this->z[i][j].setout(setprecision,setw);
            }
        }
        return *this;
    }

    template<typename type>
    complex<type> & function<type>::operator () (const int i, const int j)
    {
        return z[i][j];
    }

    template<typename type>
    void function<type>::setf(int i, int j, const complex<type> & other)
    {
        this->z[i][j]=other;
    }

    template<typename type>
    void function<type>::setx(int i, int j, const type & other)
    {
        this->x[i][j]=other;
    }

    template<typename type>
    void function<type>::sety(int i, int j, const type & other)
    {
        this->y[i][j]=other;
    }

    template<typename type>
    int function<type>::getnx() const
    {
        return nx;
    }

    template<typename type>
    int function<type>::getny() const
    {
        return ny;
    }

    template<typename type>
    type function<type>::getx(int i, int j) const
    {
        return x[i][j];
    }

    template<typename type>
    type function<type>::gety(int i, int j) const
    {
        return y[i][j];
    }

    template<typename type>
    complex<type> function<type>::getf(int i, int j) const
    {
        return z[i][j];
    }

    template<typename type>
    type function<type>::Re(int i, int j) const
    {
        return cn::Re(this->z[i][j]);
    }

    template<typename type>
    type function<type>::Im(int i, int j) const
    {
        return cn::Im(this->z[i][j]);
    }

    template<typename type>
    type function<type>::abs(int i, int j) const
    {
        return cn::abs(this->z[i][j]);
    }

    template<typename type>
    type function<type>::arg(int i, int j) const
    {
        return cn::arg(this->z[i][j]);
    }

    template<typename type>
    void function<type>::setout(int setprecision, int setw)
    {
        this->setprecision=setprecision;
        this->setw=setw;
        for(int i=0;i<nx;i++)
        {
            for(int j=0;j<ny;j++)
            {
                z[i][j].setout(setprecision,setw);
            }
        }
    }

    template<class type>
    std::ostream & operator<<(std::ostream &out, const function<type> &point)
    {
        for(int i=0;i<point.getnx();i++)
        {
            for(int j=0;j<point.getny();j++)
            {
                std::cout<<point.z[i][j]<<"  ";
            }
            std::cout<<std::endl;
        }
        return out;
    }

    template<class type>
    function<type> con(const function<type> & other)
    {
        function<type> res(other);
        for(int i=0;i<other.getnx();i++)
        {
            for(int j=0;j<other.getny();j++)
            {
                res.setf(i,j,con(other.getf(i,j)));
            }
        }
        return res;
    }

    template<class type>
    function<type> exp(const function<type> & other)
    {
        function<type> res(other);
        for(int i=0;i<other.getnx();i++)
        {
            for(int j=0;j<other.getny();j++)
            {
                res.setf(i,j,exp(other.getf(i,j)));
            }
        }
        return res;
    }

    template<class type>
    function<type> Ln(const function<type> & other, int k)
    {
        function<type> res(other);
        for(int i=0;i<other.getnx();i++)
        {
            for(int j=0;j<other.getny();j++)
            {
                res.setf(i,j,cn::Ln(other.getf(i,j),k));
            }
        }
        return res;
    }

    template<class type>
    function<type> pow(const function<type> & other, type deg, int k)
    {
        function<type> res(other);
        for(int i=0;i<other.getnx();i++)
        {
            for(int j=0;j<other.getny();j++)
            {
                res.setf(i,j,pow(other.getf(i,j),deg,k));
            }
        }
        return res;
    }

    template<class type>
    function<type> sin(const function<type> & other)
    {
        function<type> res(other);
        for(int i=0;i<other.getnx();i++)
        {
            for(int j=0;j<other.getny();j++)
            {
                res.setf(i,j,sin(other.getf(i,j)));
            }
        }
        return res;
    }

    template<class type>
    function<type> cos(const function<type> & other)
    {
        function<type> res(other);
        for(int i=0;i<other.getnx();i++)
        {
            for(int j=0;j<other.getny();j++)
            {
                res.setf(i,j,cos(other.getf(i,j)));
            }
        }
        return res;
    }

    template<class type>
    function<type> tg(const function<type> & other)
    {
        function<type> res(other);
        for(int i=0;i<other.getnx();i++)
        {
            for(int j=0;j<other.getny();j++)
            {
                res.setf(i,j,tg(other.getf(i,j)));
            }
        }
        return res;
    }

    template<class type>
    function<type> ctg(const function<type> & other)
    {
        function<type> res(other);
        for(int i=0;i<other.getnx();i++)
        {
            for(int j=0;j<other.getny();j++)
            {
                res.setf(i,j,ctg(other.getf(i,j)));
            }
        }
        return res;
    }

    template<class type>
    function<type> sinh(const function<type> & other)
    {
        function<type> res(other);
        for(int i=0;i<other.getnx();i++)
        {
            for(int j=0;j<other.getny();j++)
            {
                res.setf(i,j,sinh(other.getf(i,j)));
            }
        }
        return res;
    }

    template<class type>
    function<type> cosh(const function<type> & other)
    {
        function<type> res(other);
        for(int i=0;i<other.getnx();i++)
        {
            for(int j=0;j<other.getny();j++)
            {
                res.setf(i,j,cosh(other.getf(i,j)));
            }
        }
        return res;
    }

    template<class type>
    function<type> tgh(const function<type> & other)
    {
        function<type> res(other);
        for(int i=0;i<other.getnx();i++)
        {
            for(int j=0;j<other.getny();j++)
            {
                res.setf(i,j,tgh(other.getf(i,j)));
            }
        }
        return res;
    }

    template<class type>
    function<type> ctgh(const function<type> & other)
    {
        function<type> res(other);
        for(int i=0;i<other.getnx();i++)
        {
            for(int j=0;j<other.getny();j++)
            {
                res.setf(i,j,ctgh(other.getf(i,j)));
            }
        }
        return res;
    }

    template<class type>
    function<type> Arcsin(const function<type> & other, int k, int root)
    {
        function<type> res(other);
        for(int i=0;i<other.getnx();i++)
        {
            for(int j=0;j<other.getny();j++)
            {
                res.setf(i,j,Arcsin(other.getf(i,j),k,root));
            }
        }
        return res;
    }

    template<class type>
    function<type> Arccos(const function<type> & other, int k, int root)
    {
        function<type> res(other);
        for(int i=0;i<other.getnx();i++)
        {
            for(int j=0;j<other.getny();j++)
            {
                res.setf(i,j,Arccos(other.getf(i,j),k,root));
            }
        }
        return res;
    }

    template<class type>
    function<type> Arctg(const function<type> & other, int k)
    {
        function<type> res(other);
        for(int i=0;i<other.getnx();i++)
        {
            for(int j=0;j<other.getny();j++)
            {
                res.setf(i,j,Arctg(other.getf(i,j),k));
            }
        }
        return res;
    }

    template<class type>
    function<type> Arcctg(const function<type> & other, int k)
    {
        function<type> res(other);
        for(int i=0;i<other.getnx();i++)
        {
            for(int j=0;j<other.getny();j++)
            {
                res.setf(i,j,Arcctg(other.getf(i,j),k));
            }
        }
        return res;
    }

    template<class type>
    function<type> Arcsinh(const function<type> & other, int k, int root)
    {
        function<type> res(other);
        for(int i=0;i<other.getnx();i++)
        {
            for(int j=0;j<other.getny();j++)
            {
                res.setf(i,j,Arcsinh(other.getf(i,j),k,root));
            }
        }
        return res;
    }

    template<class type>
    function<type> Arccosh(const function<type> & other, int k, int root)
    {
        function<type> res(other);
        for(int i=0;i<other.getnx();i++)
        {
            for(int j=0;j<other.getny();j++)
            {
                res.setf(i,j,Arccosh(other.getf(i,j),k,root));
            }
        }
        return res;
    }

    template<class type>
    function<type> Arctgh(const function<type> & other, int k)
    {
        function<type> res(other);
        for(int i=0;i<other.getnx();i++)
        {
            for(int j=0;j<other.getny();j++)
            {
                res.setf(i,j,Arctgh(other.getf(i,j),k));
            }
        }
        return res;
    }

    template<class type>
    function<type> Arcctgh(const function<type> & other, int k)
    {
        function<type> res(other);
        for(int i=0;i<other.getnx();i++)
        {
            for(int j=0;j<other.getny();j++)
            {
                res.setf(i,j,Arcctgh(other.getf(i,j),k));
            }
        }
        return res;
    }
}


#endif
