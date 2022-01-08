#ifndef _FUNCTION_H_
#define _FUNCTION_H_

#define x_arg 1
#define y_arg 2

#include "complex.h"
#include "matrix.h"

enum domainState
{
    nonNegative = 0,
    positive = 1
};

namespace X
{
    template<typename type>
    struct point
    {
        point(type x, type f, int ix):x(x),f(f),ix(ix){}
        type x,f;
        int ix;
    };

    template <typename type>
    class function;

    template<typename T1,typename T2> function<T2> operator + (const T1 & C,const function<T2> & other);
    template<typename T1,typename T2> function<T2> operator - (const T1 & C,const function<T2> & other);
    template<typename T1,typename T2> function<T2> operator * (const T1 & C,const function<T2> & other);
    template<typename T1,typename T2> function<T2> operator / (const T1 & C,const function<T2> & other);

    template<typename T> std::ostream & operator<<(std::ostream & out, const function<T> & point);

    template<typename T> function<T> exp(const function<T> & other);
    template<typename T> function<T> ln(const function<T> & other);
    template<typename T> function<T> sqrt(const function<T> & other);
    template<typename T> function<T> pow(const function<T> & other, const T & deg);
    template<typename T> function<T> sin(const function<T> & other);
    template<typename T> function<T> cos(const function<T> & other);
    template<typename T> function<T> tg(const function<T> & other);
    template<typename T> function<T> ctg(const function<T> & other);
    template<typename T> function<T> sinh(const function<T> & other);
    template<typename T> function<T> cosh(const function<T> & other);
    template<typename T> function<T> tgh(const function<T> & other);
    template<typename T> function<T> ctgh(const function<T> & other);
    template<typename T> function<T> arcsin(const function<T> & other);
    template<typename T> function<T> arccos(const function<T> & other);
    template<typename T> function<T> arctg(const function<T> & other);
    template<typename T> function<T> arcctg(const function<T> & other);
    template<typename T> function<T> arcsinh(const function<T> & other);
    template<typename T> function<T> arccosh(const function<T> & other);
    template<typename T> function<T> arctgh(const function<T> & other);
    template<typename T> function<T> arcctgh(const function<T> & other);

    template<typename T> T integralS1(int x1, int x2, const function<T> & other);
    template<typename T> T integralS2(int x1, int x2, const function<T> & other);

    template<typename T> function<T> spline(function<T> f, int N);

    template <typename type>
    class function
    {
        private:
            type *x;
            type *f;
            int nx;
            int setprecision=1, setw=2;
            int SimpsonIntegralType = 2;

        public:
            function();
            function(type x0, type xn, int nx);
            function(type *x, int nx);
            function(type *x, type *f, int nx);
            function(int nx);
            function(const function<type> & other);
            ~function();

            type getf(int i) const;
            type getx(int i) const;
            int getn() const;
            void getx();
            void setf(int i, type f);
            void setx(int i, type x);
            void setout(int setprecision, int setw);

            int zeropoint() const;
            bool domainIsEqual(const function<type> & other) const;
            bool positive(domainState cond) const;

            bool operator == (const function<type> & other) const;
            bool operator != (const function<type> & other) const;
            function &operator = (const function<type> & other);
            function operator + (const function<type> & other) const;
            function operator - (const function<type> & other) const;
            function operator + (const type & other) const;
            function operator - (const type & other) const;
            function operator - () const;
            function operator * (const function<type> & other) const;
            function operator / (const function<type> & other) const;
            function operator * (const type & C) const;
            function operator / (const type & C) const;
            function operator ^ (const int & x) const;                                                                  //плохой оператор возведения в целую степень, поскольку имеет низкий приоритет
            template<typename T1,typename T2> friend function<T2> operator + (const T1 & C,const function<T2> & other);
            template<typename T1,typename T2> friend function<T2> operator - (const T1 & C,const function<T2> & other);
            template<typename T1,typename T2> friend function<T2> operator * (const T1 & C,const function<T2> & other);
            template<typename T1,typename T2> friend function<T2> operator / (const T1 & C,const function<T2> & other);
            type & operator () (const int i);

            template<typename T> friend std::ostream & operator<<(std::ostream & out, const function<T> & point);

            template<typename T> friend function<T> exp(const function<T> & other);
            template<typename T> friend function<T> ln(const function<T> & other);
            template<typename T> friend function<T> sqrt(const function<T> & other);
            template<typename T> friend function<T> pow(const function<T> & other, const T & deg);
            template<typename T> friend function<T> sin(const function<T> & other);
            template<typename T> friend function<T> cos(const function<T> & other);
            template<typename T> friend function<T> tg(const function<T> & other);
            template<typename T> friend function<T> ctg(const function<T> & other);
            template<typename T> friend function<T> sinh(const function<T> & other);
            template<typename T> friend function<T> cosh(const function<T> & other);
            template<typename T> friend function<T> tgh(const function<T> & other);
            template<typename T> friend function<T> ctgh(const function<T> & other);
            template<typename T> friend function<T> arcsin(const function<T> & other);
            template<typename T> friend function<T> arccos(const function<T> & other);
            template<typename T> friend function<T> arctg(const function<T> & other);
            template<typename T> friend function<T> arcctg(const function<T> & other);
            template<typename T> friend function<T> arcsinh(const function<T> & other);
            template<typename T> friend function<T> arccosh(const function<T> & other);
            template<typename T> friend function<T> arctgh(const function<T> & other);
            template<typename T> friend function<T> arcctgh(const function<T> & other);

            type integralS1(int x1, int x2) const;
            type integralS2(int x1, int x2) const;
            template<typename T> friend T integralS1(int x1, int x2, const function<T> & other);
            template<typename T> friend T integralS2(int x1, int x2, const function<T> & other);

            void normalize();
            type average() const;
            type dispersion() const;
            type moment() const;

            void normalizeS1();
            type averageS1() const;
            type dispersionS1() const;

            type max() const;
            type min() const;
            int maxIndex() const;
            int minIndex() const;
            point<type> maxPoint() const;
            point<type> minPoint() const;

            function & spline(int N);
            template<typename T> friend function<T> spline(function<T> f, int N);

            const char * savef(const char * name) const;
            const char * savex(const char * name) const;
    };
}

namespace XY
{
    template<typename type>
    struct point
    {
        point(type x, type y, type f, int ix, int iy):x(x),y(y),f(f),ix(ix),iy(iy){}
        type x,y,f;
        int ix,iy;
    };

    template <typename type>
    class function;

    template<typename T1,typename T2> function<T2> operator + (const T1 & C,const function<T2> & other);
    template<typename T1,typename T2> function<T2> operator - (const T1 & C,const function<T2> & other);
    template<typename T1,typename T2> function<T2> operator * (const T1 & C,const function<T2> & other);
    template<typename T1,typename T2> function<T2> operator / (const T1 & C,const function<T2> & other);

    template<typename T> std::ostream & operator<<(std::ostream & out, const function<T> & point);

    template<typename T> function<T> exp(const function<T> & other);
    template<typename T> function<T> ln(const function<T> & other);
    template<typename T> function<T> sqrt(const function<T> & other);
    template<typename T> function<T> pow(const function<T> & other, const T & deg);
    template<typename T> function<T> sin(const function<T> & other);
    template<typename T> function<T> cos(const function<T> & other);
    template<typename T> function<T> tg(const function<T> & other);
    template<typename T> function<T> ctg(const function<T> & other);
    template<typename T> function<T> sinh(const function<T> & other);
    template<typename T> function<T> cosh(const function<T> & other);
    template<typename T> function<T> tgh(const function<T> & other);
    template<typename T> function<T> ctgh(const function<T> & other);
    template<typename T> function<T> arcsin(const function<T> & other);
    template<typename T> function<T> arccos(const function<T> & other);
    template<typename T> function<T> arctg(const function<T> & other);
    template<typename T> function<T> arcctg(const function<T> & other);
    template<typename T> function<T> arcsinh(const function<T> & other);
    template<typename T> function<T> arccosh(const function<T> & other);
    template<typename T> function<T> arctgh(const function<T> & other);
    template<typename T> function<T> arcctgh(const function<T> & other);

    template<typename T> function<cn::complex<T>> exp(const function<cn::complex<T>> & other);

    template<typename T> T integralS11(int x1, int x2, int y1, int y2, const function<T> & other);
    template<typename T> T integralS22(int x1, int x2, int y1, int y2, const function<T> & other);
    template<typename T> T integralS21(int x1, int x2, int y1, int y2, const function<T> & other);
    template<typename T> T integralS12(int x1, int x2, int y1, int y2, const function<T> & other);

    template<typename T> function<T> spline(const function<T> & f, int Nx, int Ny);

    template <typename type>
    class function
    {
    private:
        type **x;
        type **y;
        type **f;
        int nx;
        int ny;
        int setprecision=1, setw=2;
    public:
        function();
        function(type x0, type xn, type y0, type yn, int nx, int ny, int isArg=0);
        function(type **x, type **y, int nx, int ny, int isArg=0);
        function(type **x, type **y, type **f, int nx, int ny);
        function(int nx, int ny);
        function(const function<type> & other);
        function(function<type> & other, int isArg=0);
        ~function();

        type getf(int i, int j) const;
        type getx(int i, int j) const;
        type gety(int i, int j) const;
        int getnx() const;
        int getny() const;
        function getx() const;
        function gety() const;
        void setf(int i, int j, type f);
        void setx(int i, int j, type x);
        void sety(int i, int j, type x);
        void setout(int setprecision, int setw);

        int zeropoint() const;
        bool domainIsEqual(const function<type> & other) const;
        bool positive(domainState cond) const;

        bool operator == (const function<type> & other) const;
        bool operator != (const function<type> & other) const;
        function &operator = (const function<type> & other);
        function operator + (const function<type> & other) const;
        function operator - (const function<type> & other) const;
        function operator + (const type & other) const;
        function operator - (const type & other) const;
        function operator - () const;
        function operator * (const function<type> & other) const;
        function operator / (const function<type> & other) const;
        function operator * (const type & C) const;
        function operator / (const type & C) const;
        template<typename T1,typename T2> friend function<T2> operator + (const T1 & C,const function<T2> & other);
        template<typename T1,typename T2> friend function<T2> operator - (const T1 & C,const function<T2> & other);
        template<typename T1,typename T2> friend function<T2> operator * (const T1 & C,const function<T2> & other);
        template<typename T1,typename T2> friend function<T2> operator / (const T1 & C,const function<T2> & other);
        type & operator () (const int i, const int j);

        template<typename T> friend std::ostream & operator<<(std::ostream & out, const function<T> & point);

        template<typename T> friend function<T> exp(const function<T> & other);
        template<typename T> friend function<T> ln(const function<T> & other);
        template<typename T> friend function<T> sqrt(const function<T> & other);
        template<typename T> friend function<T> pow(const function<T> & other, const T & deg);
        template<typename T> friend function<T> sin(const function<T> & other);
        template<typename T> friend function<T> cos(const function<T> & other);
        template<typename T> friend function<T> tg(const function<T> & other);
        template<typename T> friend function<T> ctg(const function<T> & other);
        template<typename T> friend function<T> sinh(const function<T> & other);
        template<typename T> friend function<T> cosh(const function<T> & other);
        template<typename T> friend function<T> tgh(const function<T> & other);
        template<typename T> friend function<T> ctgh(const function<T> & other);
        template<typename T> friend function<T> arcsin(const function<T> & other);
        template<typename T> friend function<T> arccos(const function<T> & other);
        template<typename T> friend function<T> arctg(const function<T> & other);
        template<typename T> friend function<T> arcctg(const function<T> & other);
        template<typename T> friend function<T> arcsinh(const function<T> & other);
        template<typename T> friend function<T> arccosh(const function<T> & other);
        template<typename T> friend function<T> arctgh(const function<T> & other);
        template<typename T> friend function<T> arcctgh(const function<T> & other);

        template<typename T> friend function<cn::complex<T>> exp(const function<cn::complex<T>> & other);

        type integralS11(int x1, int x2, int y1, int y2) const;
        type integralS22(int x1, int x2, int y1, int y2) const;
        type integralS21(int x1, int x2, int y1, int y2) const;
        type integralS12(int x1, int x2, int y1, int y2) const;
        template<typename T> friend T integralS11(int x1, int x2, int y1, int y2, const function<T> & other);
        template<typename T> friend T integralS22(int x1, int x2, int y1, int y2, const function<T> & other);
        template<typename T> friend T integralS21(int x1, int x2, int y1, int y2, const function<T> & other);
        template<typename T> friend T integralS12(int x1, int x2, int y1, int y2, const function<T> & other);

        function & spline(int Nx, int Ny);
        template<typename T> friend function<T> spline(const function<T> & f, int Nx, int Ny);

        point<type> maxPoint() const;
        point<type> minPoint() const;

        const char * savef(const char * name) const;
        const char * savex(const char * name) const;
        const char * savey(const char * name) const;
    };
}

namespace XYZ
{
    template <typename type, int nx, int ny>
    class function
    {
    public:
        function(type ***x, type ***y, type ***z);
        function(type ***x, type ***y, type ***z, type ***f);
        ~function();
    };
}

namespace X {
    template <>
    class function<char>{};
}

namespace XY {
    template <>
    class function<char>{};
}

//============================================================================================================================================//
//============================================================================================================================================//
//============================================================================================================================================//

namespace X {

    template<typename type>
    X::function<type>::function()
    {
        x=nullptr;
        f=nullptr;
    }

    template<typename type>
    X::function<type>::function(type x0, type xn, int nx)
    {
        this->x=new type[nx];
        this->f=new type[nx];
        this->nx=nx;
        type h = (xn-x0)/(nx-1);
        for(int i=0;i<nx;i++)
        {
            this->x[i]=x0+h*i;
            this->f[i]=x[i];
        }
    }

    template<typename type>
    X::function<type>::function(type *x, int nx)
    {
        this->x=new type[nx];
        this->f=new type[nx];
        this->nx=nx;
        for(int i=0;i<nx;i++)
        {
            this->x[i]=x[i];
            this->f[i]=x[i];
        }
    }

    template<typename type>
    X::function<type>::function(type *x, type *f, int nx)
    {
        this->x=new type[nx];
        this->f=new type[nx];
        this->nx=nx;
        for(int i=0;i<nx;i++)
        {
            this->x[i]=x[i];
            this->f[i]=f[i];
        }
    }

    template<typename type>
    X::function<type>::function(int nx)
    {
        this->x=new type[nx];
        this->f=new type[nx];
        this->nx=nx;
    }

    template<typename type>
    X::function<type>::function(const X::function<type> &other)
    {
        this->setprecision=other.setprecision;
        this->setw=other.setw;
        this->nx=other.nx;
        this->x = new type[nx];
        this->f = new type[nx];
        for(int i=0;i<nx;i++)
        {
            this->x[i]=other.x[i];
            this->f[i]=other.f[i];
        }
    }

    template<typename type>
    X::function<type>::~function()
    {
        if(this->x!=nullptr)
        {
            delete[] this->x;
        }
        if(this->f!=nullptr)
        {
            delete[] this->f;
        }
    }

    template<typename type>
    type function<type>::getx(int i) const
    {
        return x[i];
    }

    template<typename type>
    int function<type>::getn() const
    {
        return nx;
    }

    template<typename type>
    void function<type>::getx()
    {
        for(int i=0;i<nx;i++)
        {
            std::cout <<std::setw(setw)<<std::fixed<<std::setprecision(setprecision)<<x[i]<<"  ";
        }
    }

    template<typename type>
    type function<type>::getf(int i) const
    {
        return f[i];
    }

    template<typename type>
    void function<type>::setx(int i, type x)
    {
        this->x[i]=x;
    }

    template<typename type>
    void function<type>::setf(int i, type f)
    {
        this->f[i]=f;
    }

    template<typename type>
    void function<type>::setout(int setprecision, int setw)
    {
        this->setprecision=setprecision;
        this->setw=setw;
    }

    template<typename type>
    function<type> &function<type>::operator =(const function<type> &other)
    {
        if(this->x!=nullptr)
        {
            delete[] this->x;
        }
        if(this->f!=nullptr)
        {
            delete[] this->f;
        }
        this->setprecision=other.setprecision;
        this->setw=other.setw;
        this->nx=other.nx;
        this->x=new type[nx];
        this->f=new type[nx];
        for(int i=0;i<nx;i++)
        {
            this->x[i]=other.x[i];
            this->f[i]=other.f[i];
        }
        return *this;
    }

    template<typename type>
    bool X::function<type>::operator == (const function<type> & other) const
    {
        for(int i=0;i<nx;i++)
        {
            if((this->f[i]==other.f[i])&&(this->x[i]==other.x[i])){}
            else{return false;}
        }
        return true;
    }

    template<typename type>
    bool X::function<type>::operator != (const function<type> & other) const
    {
        for(int i=0;i<nx;i++)
        {
            if((this->f[i]==other.f[i])&&(this->x[i]==other.x[i])){}
            else{return true;}
        }
        return false;
    }

    template<typename type>
    bool X::function<type>::domainIsEqual(const function &other) const
    {
        for(int i=0;i<nx;i++)
        {
            if(this->x[i]==other.x[i]){}
            else{return false;}
        }
        return true;
    }

    template<typename type>
    int function<type>::zeropoint() const
    {
        for(int i=0;i<nx;i++)
        {
            if(this->f[i]==0){return i;}
            else{}
        }
        return -1;
    }

    template<typename type>
    bool X::function<type>::positive(domainState cond) const
    {
        switch (cond) {
            case domainState::nonNegative:
            {
                for(int i=0;i<nx;i++)
                {
                    if(this->f[i]>=0){}
                    else{return false;}
                }
                break;
            }
            case domainState::positive:
            {
                for(int i=0;i<nx;i++)
                {
                    if(this->f[i]>0){}
                    else{return false;}
                }
                break;
            }
        }
        return true;
    }

    template<typename type>
    X::function<type> X::function<type>::operator +(const function &other) const
    {
        function<type> res(*this);
        #ifdef domaincheck
            if(domainIsEqual(other)){
                for(int i=0;i<nx;i++)
                {
                    res.f[i]+=other.f[i];
                }
            }
        #else
            for(int i=0;i<nx;i++)
            {
                res.f[i]+=other.f[i];
            }
        #endif
        return res;
    }

    template<typename type>
    X::function<type> X::function<type>::operator -(const function &other) const
    {
        function<type> res(*this);
        #ifdef domaincheck
            if(domainIsEqual(other)){
                for(int i=0;i<nx;i++)
                {
                    res.f[i]-=other.f[i];
                }
            }
        #else
            for(int i=0;i<nx;i++)
            {
                res.f[i]-=other.f[i];
            }
        #endif
        return res;
    }

    template<typename type>
    X::function<type> X::function<type>::operator +(const type &other) const
    {
        function<type> res(*this);
        #ifdef domaincheck
            if(domainIsEqual(other)){
                for(int i=0;i<nx;i++)
                {
                    res.f[i]+=other;
                }
            }
        #else
            for(int i=0;i<nx;i++)
            {
                res.f[i]+=other;
            }
        #endif
        return res;
    }

    template<typename type>
    X::function<type> X::function<type>::operator -(const type &other) const
    {
        function<type> res(*this);
        #ifdef domaincheck
            if(domainIsEqual(other)){
                for(int i=0;i<nx;i++)
                {
                    res.f[i]-=other;
                }
            }
        #else
            for(int i=0;i<nx;i++)
            {
                res.f[i]-=other;
            }
        #endif
        return res;
    }

    template<typename type>
    X::function<type> X::function<type>::operator -() const
    {
        function<type> res(*this);
        for(int i=0;i<nx;i++)
        {
            res.f[i]=res.f[i]*(-1);
        }
        return res;
    }

    template<typename type>
    X::function<type> X::function<type>::operator *(const function &other) const
    {
        function<type> res(*this);
        #ifdef domaincheck
            if(domainIsEqual(other)){
                for(int i=0;i<nx;i++)
                {
                    res.f[i]*=other.f[i];
                }
            }
        #else
            for(int i=0;i<nx;i++)
            {
                res.f[i]*=other.f[i];
            }
        #endif
        return res;
    }

    template<typename type>
    X::function<type> X::function<type>::operator /(const function &other) const
    {
        function<type> res(*this);
        #ifdef domaincheck
            if(domainIsEqual(other)&&(zeropoint()==-1)){
                for(int i=0;i<nx;i++)
                {
                    res.f[i]/=other.f[i];
                }
            }
        #else
            for(int i=0;i<nx;i++)
            {
                res.f[i]/=other.f[i];
            }
        #endif
        return res;
    }

    template<typename type>
    X::function<type> X::function<type>::operator *(const type &C) const
    {
        function<type> res(*this);
        for(int i=0;i<nx;i++)
        {
            res.f[i]*=C;
        }
        return res;
    }

    template<typename type>
    X::function<type> X::function<type>::operator /(const type & C) const
    {
        function<type> res(*this);
        for(int i=0;i<nx;i++)
        {
            res.f[i]=res.f[i]/C;
        }
        return res;
    }

    template<typename type>
    function<type> function<type>::operator ^(const int & x) const
    {
        function<type> res(*this);
        for(int i=0;i<nx;i++)
        {
            for(int d=1; d<x;d++)
            {
                res.f[i] *= f[i];
            }
        }
        return res;
    }

    template<typename T1,typename T2>
    function<T2> operator + (const T1 & C,const function<T2> & other)
    {
        function<T2> res(other);
        for(int i=0;i<other.nx;i++)
        {
            res.f[i]=C;
            res.f[i]+=other.f[i];
        }
        return res;
    }

    template<typename T1,typename T2>
    function<T2> operator - (const T1 & C,const function<T2> & other)
    {
        function<T2> res(other);
        for(int i=0;i<other.nx;i++)
        {
            res.f[i]=C;
            res.f[i]-=other.f[i];
        }
        return res;
    }

    template<typename T1,typename T2>
    function<T2> operator * (const T1 & C,const function<T2> & other)
    {
        function<T2> res(other);
        for(int i=0;i<other.nx;i++)
        {
            res.f[i]=C;
            res.f[i]*=other.f[i];
        }
        return res;
    }

    template<typename T1,typename T2>
    function<T2> operator / (const T1 & C,const function<T2> & other)
    {
        function<T2> res(other);
        for(int i=0;i<other.nx;i++)
        {
            res.f[i]=C;
            res.f[i]/=other.f[i];
        }
        return res;
    }

    template<typename type>
    type & function<type>::operator ()(const int i)
    {
        return f[i];
    }

    template<typename T>
    std::ostream & operator<<(std::ostream &out, const X::function<T> &point)
    {
        for(int i=0;i<point.nx;i++)
        {
            out <<std::setw(point.setw)<<std::fixed<<std::setprecision(point.setprecision)<<point.f[i]<<"  ";
        }
        return out;
    }

    template<typename T>
    function<T> exp(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            res.setf(i,
                     std::exp(other.getf(i)));
        }
        return res;
    }

    template<typename T>
    function<T> ln(const function<T> & other)
    {
        function<T> res(other);
        #ifdef domaincheck
            if(res.positive(domainState::positive))
            {
                for(int i=0;i<other.nx;i++)
                {
                    res.setf(i,
                             std::log(other.getf(i)));
                }
            }else
            {
                std::cout<<"negative or zero argument, was returned argument"<<std::endl;
            }
        #else
            for(int i=0;i<other.nx;i++)
            {
                res.setf(i,
                         std::log(other.getf(i)));
            }
        #endif
        return res;
    }

    template<typename T>
    function<T> sqrt(const function<T> & other)
    {
        function<T> res(other);
        #ifdef domaincheck
            if(res.positive(domainState::nonNegative))
            {
                for(int i=0;i<other.nx;i++)
                {
                    res.setf(i,
                             std::sqrt(other.getf(i)));
                }
            }else
            {
                std::cout<<"negative argument, was returned argument"<<std::endl;
            }
        #else
            for(int i=0;i<other.nx;i++)
            {
                res.setf(i,
                         std::sqrt(other.getf(i)));
            }
        #endif
        return res;
    }

    template<typename T>
    function<T> pow(const function<T> & other, const T & deg)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            res.setf(i,
                     std::pow(other.getf(i),deg));
        }
        return res;
    }

    template<typename T>
    function<T> sin(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            res.setf(i,
                     std::sin(other.getf(i)));
        }
        return res;
    }

    template<typename T>
    function<T> cos(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            res.setf(i,
                     std::cos(other.getf(i)));
        }
        return res;
    }

    template<typename T>
    function<T> tg(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            res.setf(i,
                     std::tan(other.getf(i)));
        }
        return res;
    }

    template<typename T>
    function<T> ctg(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            res.setf(i,
                     1/std::tan(other.getf(i)));
        }
        return res;
    }

    template<typename T>
    function<T> sinh(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            res.setf(i,
                     std::sinh(other.getf(i)));
        }
        return res;
    }

    template<typename T>
    function<T> cosh(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            res.setf(i,
                     std::cosh(other.getf(i)));
        }
        return res;
    }

    template<typename T>
    function<T> tgh(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            res.setf(i,
                     std::tanh(other.getf(i)));
        }
        return res;
    }

    template<typename T>
    function<T> ctgh(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            res.setf(i,
                     1/std::tanh(other.getf(i)));
        }
        return res;
    }

    template<typename T>
    function<T> arcsin(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            res.setf(i,
                     std::asin(other.getf(i)));
        }
        return res;
    }

    template<typename T>
    function<T> arccos(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            res.setf(i,
                     std::acos(other.getf(i)));
        }
        return res;
    }

    template<typename T>
    function<T> arctg(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            res.setf(i,
                     std::atan(other.getf(i)));
        }
        return res;
    }

    template<typename T>
    function<T> arcctg(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            res.setf(i,
                     1/std::atan(other.getf(i)));
        }
        return res;
    }

    template<typename T>
    function<T> arcsinh(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            res.setf(i,
                     std::asinh(other.getf(i)));
        }
        return res;
    }

    template<typename T>
    function<T> arccosh(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            res.setf(i,
                     std::acosh(other.getf(i)));
        }
        return res;
    }

    template<typename T>
    function<T> arctgh(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            res.setf(i,
                     std::atanh(other.getf(i)));
        }
        return res;
    }

    template<typename T>
    function<T> arcctgh(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            res.setf(i,
                     1/std::atanh(other.getf(i)));
        }
        return res;
    }

    template<typename type>
    type function<type>::integralS1(int x1, int x2) const
    {
        type I=0;
        int N=x2-x1;
        type *h=new type[N];
        for(int i=0;i<N;i++)
        {
            h[i]=x[x1+i+1]-x[x1+i];
        }

        for(int i=0;i<(N-1)/2;i++)
        {
            type S=(-h[2*i+1]/h[2*i]+2)*f[x1+2*i]+
                    (h[2*i]+h[2*i+1])*(h[2*i]+h[2*i+1])/h[2*i]/h[2*i+1]*f[x1+2*i+1]+
                    (-h[2*i]/h[2*i+1]+2)*f[x1+2*i+2];
            S=S*(h[2*i]+h[2*i+1])/6;
            I=I+S;
        }
        type a = (h[N-1]*h[N-1]*2+h[N-1]*h[N-2]*3)/6/(h[N-2]+h[N-1]);
        type b = (h[N-1]*h[N-1]+h[N-1]*h[N-2]*3)/6/(h[N-2]);
        type c = h[N-1]*h[N-1]*h[N-1]/6/h[N-2]/(h[N-2]+h[N-1]);
        I+=(a*f[x1+N]+b*f[x1+N-1]-c*f[x1+N-2]);

        delete [] h;

        return I;
    }

    template<typename type>
    type function<type>::integralS2(int x1, int x2) const
    {
        type I=0;
        int N=x2-x1;
        type *h=new type[N];
        for(int i=0;i<N;i++)
        {
            h[i]=x[x1+i+1];
            h[i]-=x[x1+i];
        }

        for(int i=0;i<N/2;i++)
        {
            type S=(-h[2*i+1]/h[2*i]+2); S*=f[x1+2*i];
            type arg = (h[2*i]+h[2*i+1]); arg *= h[2*i]+h[2*i+1]; arg/=h[2*i]; arg/=h[2*i+1]; arg*=f[x1+2*i+1];
            S+=arg;
            arg = (-h[2*i]/h[2*i+1]+2); arg *=f[x1+2*i+2];
            S+=arg;
            S*=(h[2*i]+h[2*i+1]);S/=6;
            I+=S;
        }

        delete [] h;

        return I;
    }

    template<typename T>
    T integralS1(int x1, int x2, const function<T> & other)
    {
        return other.integralS1(x1,x2);
    }

    template<typename T>
    T integralS2(int x1, int x2, const function<T> & other)
    {
        return other.integralS2(x1,x2);
    }

    template<typename type>
    void function<type>::normalize()
    {
        type N = this->integralS2(0,nx-1);
        for(int i=0;i<nx;i++)
        {
            f[i]=f[i]/N;
        }
    }

    template<typename type>
    type function<type>::average() const
    {
        function<type> buf(*this);
        for(int i=0;i<nx;i++)
        {
            buf(i)=x[i]*f[i];
        }
        return buf.integralS2(0,nx-1);
    }

    template<typename type>
    type function<type>::dispersion() const
    {
        function<type> buf(*this);
        type E = buf.average();
        for(int i=0;i<nx;i++)
        {
            buf(i)=(x[i]-E)*(x[i]-E)*f[i];
        }
        return buf.integralS2(0,nx-1);
    }

    template<typename type>
    type function<type>::moment() const
    {
        function<type> buf(*this);
        type E = buf.average();
        for(int i=0;i<nx;i++)
        {
            buf(i)=(x[i]-E)*(x[i]-E)*(x[i]-E)*(x[i]-E)*f[i];
        }
        return buf.integralS2(0,nx-1);
    }

    template<typename type>
    void function<type>::normalizeS1()
    {
        type N = this->integralS1(0,nx-1);
        for(int i=0;i<nx;i++)
        {
            f[i]=f[i]/N;
        }
    }

    template<typename type>
    type function<type>::averageS1() const
    {
        function<type> buf(*this);
        for(int i=0;i<nx;i++)
        {
            buf(i)=x[i]*f[i];
        }
        return buf.integralS1(0,nx-1);
    }

    template<typename type>
    type function<type>::dispersionS1() const
    {
        function<type> buf(*this);
        type E = buf.average();
        for(int i=0;i<nx;i++)
        {
            buf(i)=(x[i]-E)*(x[i]-E)*f[i];
        }
        return buf.integralS1(0,nx-1);
    }

    template<typename type>
    type function<type>::max() const
    {
        type max = f[0];

        for(int i=0;i<nx;i++)
        {
            if(f[i]>max)
            {
                max=f[i];
            }
        }
        return max;
    }

    template<typename type>
    type function<type>::min() const
    {
        type min = f[0];

        for(int i=0;i<nx;i++)
        {
            if(f[i]<min)
            {
                min=f[i];
            }
        }
        return min;
    }

    template<typename type>
    int function<type>::maxIndex() const
    {
        type max = f[0];
        int point;

        for(int i=0;i<nx;i++)
        {
            if(f[i]>max)
            {
                max=f[i];
                point=i;
            }
        }
        return point;
    }

    template<typename type>
    int function<type>::minIndex() const
    {
        type min = f[0];
        int point;

        for(int i=0;i<nx;i++)
        {
            if(f[i]<min)
            {
                min=f[i];
                point=i;
            }
        }
        return point;
    }


    template<typename type>
    point<type> function<type>::maxPoint() const
    {
        point<type> max(x[0],f[0],0);

        for(int i=0;i<nx;i++)
        {
            if(f[i]>max.f)
            {
                max.f=f[i];
                max.x=x[i];
                max.ix=i;
            }
        }
        return max;
    }

    template<typename type>
    point<type> function<type>::minPoint() const
    {
        point<type> min(x[0],f[0],0);

        for(int i=0;i<nx;i++)
        {
            if(f[i]<min.f)
            {
                min.f=f[i];
                min.x=x[i];
                min.ix=i;
            }
        }
        return min;
    }

    template<typename type>
    const char *function<type>::savef(const char *name) const
    {
        FILE *fp;
        fp=freopen(name,"w",stdout);
        for(int i=0;i<nx;i++)
        {
            std::cout <<std::setw(setw)<<std::fixed<<std::setprecision(setprecision)<<f[i]<<"  ";
        }
        fclose(fp);
        freopen("CON","w",stdout);
        return name;
    }

    template<typename type>
    const char *function<type>::savex(const char *name) const
    {
        FILE *fp;
        fp=freopen(name,"w",stdout);
        for(int i=0;i<nx;i++)
        {
            std::cout <<std::setw(setw)<<std::fixed<<std::setprecision(setprecision)<<x[i]<<"  ";
        }
        fclose(fp);
        freopen("CON","w",stdout);
        return name;
    }

    template<typename type>
    function<type> & function<type>::spline(int N)
    {
        type *a = new type[nx];
        type *b = new type[nx];
        type *c = new type[nx];
        type *d = new type[nx];
        type *h = new type[nx];
        for(int i=1;i<nx;i++)
        {
            h[i]=x[i]-x[i-1];
            a[i]=f[i];
        }
        a[0]=f[0];
        c[0]=0;
        c[nx-1]=0;
        matrix<type> C(f-2,f-1);

        for(int i=2;i<f-2;i++)
        {
            C(i-1,i-2)=h[i];
            C(i-1,i-1)=(h[i]+h[i+1])*2.0;
            C(i-1,i)=h[i+1];
            C(i-1,f-2)=((a[i+1]-a[i])/h[i+1]-(a[i]-a[i-1])/h[i])*3.0;
        }
        C(0,0)=(h[1]+h[2])*2.0;
        C(0,1)=h[2];
        C(0,f-2)=((a[2]-a[1])/h[2]-(a[2]-a[1])/h[1])*3.0;
        C(f-3,f-4)=h[f-2];
        C(f-3,f-3)=(h[f-2]+h[f-1])*2.0;
        C(f-3,f-2)=((a[f-1]-a[f-2])/h[f-1]-(a[f-2]-a[f-3])/h[f-2])*3.0;
        C.tridiagonal();
        for(int i=1;i<f-1;i++)
        {
            c[i]=C(i-1,f-2);
        }
        for(int i=1;i<f;i++)
        {
            d[i]=(c[i]-c[i-1])/(h[i]*3);
            b[i]=(a[i]-a[i-1])/h[i]+(c[i]*2+c[i-1])*h[i]/3;
        }
        delete [] h;

        type *xnew = new type[N*(f-1)+1];
        xnew[0]=x[0];
        int k=1;
        for(int i=1;i<f;i++)
        {
            for(int j=0;j<N;j++)
            {
                xnew[k]=x[i-1]+(x[i]-x[i-1])*(j+1)/N;
                k++;
            }
        }
        function<type> res(xnew,N*(nx-1)+1);
        delete [] xnew;

        res(0)=a[0];
        k=1;
        for(int i=1;i<nx;i++)
        {
            for(int j=0;j<N;j++)
            {
                res(k)=a[i]+b[i]*(res.getx(k)-x[i])+c[i]*(res.getx(k)-x[i])*(res.getx(k)-x[i])+d[i]*(res.getx(k)-x[i])*(res.getx(k)-x[i])*(res.getx(k)-x[i]);
                k++;
            }
        }
        delete [] a;
        delete [] b;
        delete [] c;
        delete [] d;
        *this = res;
        return *this;
    }

    template<typename T>
    function<T> spline(function<T> f, int N)
    {
        T *a = new T[f.nx];
        T *b = new T[f.nx];
        T *c = new T[f.nx];
        T *d = new T[f.nx];
        T *h = new T[f.nx];
        for(int i=1;i<f.nx;i++)
        {
            h[i]=f.x[i]-f.x[i-1];
            a[i]=f.f[i];
        }
        a[0]=f.f[0];
        c[0]=0;
        c[f.nx-1]=0;
        matrix<T> C(f.nx-2,f.nx-1);

        for(int i=2;i<f.nx-2;i++)
        {
            C(i-1,i-2)=h[i];
            C(i-1,i-1)=(h[i]+h[i+1])*2.0;
            C(i-1,i)=h[i+1];
            C(i-1,f.nx-2)=((a[i+1]-a[i])/h[i+1]-(a[i]-a[i-1])/h[i])*3.0;
        }
        C(0,0)=(h[1]+h[2])*2.0;
        C(0,1)=h[2];
        C(0,f.nx-2)=((a[2]-a[1])/h[2]-(a[2]-a[1])/h[1])*3.0;
        C(f.nx-3,f.nx-4)=h[f.nx-2];
        C(f.nx-3,f.nx-3)=(h[f.nx-2]+h[f.nx-1])*2.0;
        C(f.nx-3,f.nx-2)=((a[f.nx-1]-a[f.nx-2])/h[f.nx-1]-(a[f.nx-2]-a[f.nx-3])/h[f.nx-2])*3.0;
        C.tridiagonal();
        for(int i=1;i<f.nx-1;i++)
        {
            c[i]=C(i-1,f.nx-2);
        }
        for(int i=1;i<f.nx;i++)
        {
            d[i]=(c[i]-c[i-1])/(h[i]*3);
            b[i]=(a[i]-a[i-1])/h[i]+(c[i]*2+c[i-1])*h[i]/3;
        }
        delete [] h;

        T *xnew = new T[N*(f.nx-1)+1];
        xnew[0]=f.x[0];
        int k=1;
        for(int i=1;i<f.nx;i++)
        {
            for(int j=0;j<N;j++)
            {
                xnew[k]=f.x[i-1]+(f.x[i]-f.x[i-1])*(j+1)/N;
                k++;
            }
        }
        function<T> res(xnew,N*(f.nx-1)+1);
        delete [] xnew;

        res(0)=a[0];
        k=1;
        for(int i=1;i<f.nx;i++)
        {
            for(int j=0;j<N;j++)
            {
                res(k)=a[i]+b[i]*(res.getx(k)-f.getx(i))+c[i]*(res.getx(k)-f.getx(i))*(res.getx(k)-f.getx(i))+d[i]*(res.getx(k)-f.getx(i))*(res.getx(k)-f.getx(i))*(res.getx(k)-f.getx(i));
                k++;
            }
        }
        delete [] a;
        delete [] b;
        delete [] c;
        delete [] d;
        return res;
    }
}

namespace XY {

    template<typename type>
    function<type>::function()
    {
        x=nullptr;
        y=nullptr;
        f=nullptr;
    }

    template<typename type>
    function<type>::function(type x0, type xn, type y0, type yn, int nx, int ny, int isArg)
    {
        type hx = (xn-x0)/(nx-1);
        type hy = (yn-y0)/(ny-1);
        this->nx=nx;
        this->ny=ny;
        this->x=new type*[nx];
        this->y=new type*[nx];
        this->f=new type*[nx];
        for(int i=0;i<nx;i++)
        {
            this->x[i]=new type[ny];
            this->y[i]=new type[ny];
            this->f[i]=new type[ny];
            for(int j=0;j<ny;j++)
            {
                this->x[i][j]=x0+hx*i;
                this->y[i][j]=y0+hy*j;
            }
        }
        switch (isArg) {
            case 0:
            {
                for(int i=0;i<nx;i++)
                {
                    for(int j=0;j<ny;j++)
                    {
                        this->f[i][j]=0;
                    }
                }
                break;
            }
            case 1:
            {
                for(int i=0;i<nx;i++)
                {
                    for(int j=0;j<ny;j++)
                    {
                        this->f[i][j]=x[i][j];
                    }
                }
                break;
            }
            case 2:
            {
                for(int i=0;i<nx;i++)
                {
                    for(int j=0;j<ny;j++)
                    {
                        this->f[i][j]=y[i][j];
                    }
                }
                break;
            }
        }
    }

    template<typename type>
    function<type>::function(type **x, type **y, int nx, int ny, int isArg)
    {
        this->nx=nx;
        this->ny=ny;
        this->x=new type*[nx];
        this->y=new type*[nx];
        this->f=new type*[nx];
        for(int i=0;i<nx;i++)
        {
            this->x[i]=new type[ny];
            this->y[i]=new type[ny];
            this->f[i]=new type[ny];
            for(int j=0;j<ny;j++)
            {
                this->x[i][j]=x[i][j];
                this->y[i][j]=y[i][j];
            }
        }
        switch (isArg) {
            case 0:
            {
                for(int i=0;i<nx;i++)
                {
                    for(int j=0;j<ny;j++)
                    {
                        this->f[i][j]=0;
                    }
                }
                break;
            }
            case 1:
            {
                for(int i=0;i<nx;i++)
                {
                    for(int j=0;j<ny;j++)
                    {
                        this->f[i][j]=x[i][j];
                    }
                }
                break;
            }
            case 2:
            {
                for(int i=0;i<nx;i++)
                {
                    for(int j=0;j<ny;j++)
                    {
                        this->f[i][j]=y[i][j];
                    }
                }
            }
            break;
        }
    }

    template<typename type>
    function<type>::function(type **x, type **y, type **f, int nx, int ny)
    {
        this->nx=nx;
        this->ny=ny;
        this->x=new type*[nx];
        this->y=new type*[nx];
        this->f=new type*[nx];
        for(int i=0;i<nx;i++)
        {
            this->x[i]=new type[ny];
            this->y[i]=new type[ny];
            this->f[i]=new type[ny];
            for(int j=0;j<ny;j++)
            {
                this->x[i][j]=x[i][j];
                this->y[i][j]=y[i][j];
                this->f[i][j]=f[i][j];
            }
        }
    }

    template<typename type>
    function<type>::function(int nx, int ny)
    {
        this->nx=nx;
        this->ny=ny;
        this->x=new type*[nx];
        this->y=new type*[nx];
        this->f=new type*[nx];
        for(int i=0;i<nx;i++)
        {
            this->x[i]=new type[ny];
            this->y[i]=new type[ny];
            this->f[i]=new type[ny];
        }
    }

    template<typename type>
    function<type>::function(const function<type> & other)
    {
        this->setprecision=other.setprecision;
        this->setw=other.setw;
        this->nx=other.nx;
        this->ny=other.ny;
        this->x=new type*[nx];
        this->y=new type*[nx];
        this->f=new type*[nx];
        for(int i=0;i<nx;i++)
        {
            this->x[i]=new type[ny];
            this->y[i]=new type[ny];
            this->f[i]=new type[ny];
            for(int j=0;j<ny;j++)
            {
                this->x[i][j]=other.x[i][j];
                this->y[i][j]=other.y[i][j];
                this->f[i][j]=other.f[i][j];
            }
        }
    }

    template<typename type>
    function<type>::function(function<type> & other, int isArg)
    {
        this->setprecision=other.setprecision;
        this->setw=other.setw;
        this->nx=other.nx;
        this->ny=other.ny;
        this->x=new type*[nx];
        this->y=new type*[nx];
        this->f=new type*[nx];
        for(int i=0;i<nx;i++)
        {
            this->x[i]=new type[ny];
            this->y[i]=new type[ny];
            this->f[i]=new type[ny];
            for(int j=0;j<ny;j++)
            {
                this->x[i][j]=other.x[i][j];
                this->y[i][j]=other.y[i][j];
            }
        }
        switch (isArg) {
            case 0:
            {
                for(int i=0;i<nx;i++)
                {
                    for(int j=0;j<ny;j++)
                    {
                        this->f[i][j]=other.f[i][j];
                    }
                }
                break;
            }
            case 1:
            {
                for(int i=0;i<nx;i++)
                {
                    for(int j=0;j<ny;j++)
                    {
                        this->f[i][j]=x[i][j];
                    }
                }
                break;
            }
            case 2:
            {
                for(int i=0;i<nx;i++)
                {
                    for(int j=0;j<ny;j++)
                    {
                        this->f[i][j]=y[i][j];
                    }
                }
            }
            break;
        }
    }

    template<typename type>
    function<type>::~function()
    {
        if(x!=nullptr)
        {
            for(int i=0;i<nx;i++)
            {
                delete [] x[i];
            }
            delete[] x;
        }
        if(y!=nullptr)
        {
            for(int i=0;i<nx;i++)
            {
                delete [] y[i];
            }
            delete[] y;
        }
        if(f!=nullptr)
        {
            for(int i=0;i<nx;i++)
            {
                delete [] f[i];
            }
            delete[] f;
        }
    }

    template<typename type>
    type function<type>::getf(int i, int j) const
    {
        return f[i][j];
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
    function<type> function<type>::getx() const
    {
        function<type> res(x,y,nx,ny,x_arg);
        return res;
    }

    template<typename type>
    function<type> function<type>::gety() const
    {
        function<type> res(x,y,nx,ny,y_arg);
        return res;
    }

    template<typename type>
    void function<type>::setf(int i, int j, type f)
    {
        this->f[i][j]=f;
    }

    template<typename type>
    void function<type>::setx(int i, int j, type x)
    {
        this->x[i][j]=x;
    }

    template<typename type>
    void function<type>::sety(int i, int j, type x)
    {
        this->y[i][j]=x;
    }

    template<typename type>
    void function<type>::setout(int setprecision, int setw)
    {
        this->setprecision=setprecision;
        this->setw=setw;
    }

    template<typename type>
    bool function<type>::domainIsEqual(const function &other) const
    {
        for(int i=0;i<nx;i++)
        {
            for(int j=0;j<ny;j++)
            {
                if(!(this->x[i][j]==other.x[i][j]&&this->y[i][j]==other.y[i][j]))
                {
                    return false;
                }
            }
        }
        return true;
    }

    template<typename type>
    int function<type>::zeropoint() const
    {
        for(int i=0;i<nx;i++)
        {
            for(int j=0;j<ny;j++)
            {
                if(this->f[i][j]==0){return 0;}
            }
        }
        return -1;
    }

    template<typename type>
    bool function<type>::positive(domainState cond) const
    {
        switch (cond) {
            case domainState::nonNegative:
            {
                for(int i=0;i<nx;i++)
                {
                    for(int j=0;j<ny;j++)
                    {
                        if(!(this->f[i][j]>=0)){return false;}
                    }
                }
                break;
            }
            case domainState::positive:
            {
                for(int i=0;i<nx;i++)
                {
                    for(int j=0;j<ny;j++)
                    {
                        if(!(this->f[i][j]>0)){return false;}
                    }
                }
                break;
            }
        }
        return true;
    }

    template<typename type>
    function<type> &function<type>::operator =(const function<type> &other)
    {
        if(x!=nullptr)
        {
            for(int i=0;i<nx;i++)
            {
                delete [] x[i];
            }
            delete[] x;
        }
        if(y!=nullptr)
        {
            for(int i=0;i<nx;i++)
            {
                delete [] y[i];
            }
            delete[] y;
        }
        if(f!=nullptr)
        {
            for(int i=0;i<nx;i++)
            {
                delete [] f[i];
            }
            delete[] f;
        }

        this->setprecision=other.setprecision;
        this->setw=other.setw;

        this->nx=other.nx;
        this->ny=other.ny;
        this->x=new type*[nx];
        this->y=new type*[nx];
        this->f=new type*[nx];
        for(int i=0;i<nx;i++)
        {
            this->x[i]=new type[ny];
            this->y[i]=new type[ny];
            this->f[i]=new type[ny];
            for(int j=0;j<ny;j++)
            {
                this->x[i][j]=other.x[i][j];
                this->y[i][j]=other.y[i][j];
                this->f[i][j]=other.f[i][j];
            }
        }
        return *this;
    }

    template<typename type>
    bool function<type>::operator == (const function<type> & other) const
    {
        #ifdef domaincheck
            if(domainIsEqual(other))
            {
                for(int i=0;i<nx;i++)
                {
                    for(int j=0;j<ny;j++)
                    {
                        if(!(this->f[i][j]==other.f[i][j])){return false;}
                    }
                }
                return true;
            }
            else
            {
                return false;
            }
        #else
            for(int i=0;i<nx;i++)
            {
                for(int j=0;j<ny;j++)
                {
                    if(!(this->f[i][j]==other.f[i][j])){return false;}
                }
            }
            return true;
        #endif
    }

    template<typename type>
    bool function<type>::operator != (const function<type> & other) const
    {
        #ifdef domaincheck
            if(domainIsEqual(other))
            {
                for(int i=0;i<nx;i++)
                {
                    for(int j=0;j<ny;j++)
                    {
                        if(!(this->f[i][j]==other.f[i][j])){return true;}
                    }
                }
                return false;
            }
            else
            {
                return true;
            }
        #else
            for(int i=0;i<nx;i++)
            {
                for(int j=0;j<ny;j++)
                {
                    if(!(this->f[i][j]==other.f[i][j])){return true;}
                }
            }
            return false;
        #endif
    }

    template<typename type>
    function<type> function<type>::operator +(const function &other) const
    {
        function<type> res(*this);
        #ifdef domaincheck
            if(domainIsEqual(other)){
                for(int i=0;i<nx;i++)
                {
                    for(int j=0;j<ny;j++)
                    {
                        res.f[i][j]+=other.f[i][j];
                    }
                }
            }
        #else
            for(int i=0;i<nx;i++)
            {
                for(int j=0;j<ny;j++)
                {
                    res.f[i][j]+=other.f[i][j];
                }
            }
        #endif
        return res;
    }

    template<typename type>
    function<type> function<type>::operator -(const function &other) const
    {
        function<type> res(*this);
        #ifdef domaincheck
            if(domainIsEqual(other)){
                for(int i=0;i<nx;i++)
                {
                    for(int j=0;j<ny;j++)
                    {
                        res.f[i][j]-=other.f[i][j];
                    }
                }
            }
        #else
            for(int i=0;i<nx;i++)
            {
                for(int j=0;j<ny;j++)
                {
                    res.f[i][j]-=other.f[i][j];
                }
            }
        #endif
        return res;
    }

    template<typename type>
    function<type> function<type>::operator +(const type &other) const
    {
        function<type> res(*this);
        for(int i=0;i<nx;i++)
        {
            for(int j=0;j<ny;j++)
            {
                res.f[i][j]+=other;
            }
        }
        return res;
    }

    template<typename type>
    function<type> function<type>::operator -(const type &other) const
    {
        function<type> res(*this);
        for(int i=0;i<nx;i++)
        {
            for(int j=0;j<ny;j++)
            {
                res.f[i][j]-=other;
            }
        }
        return res;
    }

    template<typename type>
    function<type> function<type>::operator -() const
    {
        function<type> res(*this);
        for(int i=0;i<nx;i++)
        {
            for(int j=0;j<ny;j++)
            {
                res.f[i][j]*=(-1);
            }
        }
        return res;
    }

    template<typename type>
    function<type> function<type>::operator *(const function &other) const
    {
        function<type> res(*this);
        #ifdef domaincheck
            if(domainIsEqual(other)){
                for(int i=0;i<nx;i++)
                {
                    for(int j=0;j<ny;j++)
                    {
                        res.f[i][j]*=other.f[i][j];
                    }
                }
            }
        #else
                for(int i=0;i<nx;i++)
                {
                    for(int j=0;j<ny;j++)
                    {
                        res.f[i][j]*=other.f[i][j];
                    }
                }
        #endif
        return res;
    }

    template<typename type>
    function<type> function<type>::operator /(const function &other) const
    {
        function<type> res(*this);
        #ifdef domaincheck
            if(domainIsEqual(other)&&(zeropoint()==-1)){
                for(int i=0;i<nx;i++)
                {
                    for(int j=0;j<ny;j++)
                    {
                        res.f[i][j]/=other.f[i][j];
                    }
                }
            }
        #else
            for(int i=0;i<nx;i++)
            {
                for(int j=0;j<ny;j++)
                {
                    res.f[i][j]/=other.f[i][j];
                }
            }
        #endif
        return res;
    }

    template<typename type>
    function<type> function<type>::operator *(const type &C) const
    {
        function<type> res(*this);
        for(int i=0;i<nx;i++)
        {
            for(int j=0;j<ny;j++)
            {
                res.f[i][j]*=C;
            }
        }
        return res;
    }

    template<typename type>
    function<type> function<type>::operator /(const type & C) const
    {
        function<type> res(*this);
        for(int i=0;i<nx;i++)
        {
            for(int j=0;j<ny;j++)
            {
                res.f[i][j]/=C;
            }
        }
        return res;
    }

    template<typename T1,typename T2>
    function<T2> operator + (const T1 & C,const function<T2> & other)
    {
        function<T2> res(other);
        for(int i=0;i<other.nx;i++)
        {
            for(int j=0;j<other.ny;j++)
            {
                res.f[i][j] = C;
                res.f[i][j] += other.f[i][j];
            }
        }
        return res;
    }

    template<typename T1,typename T2>
    function<T2> operator - (const T1 & C,const function<T2> & other)
    {
        function<T2> res(other);
        for(int i=0;i<other.nx;i++)
        {
            for(int j=0;j<other.ny;j++)
            {
                res.f[i][j] = C;
                res.f[i][j] -= other.f[i][j];
            }
        }
        return res;
    }

    template<typename T1,typename T2>
    function<T2> operator * (const T1 & C,const function<T2> & other)
    {
        function<T2> res(other);
        for(int i=0;i<other.nx;i++)
        {
            for(int j=0;j<other.ny;j++)
            {
                res.f[i][j] = C;
                res.f[i][j] *= other.f[i][j];
            }
        }
        return res;
    }

    template<typename T1,typename T2>
    function<T2> operator / (const T1 & C,const function<T2> & other)
    {
        function<T2> res(other);
        for(int i=0;i<other.nx;i++)
        {
            for(int j=0;j<other.ny;j++)
            {
                res.f[i][j] = C;
                res.f[i][j] /= other.f[i][j];
            }
        }
        return res;
    }

    template<typename type>
    type & function<type>::operator ()(const int i, const int j)
    {
        return f[i][j];
    }

    template<typename T>
    std::ostream & operator<<(std::ostream &out, const XY::function<T> &point)
    {
        for(int i=0;i<point.nx;i++)
        {
            for(int j=0;j<point.ny;j++)
            {
                out <<std::setw(point.setw)<<std::fixed<<std::setprecision(point.setprecision)<<point.f[i][j]<<"  ";
            }
            out<<std::endl;
        }
        return out;
    }

    template<typename T>
    function<T> exp(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            for(int j=0;j<other.ny;j++)
            {
                res.setf(i,j,
                         std::exp(other.getf(i,j)));
            }
        }
        return res;
    }

    template<typename T>
    function<T> ln(const function<T> & other)
    {
        function<T> res(other);
        #ifdef domaincheck
            if(res.positive(domainState::positive))
            {
                for(int i=0;i<other.nx;i++)
                {
                    for(int j=0;j<other.ny;j++)
                    {
                        res.setf(i,j,
                                 std::log(other.getf(i,j)));
                    }
                }
            }else
            {
                std::cout<<"negative or zero argument, was returned argument"<<std::endl;
            }
        #else
            for(int i=0;i<other.nx;i++)
            {
                for(int j=0;j<other.ny;j++)
                {
                    res.setf(i,j,
                             std::log(other.getf(i,j)));
                }
            }
        #endif
        return res;
    }

    template<typename T>
    function<T> sqrt(const function<T> & other)
    {
        function<T> res(other);
        #ifdef domaincheck
            if(res.positive(domainState::positive))
            {
                for(int i=0;i<other.nx;i++)
                {
                    for(int j=0;j<other.ny;j++)
                    {
                        res.setf(i,j,
                                 std::sqrt(other.getf(i,j)));
                    }
                }
            }else
            {
                std::cout<<"negative or zero argument, was returned argument"<<std::endl;
            }
        #else
            for(int i=0;i<other.nx;i++)
            {
                for(int j=0;j<other.ny;j++)
                {
                    res.setf(i,j,
                             std::sqrt(other.getf(i,j)));
                }
            }
        #endif
        return res;
    }

    template<typename T>
    function<T> pow(const function<T> & other, const T & deg)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            for(int j=0;j<other.ny;j++)
            {
                res.setf(i,j,
                         std::pow(other.getf(i,j),deg));
            }
        }
        return res;
    }

    template<typename T>
    function<T> sin(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            for(int j=0;j<other.ny;j++)
            {
                res.setf(i,j,
                         std::sin(other.getf(i,j)));
            }
        }
        return res;
    }

    template<typename T>
    function<T> cos(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            for(int j=0;j<other.ny;j++)
            {
                res.setf(i,j,
                         std::cos(other.getf(i,j)));
            }
        }
        return res;
    }

    template<typename T>
    function<T> tg(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            for(int j=0;j<other.ny;j++)
            {
                res.setf(i,j,
                         std::tan(other.getf(i,j)));
            }
        }
        return res;
    }

    template<typename T>
    function<T> ctg(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            for(int j=0;j<other.ny;j++)
            {
                res.setf(i,j,
                         1.0/std::tan(other.getf(i,j)));
            }
        }
        return res;
    }

    template<typename T>
    function<T> sinh(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            for(int j=0;j<other.ny;j++)
            {
                res.setf(i,j,
                         std::sinh(other.getf(i,j)));
            }
        }
        return res;
    }

    template<typename T>
    function<T> cosh(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            for(int j=0;j<other.ny;j++)
            {
                res.setf(i,j,
                         std::cosh(other.getf(i,j)));
            }
        }
        return res;
    }

    template<typename T>
    function<T> tgh(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            for(int j=0;j<other.ny;j++)
            {
                res.setf(i,j,
                         std::tanh(other.getf(i,j)));
            }
        }
        return res;
    }

    template<typename T>
    function<T> ctgh(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            for(int j=0;j<other.ny;j++)
            {
                res.setf(i,j,
                         1.0/std::tanh(other.getf(i,j)));
            }
        }
        return res;
    }

    template<typename T>
    function<T> arcsin(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            for(int j=0;j<other.ny;j++)
            {
                res.setf(i,j,
                         std::asin(other.getf(i,j)));
            }
        }
        return res;
    }

    template<typename T>
    function<T> arccos(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            for(int j=0;j<other.ny;j++)
            {
                res.setf(i,j,
                         std::acos(other.getf(i,j)));
            }
        }
        return res;
    }

    template<typename T>
    function<T> arctg(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            for(int j=0;j<other.ny;j++)
            {
                res.setf(i,j,
                         std::atan(other.getf(i,j)));
            }
        }
        return res;
    }

    template<typename T>
    function<T> arcctg(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            for(int j=0;j<other.ny;j++)
            {
                res.setf(i,j,
                         1.0/std::atan(other.getf(i,j)));
            }
        }
        return res;
    }

    template<typename T>
    function<T> arcsinh(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            for(int j=0;j<other.ny;j++)
            {
                res.setf(i,j,
                         std::asinh(other.getf(i,j)));
            }
        }
        return res;
    }

    template<typename T>
    function<T> arccosh(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            for(int j=0;j<other.ny;j++)
            {
                res.setf(i,j,
                         std::acosh(other.getf(i,j)));
            }
        }
        return res;
    }

    template<typename T>
    function<T> arctgh(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            for(int j=0;j<other.ny;j++)
            {
                res.setf(i,j,
                         std::atanh(other.getf(i,j)));
            }
        }
        return res;
    }

    template<typename T>
    function<T> arcctgh(const function<T> & other)
    {
        function<T> res(other);
        for(int i=0;i<other.nx;i++)
        {
            for(int j=0;j<other.ny;j++)
            {
                res.setf(i,j,
                         1.0/std::atanh(other.getf(i,j)));
            }
        }
        return res;
    }

    template<typename type>
    type function<type>::integralS11(int x1, int x2, int y1, int y2) const
    {
        type *X=new type[x2-x1+1];
        for(int i=x1;i<x2-x1+1;i++)
        {
            X[i]=this->x[i][0];
        }
        X::function<type> fx(X,x2-x1+1);
        delete [] X;
        for(int i=x1;i<x2-x1+1;i++)
        {
            X::function<type> fy(y[i],f[i],y2-y1+1);
            fx(i)=fy.integralS1(y1,y2);
        }
        return fx.integralS1(x1,x2);
    }

    template<typename type>
    type function<type>::integralS22(int x1, int x2, int y1, int y2) const
    {
        type *X=new type[x2-x1+1];
        for(int i=x1;i<x2-x1+1;i++)
        {
            X[i]=this->x[i][0];
        }
        X::function<type> fx(X,x2-x1+1);
        delete [] X;
        for(int i=x1;i<x2-x1+1;i++)
        {
            X::function<type> fy(y[i],f[i],y2-y1+1);
            fx(i)=fy.integralS2(y1,y2);
        }
        return fx.integralS2(x1,x2);
    }

    template<typename type>
    type function<type>::integralS12(int x1, int x2, int y1, int y2) const
    {
        type *X=new type[x2-x1+1];
        for(int i=x1;i<x2-x1+1;i++)
        {
            X[i]=this->x[i][0];
        }
        X::function<type> fx(X,x2-x1+1);
        delete [] X;
        for(int i=x1;i<x2-x1+1;i++)
        {
            X::function<type> fy(y[i],f[i],y2-y1+1);
            fx(i)=fy.integralS2(y1,y2);
        }
        return fx.integralS1(x1,x2);
    }

    template<typename type>
    type function<type>::integralS21(int x1, int x2, int y1, int y2) const
    {
        type *X=new type[x2-x1+1];
        for(int i=x1;i<x2-x1+1;i++)
        {
            X[i]=this->x[i][0];
        }
        X::function<type> fx(X,x2-x1+1);
        delete [] X;
        for(int i=x1;i<x2-x1+1;i++)
        {
            X::function<type> fy(y[i],f[i],y2-y1+1);
            fx(i)=fy.integralS1(y1,y2);
        }
        return fx.integralS2(x1,x2);
    }

    template<typename T>
    T integralS22(int x1, int x2, int y1, int y2, const function<T> & other)
    {
        T res = other.integralS22(x1,x2,y1,y2);
        return res;
    }

    template<typename T>
    T integralS11(int x1, int x2, int y1, int y2, const function<T> & other)
    {
        T res = other.integralS11(x1,x2,y1,y2);
        return res;
    }

    template<typename T>
    T integralS21(int x1, int x2, int y1, int y2, const function<T> & other)
    {
        T res = other.integralS21(x1,x2,y1,y2);
        return res;
    }

    template<typename T>
    T integralS12(int x1, int x2, int y1, int y2, const function<T> & other)
    {
        T res = other.integralS12(x1,x2,y1,y2);
        return res;
    }

    template<typename type>
    function<type> & function<type>::spline(int Nx, int Ny)
    {
        int nxnew = Nx*(nx-1)+1;
        function<type> buf(nxnew,ny);

        for(int j=0;j<ny;j++)
        {
            type *x1 = new type[nx];
            type *f1 = new type[nx];
            for(int i=0;i<nx;i++)
            {
                x1[i]=this->x[i][j];
                f1[i]=this->f[i][j];
            }
            X::function<type> F(x1,f1,nx);
            delete [] x1; delete [] f1;
            F = X::spline(F,Nx);

            for(int i=0;i<nxnew;i++)
            {
                buf.setx(i,j,F.getx(i));
                buf.sety(i,j,this->y[0][j]);
                buf.setf(i,j,F.getf(i));
            }
        }

        int nynew = Ny*(ny-1)+1;
        function<type> res(nxnew,nynew);

        for(int j=0;j<nxnew;j++)
        {
            type *y1 = new type[ny];
            type *f1 = new type[ny];
            for(int i=0;i<ny;i++)
            {
                y1[i]=buf.gety(j,i);
                f1[i]=buf.getf(j,i);
            }
            X::function<type> F(y1,f1,ny);
            delete [] y1; delete [] f1;
            F = X::spline(F,Ny);

            for(int i=0;i<nynew;i++)
            {
                res.setx(j,i, buf.getx(j,0));
                res.sety(j,i, F.getx(i));
                res.setf(j,i, F.getf(i));
            }
        }
        *this = res;
        return *this;
    }

    template<typename T>
    function<T> spline(const function<T> & f, int Nx, int Ny)
    {
        int nx = Nx*(f.getnx()-1)+1;
        int ny = f.getny();
        function<T> buf(nx,ny);

        for(int j=0;j<f.getny();j++)
        {
            T *x = new T[f.getnx()];
            T *f1 = new T[f.getnx()];
            for(int i=0;i<f.getnx();i++)
            {
                x[i]=f.getx(i,j);
                f1[i]=f.getf(i,j);
            }
            X::function<T> F(x,f1,f.getnx());
            delete [] x; delete [] f1;
            F = X::spline(F,Nx);

            for(int i=0;i<nx;i++)
            {
                buf.setx(i,j,F.getx(i));
                buf.sety(i,j,f.gety(0,j));
                buf.setf(i,j,F.getf(i));
            }
        }

        ny = Ny*(f.getny()-1)+1;
        function<T> res(nx,ny);

        for(int j=0;j<nx;j++)
        {
            T *y = new T[f.getny()];
            T *f1 = new T[f.getny()];
            for(int i=0;i<f.getny();i++)
            {
                y[i]=buf.gety(j,i);
                f1[i]=buf.getf(j,i);
            }
            X::function<T> F(y,f1,f.getny());
            delete [] y; delete [] f1;
            F = X::spline(F,Ny);

            for(int i=0;i<ny;i++)
            {
                res.setx(j,i, buf.getx(j,0));
                res.sety(j,i, F.getx(i));
                res.setf(j,i, F.getf(i));
            }
        }
        return res;
    }
    template<typename type>
    point<type> function<type>::maxPoint() const
    {
        point<type> max(x[0][0],y[0][0],f[0][0],0,0);

        for(int i=0;i<nx;i++)
        {
            for(int j=0;j<ny;j++)
            {
                if(f[i][j]>max.f)
                {
                    max.f=f[i][j];
                    max.x=x[i][j];
                    max.x=y[i][j];
                    max.ix=i;
                    max.iy=j;
                }
            }
        }
        return max;
    }

    template<typename type>
    point<type> function<type>::minPoint() const
    {
        point<type> min(x[0][0],y[0][0],f[0][0],0,0);

        for(int i=0;i<nx;i++)
        {
            for(int j=0;j<ny;j++)
            {
                if(f[i][j]<min.f)
                {
                    min.f=f[i][j];
                    min.x=x[i][j];
                    min.x=y[i][j];
                    min.ix=i;
                    min.iy=j;
                }
            }
        }
        return min;
    }

    template<typename type>
    const char *function<type>::savef(const char *name) const
    {
        FILE *fp;
        fp=freopen(name,"w",stdout);
        for(int i=0;i<nx;i++)
        {
            for(int j=0;j<ny;j++)
            {
                std::cout <<std::setw(setw)<<std::fixed<<std::setprecision(setprecision)<<f[i][j]<<"  ";
            }
            std::cout<<std::endl;
        }
        fclose(fp);
        freopen("CON","w",stdout);
        return name;
    }

    template<typename type>
    const char *function<type>::savex(const char *name) const
    {
        FILE *fp;
        fp=freopen(name,"w",stdout);
        for(int i=0;i<nx;i++)
        {
            for(int j=0;j<ny;j++)
            {
                std::cout <<std::setw(setw)<<std::fixed<<std::setprecision(setprecision)<<x[i][j]<<"  ";
            }
            std::cout<<std::endl;
        }
        fclose(fp);
        freopen("CON","w",stdout);
        return name;
    }

    template<typename type>
    const char *function<type>::savey(const char *name) const
    {
        FILE *fp;
        fp=freopen(name,"w",stdout);
        for(int i=0;i<nx;i++)
        {
            for(int j=0;j<ny;j++)
            {
                std::cout <<std::setw(setw)<<std::fixed<<std::setprecision(setprecision)<<y[i][j]<<"  ";
            }
            std::cout<<std::endl;
        }
        fclose(fp);
        freopen("CON","w",stdout);
        return name;
    }
}

#endif
