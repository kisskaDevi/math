#ifndef _PSIFUNCTION_H_
#define _PSIFUNCTION_H_

#include "function.h"

int psi_function_spline = 5;

namespace psi {

    template<typename type>
    class function;

    template<typename T> std::ostream & operator<<(std::ostream & out, const function<T> & point);
    template<typename T> function<T> exp(const function<T> & other);
    template<typename T> function<T> exp(const XY::function<cn::complex<T>> & other);
    template<typename T> function<T> sinc(const function<T> & other);
    template<typename T> function<T> sinc(const XY::function<cn::complex<T>> & other);

    template<typename T> function<T> con(function<T> other);

    template<typename T> cn::complex<T> integralS11(int x0, int xn, int y0, int yn, const psi::function<T> & other);
    template<typename T> cn::complex<T> integralS22(int x0, int xn, int y0, int yn,const psi::function<T> & other);

    template<typename T> T width(const X::function<T> & drho);

    template<typename type>
    class function
    {
    private:
        XY::function<cn::complex<type>> f;
        int setprecision=1, setw=2;

        type Schmidt1() const;
        type Schmidt2() const;
        type SchmidtRho1() const;
        type SchmidtRho2() const;
        type R1() const;
        type R2() const;
        type R3() const;
    public:
        function();
        function(const XY::function<cn::complex<type>> & f);
        function(const function<type> & f);
        ~function();
        cn::complex<type> getf(int i, int j) const;
        cn::complex<type> getx(int i, int j) const;
        cn::complex<type> gety(int i, int j) const;
        int getnx() const;
        int getny() const;
        void getx() const;
        void gety() const;
        void setf(int i, int j, const cn::complex<type> & f);
        void setx(int i, int j, const cn::complex<type> & x);
        void sety(int i, int j, const cn::complex<type> & y);
        void setout(int setprecision, int setw);
        bool operator == (const function<type> & other);
        bool operator != (const function<type> & other);
        function &operator = (const function<type> & other);
        function &operator = (const XY::function<cn::complex<type>> & other);
        function operator + (const function<type> & other);
        function operator - (const function<type> & other);
        function operator + (const cn::complex<type> & other);
        function operator - (const cn::complex<type> & other);
        function operator - ();
        function operator * (const function<type> & other);
        function operator / (const function<type> & other);
        function operator * (const cn::complex<type> & C);
        function operator / (const cn::complex<type> & C);
        cn::complex<type> & operator () (const int i, const int j);
        template<typename T> friend std::ostream & operator<<(std::ostream & out, const function<T> & point);
        template<typename T> friend function<T> exp(const function<T> & other);
        template<typename T> friend function<T> exp(const XY::function<cn::complex<T>> & other);
        template<typename T> friend function<T> sinc(const function<T> & other);
        template<typename T> friend function<T> sinc(const XY::function<cn::complex<T>> & other);

        template<typename T> friend function<T> con(function<T> other);

        cn::complex<type> integralS11(int x0, int xn, int y0, int yn);
        cn::complex<type> integralS22(int x0, int xn, int y0, int yn);
        template<typename T> friend cn::complex<T> integralS11(int x0, int xn, int y0, int yn, const psi::function<T> & other);
        template<typename T> friend cn::complex<T> integralS22(int x0, int xn, int y0, int yn,const psi::function<T> & other);

        void normalize();
        function<type> Fourier(const XY::function<cn::complex<type>> & x) const;
        function<type> rho() const;

        void normalizeS11();
        function<type> FourierS11(const XY::function<cn::complex<type>> & x) const;
        function<type> rhoS1() const;

        XY::function<type> Re() const;
        XY::function<type> Im() const;

        type entanglementSchmidt(int setType) const;
        type entanglementSchmidtRho(int setType) const;
        type entanglementR(int setType) const;
        template<typename T> friend T width(const X::function<T> & drho);

        X::function<type> majorDiagonalRho(const XY::function<cn::complex<type>> & x) const;
        X::function<type> minorDiagonalRho(const XY::function<cn::complex<type>> & x) const;

        X::function<type> majorTrDiagonalPsi() const;
        X::function<type> minorTrDiagonalPsi() const;

        type majorDiagonalWidth() const;
        type minorDiagonalWidth() const;

        X::function<type> singlePsi() const;
        X::function<type> conditionalPsi() const;

        void spline(int Nx, int Ny);

    };

}

//============================================================================================================================================//
//============================================================================================================================================//
//============================================================================================================================================//

namespace psi{

    template<typename type>
    function<type>::function()
    {}

    template<typename type>
    function<type>::function(const XY::function<cn::complex<type>> & f)
    {
        this->f=f;
    }

    template<typename type>
    function<type>::function(const function<type> &f)
    {
        this->f=f.f;
    }

    template<typename type>
    function<type>::~function()
    {}

    template<typename type>
    cn::complex<type> function<type>::getf(int i, int j) const
    {
        return f.getf(i,j);
    }

    template<typename type>
    cn::complex<type> function<type>::getx(int i, int j) const
    {
        return f.getx(i,j);
    }

    template<typename type>
    cn::complex<type> function<type>::gety(int i, int j) const
    {
        return f.gety(i,j);
    }

    template<typename type>
    int function<type>::getnx() const
    {
        return f.getnx();
    }

    template<typename type>
    int function<type>::getny() const
    {
        return f.getny();
    }

    template<typename type>
    void function<type>::getx() const
    {
        for(int i=0;i<f.getnx();i++)
        {
            for(int j=0;j<f.getny();j++)
            {
                std::cout<<f.getx(i,j)<<"   ";
            }
            std::cout<<std::endl;
        }
    }

    template<typename type>
    void function<type>::gety() const
    {
        for(int i=0;i<f.getnx();i++)
        {
            for(int j=0;j<f.getny();j++)
            {
                std::cout<<f.gety(i,j)<<"   ";
            }
            std::cout<<std::endl;
        }
    }

    template<typename type>
    void function<type>::setf(int i, int j, const cn::complex<type> & f)
    {
        this->f(i,j)=f;
    }

    template<typename type>
    void function<type>::setx(int i, int j, const cn::complex<type> & x)
    {
        f.setx(i,j,x);
    }

    template<typename type>
    void function<type>::sety(int i, int j, const cn::complex<type> & y)
    {
        f.sety(i,j,y);
    }

    template<typename type>
    void function<type>::setout(int setprecision, int setw)
    {
        this->setprecision=setprecision;
        this->setw=setw;
    }

    template<typename type>
    bool function<type>::operator == (const function<type> & other)
    {
        return *this==other;
    }

    template<typename type>
    bool function<type>::operator != (const function<type> & other)
    {
        return *this!=other;
    }

    template<typename type>
    function<type> &function<type>::operator = (const function<type> & other)
    {
        f=other.f;
        return *this;
    }

    template<typename type>
    function<type> & function<type>::operator = (const XY::function<cn::complex<type>> & other)
    {
        if(f.domainIsEqual(other))
        {
            f=other;
        }
        return *this;
    }

    template<typename type>
    function<type> function<type>::operator + (const function<type> & other)
    {
        function<type> res(*this);
        res.f = res.f + other.f;
        return res;
    }

    template<typename type>
    function<type> function<type>::operator - (const function<type> & other)
    {
        function<type> res(*this);
        res.f = res.f - other.f;
        return res;
    }

    template<typename type>
    function<type> function<type>::operator + (const cn::complex<type> & other)
    {
        function<type> res(*this);
        res.f = res.f + other;
        return res;
    }

    template<typename type>
    function<type> function<type>::operator - (const cn::complex<type> & other)
    {
        function<type> res(*this);
        res.f = res.f - other;
        return res;
    }

    template<typename type>
    function<type> function<type>::operator - ()
    {
        function<type> res(*this);
        res.f = res.f*(-1);
        return res;
    }

    template<typename type>
    function<type> function<type>::operator * (const function<type> & other)
    {
        function<type> res(*this);
        res.f = res.f * other.f;
        return res;
    }

    template<typename type>
    function<type> function<type>::operator / (const function<type> & other)
    {
        function<type> res(*this);
        res.f = res.f / other.f;
        return res;
    }

    template<typename type>
    function<type> function<type>::operator * (const cn::complex<type> & C)
    {
        function<type> res(*this);
        res.f = res.f * C;
        return res;
    }

    template<typename type>
    function<type> function<type>::operator / (const cn::complex<type> & C)
    {
        function<type> res(*this);
        res.f = res.f / C;
        return res;
    }

    template<typename type>
    cn::complex<type> & function<type>::operator () (const int i, const int j)
    {
        return f(i,j);
    }

    template<typename T>
    std::ostream & operator<<(std::ostream & out, const function<T> & point)
    {
        out<<point.f;
        return out;
    }

    template<typename T>
    function<T> exp(const function<T> & other)
    {
        function<T> res(other);

        for(int i=0;i<other.f.getnx();i++)
        {
            for(int j=0;j<other.f.getny();j++)
            {
                res(i,j) = cn::exp(other.f(i,j));
            }
        }
        return res;
    }

    template<typename T>
    function<T> exp(const XY::function<cn::complex<T>> & other)
    {
        function<T> res(other);

        for(int i=0;i<other.getnx();i++)
        {
            for(int j=0;j<other.getny();j++)
            {
                res(i,j) = cn::exp(other.getf(i,j));
            }
        }
        return res;
    }

    template<typename T>
    function<T> sinc(const function<T> & other)
    {
        function<T> res(other);
        cn::complex<T> one(1.0,0.0);
        cn::complex<T> zero(0.0,0.0);

        for(int i=0;i<other.f.getnx();i++)
        {
            for(int j=0;j<other.f.getny();j++)
            {
                if(other.f(i,j)==zero)
                {
                    res(i,j)=one;
                }
                else
                {
                    res(i,j)=cn::sin(other.f(i,j))/(other.f(i,j));
                }
            }
        }

        return res;
    }

    template<typename T>
    function<T> sinc(const XY::function<cn::complex<T>> & other)
    {
        function<T> res(other);
        cn::complex<T> one(1.0,0.0);
        cn::complex<T> zero(0.0,0.0);

        for(int i=0;i<other.getnx();i++)
        {
            for(int j=0;j<other.getny();j++)
            {
                if(other.getf(i,j)==zero)
                {
                    res(i,j)=one;
                }
                else
                {
                    res(i,j)=cn::sin(other.getf(i,j))/(other.getf(i,j));
                }
            }
        }

        return res;
    }

    template<typename T>
    function<T> con(function<T> other)
    {
        function<T> res(other);

        for(int i=0;i<other.f.getnx();i++)
        {
            for(int j=0;j<other.f.getny();j++)
            {
                res(i,j) = cn::con(other.f(i,j));
            }
        }
        return res;
    }

    template<typename type>
    cn::complex<type> function<type>::integralS11(int x0, int xn, int y0, int yn)
    {
        return f.integralS11(x0,xn,y0,yn);
    }

    template<typename type>
    cn::complex<type> function<type>::integralS22(int x0, int xn, int y0, int yn)
    {
        return f.integralS22(x0,xn,y0,yn);
    }

    template<typename type>
    void function<type>::normalize()
    {
        psi::function<type> buf(f);
        buf=buf*psi::con(buf);
        cn::complex<type> N = buf.integralS22(0,f.getnx()-1,0,f.getny()-1);
        N.set(std::sqrt(N.Real()),0.0);
        f=f/N;
    }

    template<typename type>
    void function<type>::normalizeS11()
    {
        psi::function<type> buf(f);
        buf=buf*psi::con(buf);
        cn::complex<type> N = buf.integralS11(0,f.getnx()-1,0,f.getny()-1);
        N.set(std::sqrt(N.Real()),0.0);
        f=f/N;
    }

    template<typename T>
    cn::complex<T> integralS11(int x0, int xn, int y0, int yn, psi::function<T> other)
    {
        return other.integralS11(x0,xn,y0,yn);
    }

    template<typename T>
    cn::complex<T> integralS22(int x0, int xn, int y0, int yn, psi::function<T> other)
    {
        return other.integralS22(x0,xn,y0,yn);
    }

    template<typename type>
    function<type> function<type>::Fourier(const XY::function<cn::complex<type>> &x) const
    {
        function<type> res(x);
        int nx=x.getnx();
        int ny=x.getny();
        int kx=f.getnx();
        int ky=f.getny();
        XY::function<cn::complex<type>> F(f);
        for(int i=0;i<nx;i++)
        {
            for(int j=0;j<ny;j++)
            {
                for(int i1=0;i1<kx;i1++)
                {
                    for(int j1=0;j1<ky;j1++)
                    {
                        cn::complex<type> ione(0.0,1.0);
                        cn::complex<type> arg1 = x.getx(i,j); arg1 *= f.getx(i1,j1);
                        cn::complex<type> arg2 = x.gety(i,j); arg2 *= f.gety(i1,j1);
                        arg1 += arg2;
                        ione *= arg1;
                        F(i1,j1) = f.getf(i1,j1);
                        F(i1,j1) *= ione.exp();
                        //F(i1,j1) = f.getf(i1,j1)*cn::exp(ione*(x.getx(i,j)*f.getx(i1,j1)+x.gety(i,j)*f.gety(i1,j1)));
                    }
                }
                res(i,j) = F.integralS22(0,kx-1,0,ky-1);
            }
        }
        return res;
    }

    template<typename type>
    function<type> function<type>::FourierS11(const XY::function<cn::complex<type>> &x) const
    {
        function<type> res(x);
        int nx=x.getnx();
        int ny=x.getny();
        int kx=f.getnx();
        int ky=f.getny();
        cn::complex<type> ione(0.0,1.0);
        XY::function<cn::complex<type>> F(f);
        for(int i=0;i<nx;i++)
        {
            for(int j=0;j<ny;j++)
            {
                for(int i1=0;i1<kx;i1++)
                {
                    for(int j1=0;j1<ky;j1++)
                    {
                        F(i1,j1) = f(i1,j1)*cn::exp(ione*(x.getx(i,j)*f.getx(i1,j1)+x.gety(i,j)*f.gety(i1,j1)));
                    }
                }
                res(i,j) = F.integralS11(0,kx-1,0,ky-1);
            }
        }
        return res;
    }

    template<typename type>
    function<type> function<type>::rho() const
    {
        cn::complex<type> **x1 = new cn::complex<type>*[f.getnx()];
        cn::complex<type> **x2 = new cn::complex<type>*[f.getnx()];
        cn::complex<type> *y1 = new cn::complex<type>[f.getny()];
        for(int i=0;i<f.getnx();i++)
        {
            x1[i]=new cn::complex<type>[f.getnx()];
            x2[i]=new cn::complex<type>[f.getnx()];
            for(int j=0;j<f.getnx();j++)
            {
                x1[i][j]=f.getx(i,0);
                x2[i][j]=f.getx(j,0);
            }
        }
        for(int i=0;i<f.getny();i++)
        {
            y1[i]=f.gety(0,i);
        }
        XY::function<cn::complex<type>> x(x1,x2,f.getnx(),f.getnx());
        X::function<cn::complex<type>> r(y1,f.getny());
        for(int i=0;i<f.getnx();i++)
        {
            delete [] x1[i];
            delete [] x2[i];
        }
        delete [] x1; delete [] x2; delete [] y1;

        function<type> res(x);
        int nx=f.getnx();
        int ny=f.getny();
        for(int i=0;i<nx;i++)
        {
            for(int j=0;j<nx;j++)
            {
                for(int k=0;k<ny;k++)
                {
                    r(k) = cn::con(f.getf(j,k))*f.getf(i,k);
                }
                res(i,j)= r.integralS2(0,ny-1);
            }
        }
        return res;
    }

    template<typename type>
    function<type> function<type>::rhoS1() const
    {
        cn::complex<type> **x1 = new cn::complex<type>*[f.getnx()];
        cn::complex<type> **x2 = new cn::complex<type>*[f.getnx()];
        cn::complex<type> *y1 = new cn::complex<type>[f.getny()];
        for(int i=0;i<f.getnx();i++)
        {
            x1[i]=new cn::complex<type>[f.getnx()];
            x2[i]=new cn::complex<type>[f.getnx()];
            for(int j=0;j<f.getnx();j++)
            {
                x1[i][j]=f.getx(i,0);
                x2[i][j]=f.getx(j,0);
            }
        }
        for(int i=0;i<f.getny();i++)
        {
            y1[i]=f.gety(0,i);
        }
        XY::function<cn::complex<type>> x(x1,x2,f.getnx(),f.getnx());
        X::function<cn::complex<type>> r(y1,f.getny());
        for(int i=0;i<f.getnx();i++)
        {
            delete [] x1[i];
            delete [] x2[i];
        }
        delete [] x1; delete [] x2; delete [] y1;

        function<type> res(x);
        int nx=f.getnx();
        int ny=f.getny();
        for(int i=0;i<nx;i++)
        {
            for(int j=0;j<nx;j++)
            {
                for(int k=0;k<ny;k++)
                {
                    r(k) = cn::con(f(j,k))*f(i,k);
                }
                res(i,j)= r.integralS1(0,ny-1);
            }
        }
        return res;
    }

    template<typename type>
    XY::function<type> function<type>::Re() const
    {
        type **x1 = new type*[f.getnx()];
        type **x2 = new type*[f.getnx()];
        for(int i=0;i<f.getnx();i++)
        {
            x1[i]=new type[f.getny()];
            x2[i]=new type[f.getny()];
        }
        XY::function<type> res(x1,x2,f.getnx(),f.getny());
        for(int i=0;i<f.getnx();i++)
        {
            delete [] x1[i];
            delete [] x2[i];
        }
        delete [] x1;
        delete [] x2;
        res.setout(setw,setprecision);
        for(int i=0;i<f.getnx();i++)
        {
            for(int j=0;j<f.getny();j++)
            {
                res.setx(i,j,f.getx(i,j).Real());
                res.sety(i,j,f.gety(i,j).Real());
                res(i,j)=f.getf(i,j).Real();
            }
        }
        return res;
    }

    template<typename type>
    XY::function<type> function<type>::Im() const
    {
        type **x1 = new type*[f.getnx()];
        type **x2 = new type*[f.getnx()];
        for(int i=0;i<f.getnx();i++)
        {
            x1[i]=new type[f.getny()];
            x2[i]=new type[f.getny()];
        }
        XY::function<type> res(x1,x2,f.getnx(),f.getny());
        for(int i=0;i<f.getnx();i++)
        {
            delete [] x1[i];
            delete [] x2[i];
        }
        delete [] x1;
        delete [] x2;
        res.setout(setw,setprecision);
        for(int i=0;i<f.getnx();i++)
        {
            for(int j=0;j<f.getny();j++)
            {
                res.setx(i,j,f.getx(i,j).Real());
                res.sety(i,j,f.gety(i,j).Real());
                res(i,j)=f.getf(i,j).Imaginary();
            }
        }
        return res;
    }

    template<typename type>
    type function<type>::Schmidt1() const
    {
        cn::complex<type> K(1.0,0.0);
        psi::function<type> psi(*this); psi.normalize();
        psi::function<type> rho = psi.rho();
        psi::function<type> Rho = rho*psi::con(rho);
        K /= Rho.integralS22(0,f.getnx()-1,0,f.getnx()-1);
        return K.Real();
    }

    template<typename type>
    type function<type>::Schmidt2() const
    {

        psi::function<type> psi(*this);
        psi.normalize();
        cn::complex<type> one(1.0,0.0);
        psi::function<type> rho = psi.rho();

        cn::complex<type> *x = new cn::complex<type>[f.getnx()];
        for(int i=0;i<f.getnx();i++)
        {
            x[i] = f.getx(i,0);
        }
        X::function<cn::complex<type>> Rho(x,f.getnx());
        psi::function<type> rhoBuff(rho);
        for(int i=0;i<f.getnx();i++)
        {
            for(int j=0;j<f.getnx();j++)
            {
                X::function<cn::complex<type>> buff(x,f.getnx());
                for(int k=0;k<f.getnx();k++)
                {
                   buff(k) = rho(i,k)*rho(k,j);
                }
                rhoBuff(i,j) = buff.integralS2(0,f.getnx()-1);
            }
        }
        for(int i=0;i<f.getnx();i++)
        {
            Rho(i) = rhoBuff(i,i);
        }

        cn::complex<type> K = one/Rho.integralS2(0,f.getnx()-1);

        delete [] x;

        return K.Real();
    }

    template<typename type>
    type function<type>::SchmidtRho1() const
    {
        psi::function<type> rho(*this);
        rho = rho*psi::con(rho);
        return 1.0/rho.integralS22(0,f.getnx()-1,0,f.getnx()-1).Real();
    }

    template<typename type>
    type function<type>::SchmidtRho2() const
    {
        psi::function<type> rho(*this);
        cn::complex<type> *x = new cn::complex<type>[f.getnx()];
        for(int i=0;i<f.getnx();i++)
        {
            x[i] = f.getx(i,0);
        }
        X::function<cn::complex<type>> buff(x,f.getnx());
        X::function<cn::complex<type>> Rho(x,f.getnx());
        delete [] x;

        for(int i=0;i<f.getnx();i++)
        {
            for(int k=0;k<f.getnx();k++)
            {
               buff(k) = rho(i,k)*rho(k,i);
            }
            Rho(i) = buff.integralS2(0,f.getnx()-1);
        }

        return 1.0/Rho.integralS2(0,f.getnx()-1).Real();
    }

    template<typename type>
    type function<type>::R1() const
    {
        type dplus;
        type dmin;

        type *x= new type[f.getnx()];
        type *F= new type[f.getnx()];
        type *_F= new type[f.getnx()];
        for(int i=0;i<f.getnx();i++)
        {
            x[i] = f.getx(i,i).Real();
            F[i] = f.getf(i,i).Real();
            _F[i] = f.getf(i,f.getnx()-1-i).Real();
        }
        X::function<type> drho(x,F,f.getnx());
        X::function<type> _drho(x,_F,f.getnx());
        delete [] x;
        delete [] F;
        delete [] _F;

    #ifdef psi_spline
        drho = X::spline(drho,psi_function_spline);
        _drho = X::spline(_drho,psi_function_spline);
    #endif

        type E = drho.average();
        X::point<type> Point = drho.maxPoint();
        type max = Point.f/2;
        int point = Point.ix;
        type eps = max;

        for(int i=0;i<drho.getn();i++)
        {
            if(std::abs(drho(i)-max)<eps)
            {
                eps=std::abs(drho(i)-max);
                point = i;
            }
        }

        int point2;
        if(drho(point)-max<0){ point2=point+1;}else{ point2=point-1;}
        type xcorss = (max-drho(point))/(drho(point2)-drho(point))*(drho.getx(point2)-drho.getx(point))+drho.getx(point);

    //    std::cout<<"E = "<<E<<std::endl;
    //    std::cout<<"max = "<<2*max<<std::endl;
    //    std::cout<<"point = "<<point<<std::endl;
    //    std::cout<<"point2 = "<<point2<<std::endl;
    //    std::cout<<"xcorss = "<<xcorss<<std::endl;

        dplus=std::abs(xcorss-E);

        E = _drho.average();
        Point = _drho.maxPoint();
        max = Point.f/2;
        point = Point.ix;
        eps = max;

        for(int i=0;i<_drho.getn();i++)
        {
            if(std::abs(_drho(i)-max)<eps)
            {
                eps=std::abs(_drho(i)-max);
                point = i;
            }
        }

        if(_drho(point)-max<0){ point2=point+1;}else{ point2=point-1;}
        xcorss = (max-_drho(point))/(_drho(point2)-_drho(point))*(_drho.getx(point2)-_drho.getx(point))+_drho.getx(point);

        //std::cout<<"E = "<<E<<std::endl;
        //std::cout<<"max = "<<2*max<<std::endl;
        //std::cout<<"point = "<<point<<std::endl;
        //std::cout<<"point2 = "<<point2<<std::endl;
        //std::cout<<"xcorss = "<<xcorss<<std::endl;

        dmin=std::abs(xcorss-E);

        return dplus/dmin;
    }

    template<typename type>
    type function<type>::R2() const
    {
        type *x= new type[f.getnx()];
        type *F= new type[f.getnx()];
        type *_F= new type[f.getnx()];
        for(int i=0;i<f.getnx();i++)
        {
            x[i] = f.getx(i,i).Real();
            F[i] = f.getf(i,i).Real();
            _F[i] = f.getf(i,f.getnx()-1-i).Real();
        }
        X::function<type> drho(x,F,f.getnx());
        X::function<type> _drho(x,_F,f.getnx());
        delete [] x; delete [] F; delete [] _F;

    #ifdef psi_spline
        drho = X::spline(drho,psi_function_spline);
        _drho = X::spline(_drho,psi_function_spline);
    #endif

        drho = drho*drho;
        _drho = _drho*_drho;

        drho.normalize();
        _drho.normalize();

        return std::sqrt(drho.dispersion()) / std::sqrt(_drho.dispersion());
    }

    template<typename type>
    type function<type>::R3() const
    {
        type *x= new type[f.getnx()];
        type *F= new type[f.getnx()];
        type *_F= new type[f.getnx()];
        for(int i=0;i<f.getnx();i++)
        {
            x[i] = f.getx(i,i).Real();
            F[i] = f.getf(i,i).Real();
            _F[i] = f.getf(i,f.getnx()-1-i).Real();
        }
        X::function<type> drho(x,F,f.getnx());
        X::function<type> _drho(x,_F,f.getnx());
        delete [] x; delete [] F; delete [] _F;

        drho.normalize();
        _drho.normalize();

        return std::sqrt(std::sqrt(drho.dispersion())) / std::sqrt(std::sqrt(_drho.dispersion()));
    }

    template<typename T>
    T width(const X::function<T> & drho)
    {
        T E = drho.average();
        X::point<T> Point = drho.maxPoint();
        T max = Point.f/2;
        int point = Point.ix;
        T eps = max;

        for(int i=0;i<drho.getn();i++)
        {
            if(std::abs(drho.getf(i)-max)<eps)
            {
                eps=std::abs(drho.getf(i)-max);
                point = i;
            }
        }

        int point2;
        if(drho.getf(point)-max<0){ point2=point+1;}else{ point2=point-1;}
        T xcorss = (max-drho.getf(point))/(drho.getf(point2)-drho.getf(point))*(drho.getx(point2)-drho.getx(point))+drho.getx(point);

        return std::abs(xcorss-E);
    }

    template<typename type>
    type function<type>::entanglementSchmidt(int setType) const
    {
        switch (setType) {
        case 1:
            return Schmidt1();
            break;
        case 2:
            return Schmidt2();
            break;
        }
    }

    template<typename type>
    type function<type>::entanglementSchmidtRho(int setType) const
    {
        switch (setType) {
        case 1:
            return SchmidtRho1();
            break;
        case 2:
            return SchmidtRho2();
            break;
        }
        return 0;
    }

    template<typename type>
    type function<type>::entanglementR(int setType) const
    {
        if(f.getnx()==f.getnx())
        {
            switch (setType) {
            case 1:
                return R1();
                break;
            case 2:
                return R2();
                break;
            case 3:
                return R3();
                break;
            }
        }
        return 0;
    }

    template<typename type>
    void function<type>::spline(int Nx, int Ny)
    {
        f = XY::spline(f,Nx,Ny);
    }

    template<typename type>
    X::function<type> function<type>::majorDiagonalRho(const XY::function<cn::complex<type>> & x) const
    {
        int nx=x.getnx();
        int kx=f.getnx();
        int ky=f.getny();
        type *xnew = new type[nx];
        for(int i=0;i<nx;i++)
        {
            xnew[i] = x.getx(i,i).Real();
        }
        X::function<type> res(xnew,nx);
        delete [] xnew;
        XY::function<cn::complex<type>> F(f);
        for(int i=0;i<nx;i++)
        {
            for(int i1=0;i1<kx;i1++)
            {
                for(int j1=0;j1<ky;j1++)
                {
                    cn::complex<type> ione(0.0,1.0);
                    cn::complex<type> arg1 = x.getx(i,i); arg1 *= f.getx(i1,j1);
                    cn::complex<type> arg2 = x.gety(i,i); arg2 *= f.gety(i1,j1);
                    arg1 -= arg2;
                    ione *= arg1;
                    F(i1,j1) = f.getf(i1,j1);
                    F(i1,j1) *= ione.exp();
                }
            }
            res(i) = F.integralS22(0,kx-1,0,ky-1).Real();
        }
        return res;
    }

    template<typename type>
    X::function<type> function<type>::minorDiagonalRho(const XY::function<cn::complex<type>> & x) const
    {
        int nx=x.getnx();
        int ny=nx;
        int kx=f.getnx();
        int ky=f.getny();
        type *xnew = new type[nx];
        for(int i=0;i<nx;i++)
        {
            xnew[i] = x.getx(i,ny-1-i).Real();
        }
        X::function<type> res(xnew,nx);
        delete [] xnew;
        XY::function<cn::complex<type>> F(f);
        for(int i=0;i<nx;i++)
        {
            for(int i1=0;i1<kx;i1++)
            {
                for(int j1=0;j1<ky;j1++)
                {
                    cn::complex<type> ione(0.0,1.0);
                    cn::complex<type> arg1 = x.getx(i,nx-1-i); arg1 *= f.getx(i1,j1);
                    cn::complex<type> arg2 = x.gety(i,ny-1-i); arg2 *= f.gety(i1,j1);
                    arg1 -= arg2;
                    ione *= arg1;
                    F(i1,j1) = f.getf(i1,j1);
                    F(i1,j1) *= ione.exp();
                }
            }
            res(i) = F.integralS22(0,kx-1,0,ky-1).Real();
        }
        return res;
    }

    template<typename type>
    type function<type>::majorDiagonalWidth() const
    {
        type *x= new type[f.getnx()];
        type *F= new type[f.getnx()];
        for(int i=0;i<f.getnx();i++)
        {
            x[i] = f.getx(i,i).Real();
            F[i] = f.getf(i,i).Real();
        }
        X::function<type> drho(x,F,f.getnx());
        delete [] x;
        delete [] F;

        type E = drho.average();
        X::point<type> Point = drho.maxPoint();
        type max = Point.f/2;
        int point = Point.ix;
        type eps = max;

        for(int i=0;i<drho.getn();i++)
        {
            if(std::abs(drho(i)-max)<eps)
            {
                eps=std::abs(drho(i)-max);
                point = i;
            }
        }

        int point2;
        if(drho(point)-max<0){ point2=point+1;}else{ point2=point-1;}
        type xcorss = (max-drho(point))/(drho(point2)-drho(point))*(drho.getx(point2)-drho.getx(point))+drho.getx(point);

        return std::abs(xcorss-E);
    }

    template<typename type>
    type function<type>::minorDiagonalWidth() const
    {
        type *x= new type[f.getnx()];
        type *_F= new type[f.getnx()];
        for(int i=0;i<f.getnx();i++)
        {
            x[i] = f.getx(i,i).Real();
            _F[i] = f.getf(i,f.getnx()-1-i).Real();
        }
        X::function<type> _drho(x,_F,f.getnx());
        delete [] x;
        delete [] _F;

        type E = _drho.average();
        X::point<type> Point = _drho.maxPoint();
        type max = Point.f/2;
        int point = Point.ix;
        type eps = max;

        for(int i=0;i<_drho.getn();i++)
        {
            if(std::abs(_drho(i)-max)<eps)
            {
                eps=std::abs(_drho(i)-max);
                point = i;
            }
        }

        int point2;
        if(_drho(point)-max<0){ point2=point+1;}else{ point2=point-1;}
        type xcorss = (max-_drho(point))/(_drho(point2)-_drho(point))*(_drho.getx(point2)-_drho.getx(point))+_drho.getx(point);

        return std::abs(xcorss-E);
    }

    template<typename type>
    X::function<type> function<type>::singlePsi() const
    {
        int kx=f.getnx();
        int ky=f.getny();
        int n0 = (ky-1)/2+1;
        type *xnew = new type[kx];
        for(int i=0;i<kx;i++)
        {
            xnew[i] = f.getx(i,n0).Real();
        }
        X::function<type> res(xnew,kx);
        delete [] xnew;
        X::function<type> F(xnew,kx);
        for(int i=0;i<ky;i++)
        {
            cn::complex<type> buf = f.getf((kx-1)/2+1,i)*cn::con(f.getf((kx-1)/2+1,i));
            F(i) = buf.Real();
        }
        type N = F.integralS2(0,ky);

        for(int i=0;i<kx;i++)
        {
            for(int j=0;j<ky;j++)
            {
                cn::complex<type> buf = f.getf(i,j)*cn::con(f.getf(i,j));
                F(j) = buf.Real();
            }
            res(i) = F.integralS2(0,ky-1)/N;
        }

    #ifdef psi_spline
        res = X::spline(res,psi_function_spline);
    #endif

        return res;
    }

    template<typename type>
    X::function<type> function<type>::conditionalPsi() const
    {
        int kx=f.getnx();
        int ky=f.getny();
        int n0 = (ky-1)/2;
        type *xnew = new type[kx];
        for(int i=0;i<kx;i++)
        {
            xnew[i] = f.getx(i,n0).Real();
        }
        X::function<type> res(xnew,kx);
        delete [] xnew;

        cn::complex<type> buf = f.getf(n0,n0)*cn::con(f.getf(n0,n0));
        type N = buf.Real();

        for(int i=0;i<kx;i++)
        {
            cn::complex<type> buf = f.getf(i,n0)*cn::con(f.getf(i,n0));
            res(i) = buf.Real()/N;
        }

    #ifdef psi_spline
        res = X::spline(res,psi_function_spline);
    #endif

        return res;
    }

}

#endif // PSIFUNCTION_H
