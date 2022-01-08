#ifndef _SAVEDATA_H_
#define _SAVEDATA_H_

#include <iostream>
#include <vector>
#include "kisskaMath/psifunction.h"
using namespace std;

struct fileNames
{
    fileNames(
            const char * Re,
            const char * Im,
            const char * x,
            const char * y):
        Re(Re),
        Im(Im),
        x(x),
        y(y)
    {}

    void setLabels(
            const char * lRe,
            const char * lIm,
            const char * lx,
            const char * ly)
    {
        labelRe = lRe;
        labelIm = lIm;
        labelx = lx;
        labely = ly;
    }

    void setLegend(
            const char * lRe,
            const char * lIm)
    {
        legendRe = lRe;
        legendIm = lIm;
    }

    const char * Re; const char * Im; const char * x; const char * y;
    const char * labelRe; const char * labelIm; const char * labelx; const char * labely;
    const char * legendRe; const char * legendIm;
};

template<typename type>
struct outInfo
{
    outInfo(X::function<type> psi, fileNames Names): psi(psi), Names(Names)
    {}
    X::function<type> psi;
    fileNames Names;
};

namespace  psi{
    void saveData(psi::function<double> psi, const char *name, const char *MODE, fileNames Names)
    {
        psi.setout(7,7);
        psi.Re().savef(Names.Re);
        psi.Im().savef(Names.Im);
        psi.Re().savex(Names.x);
        psi.Re().savey(Names.y);

        FILE *fp;
        fp=freopen(name,MODE,stdout);

        cout<<Names.labelx<<" =  load('";
        cout<<Names.x;
        cout<<"');"<<endl;

        cout<<Names.labely<<" =  load('";
        cout<<Names.y;
        cout<<"');"<<endl;

        cout<<Names.labelRe<<" =  load('";
        cout<<Names.Re;
        cout<<"');"<<endl;

        cout<<Names.labelIm<<" = load('";
        cout<<Names.Im;
        cout<<"');"<<endl;

        cout<<"figure(1);"<<endl;
        cout<<"surf("<<Names.labelx<<","<<Names.labely<<","<<Names.labelRe<<");"<<endl;
        cout<<"figure(2);"<<endl;
        cout<<"surf("<<Names.labelx<<","<<Names.labely<<","<<Names.labelIm<<");"<<endl;
        cout<<endl;

        fclose(fp);
        freopen("CON","w",stdout);
    }
}

namespace X
{
    template<typename type>
    void saveData(X::function<type> psi, const char *name, const char *MODE, fileNames Names)
    {
        psi.setout(7,7);
        psi.savef(Names.Re);
        psi.savex(Names.x);
        FILE *fp;
        fp=freopen(name,MODE,stdout);

        cout<<Names.labelx<<" =  load('";
        cout<<Names.x;
        cout<<"');"<<endl;

        cout<<Names.labelRe<<" =  load('";
        cout<<Names.Re;
        cout<<"');"<<endl;

        cout<<"plot("<<Names.labelx<<","<<Names.labelRe<<");"<<endl;

        fclose(fp);
        freopen("CON",MODE,stdout);
    }

    template<typename type>
    void saveData(std::vector<outInfo<type>> Info, const char *name, const char *MODE)
    {
        for(unsigned int i=0;i<Info.size();i++)
        {
            Info.at(i).psi.setout(7,7);
            Info.at(i).psi.savef(Info.at(i).Names.Re);
            Info.at(i).psi.savex(Info.at(i).Names.x);
        }
        FILE *fp;
        fp=freopen(name,MODE,stdout);
        for(unsigned int i=0;i<Info.size();i++)
        {
            cout<<Info.at(i).Names.labelx<<" =  load('";
            cout<<Info.at(i).Names.x;
            cout<<"');"<<endl;

            cout<<Info.at(i).Names.labelRe<<" =  load('";
            cout<<Info.at(i).Names.Re;
            cout<<"');"<<endl;
        }
        for(unsigned int i=0;i<Info.size();i++)
        {
            cout<<"plot("<<Info.at(i).Names.labelx<<","<<Info.at(i).Names.labelRe<<"); hold on;"<<endl;
        }
        cout<<"grid on"<<endl;
        cout<<"legend('"<<Info.at(0).Names.legendRe<<"'";
        for(unsigned int i=1;i<Info.size();i++)
        {
            cout<<",'"<<Info.at(i).Names.legendRe<<"'";
        }
        cout<<");"<<endl;


        fclose(fp);
        freopen("CON",MODE,stdout);
    }
}


namespace XY
{
    void saveData(XY::function<double> psi, const char *name, const char *MODE, fileNames Names)
    {
        psi.setout(7,7);
        psi.savef(Names.Re);
        psi.savex(Names.x);
        psi.savey(Names.y);

        FILE *fp;
        fp=freopen(name,MODE,stdout);

        cout<<Names.labelx<<" =  load('";
        cout<<Names.x;
        cout<<"');"<<endl;

        cout<<Names.labely<<" =  load('";
        cout<<Names.y;
        cout<<"');"<<endl;

        cout<<Names.labelRe<<" =  load('";
        cout<<Names.Re;
        cout<<"');"<<endl;

        cout<<"figure(1);"<<endl;
        cout<<"surf("<<Names.labelx<<","<<Names.labely<<","<<Names.labelRe<<");"<<endl;
        cout<<endl;

        fclose(fp);
        freopen("CON","w",stdout);
    }
}

#endif // SAVEDATA_H
