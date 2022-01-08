TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp

HEADERS += \
    kisskaMath/complex.h \
    kisskaMath/matrix.h\
    kisskaMath/function.h\
    kisskaMath/psifunction.h\
    timer.h\
    saveData.h
