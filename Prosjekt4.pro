TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        ising.cpp \
        main.cpp

INCLUDEPATH += C:\armadillo-9.700.2\include
DEPENDPATH += C:\armadillo-9.700.2\include

QMAKE_CXXFLAGS += -fopenmp

LIBS += \
    -LC:\armadillo-9.700.2\examples\lib_win64 \
    -llapack_win64_MT \
    -lblas_win64_MT \
    -fopenmp


HEADERS += \
    ising.h



