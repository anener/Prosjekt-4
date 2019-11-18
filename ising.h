#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>
#include "time.h"


using namespace std;
using namespace arma;

inline int periodeBoundry(int, int, int);
void inisialiser(int, mat &, double &, double &, int);
void mcMetropolis(int, int, vec &, mat &, vec &, double &, double &);
void resultat(int, int, double, vec);

void isingModel(int, int, double, double, double);
void parallel(int, int, double, double, double);




