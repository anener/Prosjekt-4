#include "ising.h"
using namespace std;

//isingModel(int spins, int mc_syklus, double temp0, double tempN, double dTemp)

int main()
{

    int spins; //antall spinns?
    int mc_syklus; //antall Monte Carlo syklus

    double temp0; //start temperaturen
    double tempN; //slutt temperaturen
    double dTemp; //temperatur steg

    cout << "Antall spins: ";
    cin >> spins;

    cout << "Temperatur start: ";
    cin >> temp0;

    cout << "Temperatur stopp: ";
    cin >> tempN;

    cout << "Temperatur steg: ";
    cin >> dTemp;

    cout << "Monte Carlo syklus: ";
    cin >> mc_syklus;

    //oppgave 4b: spins=2, temp0=1.0, tempN=1.0, dTemp=1.0
    isingModel(spins, mc_syklus, temp0, tempN, dTemp);

    parallel(spins, mc_syklus, temp0, tempN, dTemp);


    return 0;
}

