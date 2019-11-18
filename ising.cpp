#include "ising.h"
ofstream myfile;

//denne funksjonen henter naboene til punktet (x ,i+add)
inline int periodeBoundry(int i, int lim, int add) {
    return (i+lim+add) % (lim);
}


//funksjonen setter opp systemets start matrise, energi og magnetisering
void inisialiser(int spins, mat &matrise, double &E, double &M, int i) {
    if(i == 1) { //hvis jeg vil ha matrisen ordnet
        for(int x=0; x<spins; x++) {
            for(int y=0; y<spins; y++) {
                matrise(x, y) = 1.0;
            }
        }
    }
    else if(i == 0) { //hvis jeg vil ha matrisen tilfeldig
        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_real_distribution<double> RandomNumberGenerator(0.0, 1.0);

        for(int x=0; x<spins; x++) {
            for(int y=0; y<spins; y++) {
                if(RandomNumberGenerator(gen) < 0.5) {
                    matrise(x, y) = -1.0;
                }
                else{
                    matrise(x, y) = 1.0;
                }
            }
        }

        //finner start energien og magnetisering
        for(int x=0; x<spins; x++) {
            for(int y=0; y<spins; y++) {
                M += static_cast<double>(matrise(x, y));
                E -= static_cast<double>(matrise(x, y) * (matrise(periodeBoundry(x, spins, -1), y) + matrise(x, periodeBoundry(y, spins, -1))));
            }
        }

    }
}

//selve Monte Carlo algoritmen, med Metropolis algoritmen inni
void mcMetropolis(int mcSyklus, int spins, vec &expect, mat &matrise, vec &omega, double &E, double &M) {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> RandomNumberGenerator(0.0, 1.0);

    double n = static_cast<double>(spins);

    for(int mc=1; mc<mcSyklus; mc++) {

        for(int x=0; x<spins; x++) {
            for(int y=0; y<spins; y++) {
                //finner et tilfeldig punkt i matrisen
                int i = static_cast<int>(RandomNumberGenerator(gen)*n);
                int j = static_cast<int>(RandomNumberGenerator(gen)*n);

                //henter naboene i punktet og finner energien
                int nabo = matrise(i, periodeBoundry(j, spins, -1)) + matrise(periodeBoundry(i, spins, -1), j) + matrise(i, periodeBoundry(j, spins, 1)) + matrise(periodeBoundry(i, spins, 1), j);
                int dE = 2*matrise(i, j)*nabo;

                //metropolis algoritmen
                if(RandomNumberGenerator(gen) <= omega(dE+8)) {
                    matrise(i, j) *= -1.0;
                    M += static_cast<double>(2*matrise(i, j));
                    E += static_cast<double>(dE);
                }
            }
        }

        expect(0) += E; expect(1) += E*E;
        expect(2) += M; expect(3) += M*M;
        expect(4) += fabs(M);
    }
}

//funksjonen som sender resultatene til et dokument
void resultat(int spins, int mcSyklus, double t, vec expect) {
    double n = 1.0/static_cast<double>(mcSyklus);

    double E_ = expect(0)*n; double E2_ = expect(1)*n;
    double M_ = expect(2)*n; double M2_ = expect(3)*n;
    double Mabs_ = expect(4)*n;

    double E_var = (E2_ - E_*E_)/spins/spins;
    double M_var = (M2_ - Mabs_*Mabs_)/spins/spins;

    myfile << setiosflags(ios::showpoint | ios::uppercase);

    myfile << setw(15) << setprecision(8) << t;
    myfile << setw(15) << setprecision(8) << E_/spins/spins;
    myfile << setw(15) << setprecision(8) << E_var;
    myfile << setw(15) << setprecision(8) << E_var/(t*t);
    myfile << setw(15) << setprecision(8) << M_/spins/spins;
    myfile << setw(15) << setprecision(8) << M_var/(t*t);
    myfile << setw(15) << setprecision(8) << Mabs_/spins/spins << endl;

}

void isingModel(int spins, int mcSyklus, double temp0, double tempN, double tempSteg) {
    //opne fil
    string k = "IsingModel_L" + to_string(spins) + ".txt";
    myfile.open(k);

    myfile << setw(15) << setprecision(8) << "t";
    myfile << setw(15) << setprecision(8) << "E";
    myfile << setw(15) << setprecision(8) << "E_var";
    myfile << setw(15) << setprecision(8) << "C_V";
    myfile << setw(15) << setprecision(8) << "M";
    myfile << setw(15) << setprecision(8) << "X";
    myfile << setw(15) << setprecision(8) << "M_abs" << endl;

    clock_t start, stopp;
    start = clock();

    for(double t=temp0; t<=tempN; t+=tempSteg) {
        vec expect = zeros<mat>(5);
        double E = 0.0; double M = 0.0;
        mat matrise = zeros<mat>(spins, spins);

        inisialiser(spins, matrise, E, M, 0); //orderd form 0, tilfeldig form 1 (for den siste parameteret)

        vec omega = zeros<mat>(17);
        for(int i=-8; i<=8; i+=4) {
            omega(i+8) = exp(-i/t);
        }

        mcMetropolis(mcSyklus, spins, expect, matrise, omega, E, M);

        resultat(spins, mcSyklus, t, expect);
    }
    //steng fil
    stopp = clock();
    myfile << "\n";
    myfile << "Tid: " << ( (stopp - start)/CLOCKS_PER_SEC) << endl;
    myfile.close();
}

//ising modelen parallellisert
void parallel(int spins, int mcSyklus, double temp0, double tempN, double tempSteg) {
    //opne fil
    string k = "IsingParallell_L" + to_string(spins) + ".txt";
    myfile.open(k);

    myfile << setw(15) << setprecision(8) << "t";
    myfile << setw(15) << setprecision(8) << "E";
    myfile << setw(15) << setprecision(8) << "E_var";
    myfile << setw(15) << setprecision(8) << "C_V";
    myfile << setw(15) << setprecision(8) << "M";
    myfile << setw(15) << setprecision(8) << "X";
    myfile << setw(15) << setprecision(8) << "M_abs" << endl;

    clock_t start, stopp;
    start = clock();

    for(double t=temp0; t<=tempN; t+=tempSteg) {
        vec expect = zeros<mat>(5);
        double E = 0.0; double M = 0.0;
        mat matrise = zeros<mat>(spins, spins);

        inisialiser(spins, matrise, E, M, 0); //1 hvis ordnet, 0 hvis tilfeldig

        vec omega = zeros<mat>(17);
        for(int i=-8; i<=8; i+=4) {
            omega(i+8) = exp(-i/t);
        }

        vec sannsynlighet = zeros<mat>(17);

        int mc;
        double tE, tM, tE2, tM2, tMabs;
#pragma omp parallel default(shared) private(mc, E, M) reduction(+:tE, tM, tE2, tM2, tMabs)
        E = M = 0.0;
        tE = tM = tE2 = tM2 = tMabs = 0.0;

        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_real_distribution<double> RandomNumberGenerator(0.0, 1.0);

        double n = static_cast<double>(spins);
        int p = 0;

#pragma omp for
        for(mc=1; mc<mcSyklus; mc++) {

            for(int x=0; x<spins; x++) {
                for(int y=0; y<spins; y++) {
                    int i = static_cast<int>(RandomNumberGenerator(gen)*n);
                    int j = static_cast<int>(RandomNumberGenerator(gen)*n);

                    int nabo = matrise(i, periodeBoundry(j, spins, -1)) + matrise(periodeBoundry(i, spins, -1), j) + matrise(i, periodeBoundry(j, spins, 1)) + matrise(periodeBoundry(i, spins, 1), j);
                    int dE = 2*matrise(i, j)*nabo;

                    sannsynlighet(dE+8) += 1;

                    if(RandomNumberGenerator(gen) <= omega(dE+8)) {
                        matrise(i, j) *= -1.0;
                        M += static_cast<double>(2*matrise(i, j));
                        E += static_cast<double>(dE);
                        p += 1;
                    }
                 }
             }

             tE += E; tE2 += E*E; tM += M; tM2 += M*M; tMabs += fabs(M);

        }
        expect(0) = tE; expect(1) = tE2;
        expect(2) = tM; expect(3) = tM2;
        expect(4) = tMabs;
        resultat(spins, mcSyklus, t, expect);

    }
    stopp = clock();
    myfile << "\n";
    myfile << "Tid: " << ( (stopp - start)/CLOCKS_PER_SEC) << endl;
    myfile.close();
}
