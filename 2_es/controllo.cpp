#include "iostream"
#include "fstream"
#include "iomanip"
#include "cmath"
#include "integrali.h"
#include <vector>
#include <math.h>

using namespace std;

//funzione da integrare
double f (double x){
    double ris = pow(x, 7)*exp(-x);
    return ris;
}

//funzione da integrare riscalata per il metodo di Legendre
double flegendre (double x){
    double ris = (5./2.)*pow(((5./2.)*(x+1.)), 7)*exp(-(5./2.)*(x+1.));
    return ris;
}

//funzioni da integrare riscritte per poter utilizzare il metodo di Laguerre
double flaguerre1 (double x){
    double ris = pow(x, 7);
    return ris;
}

double flaguerre2(double x){
    double ris = pow(x+5, 7)*exp(-5);
    return ris;
}


int main (int argc, char** argv){
    
    long double ris_true = 672.19323731283680927664797017826546492574500186733942065874221332363336448; //risultato vero dell'integrale calcolato con Wolfram Alpha
    int N_trapezio = 2501;
    int N_simpson = 2501;
    int N_romberg = pow(2, 7)+1;

    int N = 8;

    double xmin = 0;
    double xmax = 5;
    //bool V = true;

    vector<double> ris_legendre;
    vector<double> ris_laguerre;

// calcolo dell'integrale con il metodo del trapezio e stampa su file della discrepanza con il valore vero
    fstream file1;
    file1.open("trap.dat", ios_base::out);
    for (int i=2; i<=N_trapezio; i++){
        double ris = trapezio(xmin, xmax, i, f);
        double err = abs(ris - ris_true);
        file1 << 1/static_cast<double>(i) << "\t" << setprecision(100) << err << endl;
    }
    file1.close();

// calcolo dell'integrale con il metodo di Simpson e stampa su file della discrepanza con il valore vero
    fstream file2;
    file2.open("simp.dat", ios_base::out);
    for (int i=3; i<=N_simpson; i+=2){
        double ris = simpson(xmin, xmax, i, f);
        double err = abs(ris - ris_true);
        file2 << 1/static_cast<double>(i) << "\t" << setprecision(100) << err << endl;
    }
    file2.close();

 
// calcolo dell'integrale con il metodo di Romberg ed i mwtodi di Gauss e stampa su file della discrepanza con il valore vero

    double ris_romb = romberg (xmin, xmax, N_romberg, f);

    for (int i=2; i<=N; i= i*2){
        ris_legendre.push_back(legendre(i, flegendre));
        ris_laguerre.push_back(laguerre(i, flaguerre1) - laguerre (i, flaguerre2));
    }

    fstream file3;
    file3.open("risultati.dat", ios_base::out);

    file3 << "trapezio" << "\t" << setprecision(100) << trapezio(xmin, xmax, N_trapezio, f) << "\t" << setprecision(100) << abs(trapezio(xmin, xmax, N_trapezio, f) - ris_true) << endl;

    file3 << "simpson" << "\t" << setprecision(100) << simpson(xmin, xmax, N_simpson, f) << "\t" << setprecision(100) << abs(simpson(xmin, xmax, N_simpson, f) - ris_true) << endl;

    file3 << "romberg" << "\t" << setprecision(100) << ris_romb << "\t" << setprecision(100) << abs(ris_romb - ris_true) << endl;
    for (int i=0; i<ris_legendre.size(); i++){
        file3 << "legendre " << i+1 << "\t" << setprecision(100) << ris_legendre[i] << "\t" << setprecision(100) << abs(ris_legendre[i] - ris_true) << endl;
    }
    for (int i=0; i<ris_laguerre.size(); i++){
        file3 << "laguerre " << i+1 << "\t" << setprecision(100) << ris_laguerre[i] << "\t" << setprecision(100) << abs(ris_laguerre[i] - ris_true) << endl;
    }

    file3.close();

    cout << simpson(xmin, xmax, N_simpson, f) << endl;


    return 0;
}