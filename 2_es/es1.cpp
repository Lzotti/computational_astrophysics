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
    double ris = 1.0 / (1 + pow(x, 2));
    return ris;
}

int main (int argc, char** argv){
    
    long double ris_true = M_PI_4; 
    int N_trapezio = 321;
    int N_simpson = 321;
    int N_romberg = pow(2, 7)+1;

    int N = 8;

    double xmin = 0;
    double xmax = 1;
    //bool V = true;

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

// calcolo dell'integrale con il metodo del rettangolo_left e stampa su file della discrepanza con il valore vero
    fstream file4;
    file4.open("rect_left.dat", ios_base::out);
    for (int i=2; i<=N_trapezio; i++){
        double ris = rettangolo_left(xmin, xmax, i, f);
        double err = abs(ris - ris_true);
        file4 << 1/static_cast<double>(i) << "\t" << setprecision(100) << err << endl;
    }
    file4.close();

// calcolo dell'integrale con il metodo del rettangolo_right e stampa su file della discrepanza con il valore vero
    fstream file5;
    file5.open("rect_right.dat", ios_base::out);
    for (int i=2; i<=N_trapezio; i++){
        double ris = rettangolo_right(xmin, xmax, i, f);
        double err = abs(ris - ris_true);
        file5 << 1/static_cast<double>(i) << "\t" << setprecision(100) << err << endl;
    }
    file5.close();

// calcolo dell'integrale con il metodo del rettangolo_midpoint e stampa su file della discrepanza con il valore vero
    fstream file6;
    file6.open("rect_midpoint.dat", ios_base::out);
    for (int i=2; i<=N_trapezio; i++){
        double ris = rettangolo_midpoint(xmin, xmax, i, f);
        double err = abs(ris - ris_true);
        file6 << 1/static_cast<double>(i) << "\t" << setprecision(100) << err << endl;
    }
    file6.close();

// calcolo dell'integrale con il metodo di Romberg ed i mwtodi di Gauss e stampa su file della discrepanza con il valore vero

    double ris_romb = romberg (xmin, xmax, N_romberg, f);

    fstream file3;
    file3.open("risultati.dat", ios_base::out);

    file3 << "trapezio" << "\t" << setprecision(100) << trapezio(xmin, xmax, N_trapezio, f) << "\t" << setprecision(100) << abs(trapezio(xmin, xmax, N_trapezio, f) - ris_true) << endl;

    file3 << "simpson" << "\t" << setprecision(100) << simpson(xmin, xmax, N_simpson, f) << "\t" << setprecision(100) << abs(simpson(xmin, xmax, N_simpson, f) - ris_true) << endl;

    file3 << "romberg" << "\t" << setprecision(100) << ris_romb << "\t" << setprecision(100) << abs(ris_romb - ris_true) << endl;

    file3.close();

    cout << simpson(xmin, xmax, N_simpson, f) << endl;


    return 0;
}