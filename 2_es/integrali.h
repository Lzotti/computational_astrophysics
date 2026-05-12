#ifndef integrali_h
#define integrali_h

#include "cmath"
#include "string.h"
#include "iostream"
#include <math.h>
#include "fstream"
#include <vector>

using namespace std;

//Formula del Rettangolo left point
double rettangolo_left (double xmin, double xmax, int N, double f(double)){
    double h = (xmax-xmin)/(static_cast<double>(N)); //calcolo del passo
    double sum = 0;
//applicazione della formula
    for(int i=0; i< N; i++){
        sum += f(xmin+ static_cast<double>(i)*h);
    }
    double ris = h*sum;
    return ris;
}

double rettangolo_right (double xmin, double xmax, int N, double f(double)){
    double h = (xmax-xmin)/(static_cast<double>(N)); //calcolo del passo
    double sum = 0;
//applicazione della formula
    for(int i=1; i<= N; i++){
        sum += f(xmin+ static_cast<double>(i)*h);
    }
    double ris = h*sum;
    return ris;
}

double rettangolo_midpoint (double xmin, double xmax, int N, double f(double)){
    double h = (xmax-xmin)/(static_cast<double>(N)); //calcolo del passo
    double sum = 0;
//applicazione della formula
    for(int i=0; i< N; i++){
        sum += f(xmin+ (static_cast<double>(i)+0.5)*h);
    }
    double ris = h*sum;
    return ris;
}

//----------------------------------------------------------------------------------------------------------------------

// Formula del trapezio
double trapezio ( double xmin,  double xmax, int N, double f(double)){

    double h = (xmax-xmin)/(static_cast<double>(N-1)); //calcolo del passo
    double sum = 0;
//applicazione della formula
    for(int i=1; i< N-1; i++){
        sum += f(xmin+ static_cast<double>(i)*h);
    }
    sum += (f(xmin)+f(xmax))/2;
    double ris = h*sum;
    return ris;
}
//----------------------------------------------------------------------------------------------------------------------

//Formula di Simpson
double simpson (double xmin, double xmax, int N, double f(double)){
//controllo che N sia dispari
    if( N%2==0){
        cerr << "Errore: N deve essere dispari per il metodo di Simpson" << std::endl;
        abort();    
    }
    double h = (xmax-xmin)/(static_cast<double>(N)-1); //calcolo del passo
    double partsum = 0;
    double sum = 0;
//applicazione della formula
    for(int i=1; i<N -1; i++){
        partsum = f(xmin+ static_cast<double>(i)*h);
        if( i % 2 == 0){
            partsum = partsum*2;
        }
        else{
            partsum = partsum*4;
        } 
        sum += partsum;
    }
    
    sum += f(xmin) + f(xmax);

    double ris = h*sum/3.;
    return ris;
}

//----------------------------------------------------------------------------------------------------------------------
//Formula di Romberg
double romberg (double xmin, double xmax, int N, double f(double), bool V=false){ 
//controllo che N sia del tipo 2^J + 1 per poter utilizzare il metodo
    if(((N-1) & (N-2))){
        cerr << "Errore: N deve essere del tipo 2^J + 1" << endl;
        abort();
    }

    double h, sum = 0;
// calcolo del numero di livelli J
    int J  = 0;
    int tempN = N-1;
    while (tempN != 1){
        tempN = tempN/2;
        J++;       
    }
//-----------------------------------------------------------------------
// creo una matrice e la inizializzo a 0
    vector<vector<double>> Mat(J + 1, vector<double>(J + 1, 0.0));

//-----------------------------------------------------------------------

// Calcolo del primo valore usando il metodo del trapezio
    Mat[0][0]= trapezio(xmin, xmax, 2, f); 

// Calcolo di tutti gli altri valori
    for (int n = 1; n <=J ; n++){
        for (int m = 0; m <= n; m++){
            // se m=0 uso la formula ricorsiva del trapezio
            if (m == 0){
                h = (xmax-xmin)/(pow(2,n));
                sum = 0;
                for (int k = 1; k <= (pow(2, n))/2; k++){                    
                    sum += f(xmin + (2 * k - 1) * h);
                }
                Mat[n][0] = 0.5 * Mat[n - 1][0] + h * sum;     
            }
            else{
                Mat[n][m] = (pow(4, m) * Mat[n][m - 1] - Mat[n - 1][m - 1])/(pow(4,m) - 1);
            }
        }
    }
//----------------------------------------------------------------------
// Restituisco il valore dell'integrale in base alla precisione richiesta
    if (V) {
        int n, m;
        cout << "inserisci riga e poi colonna: " << endl;
        cin >> n >> m;
        return Mat[n][m];
    }
    else{
        return Mat[J][J];
    }
}

//----------------------------------------------------------------------------------------------------------------------
//METODI DI GAUSS

//polinomi di legendre
double legendre (int N, double f(double)){

    vector<double> x(N, 0.);
    vector<double> p(N, 0.);

    // Apertura del file contenente i nodi e i pesi
    fstream Pesi_Gaussiani;
    string nomefile = "Legendre" + to_string(N) + ".dat";
    Pesi_Gaussiani.open(nomefile, ios_base::in);

    // Controllo se il file è stato aperto correttamente
    if (!Pesi_Gaussiani.is_open()) {
        cerr << "Errore: impossibile aprire il file " << nomefile << std::endl;
        abort();
    }

    // Lettura dei nodi e dei pesi dal file
    for(int i=0; i<N; i++){
        Pesi_Gaussiani >> x[i] >> p[i];
    }
    Pesi_Gaussiani.close();

   // Calcolo dell'integrale 
    double sum = 0.;
    for (int i=0; i<N; i++){
        sum +=  p[i] * f(x[i]);
    }

    return sum;
}

//----------------------------------------------------------------------------------------------------------------------

//polinomi di laguerre
double laguerre(int N, double f(double)){  

    vector<double> x(N, 0.);
    vector<double> p(N, 0.);

    // Apertura del file contenente i nodi e i pesi
    fstream Pesi_Gaussiani;
    string nomefile = "Laguerre" + to_string(N) + ".dat";
    Pesi_Gaussiani.open(nomefile, ios_base::in);

    // Controllo se il file è stato aperto correttamente
    if (!Pesi_Gaussiani.is_open()) {
        cerr << "Errore: impossibile aprire il file " << nomefile << std::endl;
        abort();
    }

    // Lettura dei nodi e dei pesi dal file
    for(int i=0; i<N; i++){
        Pesi_Gaussiani >> x[i] >> p[i];
    }
    Pesi_Gaussiani.close();

   // Calcolo dell'integrale 
    double sum = 0.;
    for (int i=0; i<N; i++){
        sum +=  p[i] * f(x[i]);
    }

    return sum;
}

//----------------------------------------------------------------------------------------------------------------------

//polinomi di hermite
double hermite(int N, double f(double)){  

    vector<double> x(N, 0.);
    vector<double> p(N, 0.);

    // Apertura del file contenente i nodi e i pesi
    fstream Pesi_Gaussiani;
    string nomefile = "Hermite" + to_string(N) + ".dat";
    Pesi_Gaussiani.open(nomefile, ios_base::in);

    // Controllo se il file è stato aperto correttamente
    if (!Pesi_Gaussiani.is_open()) {
        cerr << "Errore: impossibile aprire il file " << nomefile << std::endl;
        abort();
    }

    // Lettura dei nodi e dei pesi dal file
    for(int i=0; i<N; i++){
        Pesi_Gaussiani >> x[i] >> p[i];
    }
    Pesi_Gaussiani.close();

   // Calcolo dell'integrale 
    double sum = 0.;
    for (int i=0; i<N; i++){
        sum +=  p[i] * f(x[i]);
    }

    return sum;
}

//----------------------------------------------------------------------------------------------------------------------
//Formule aperte

//Formula del midpoint a 5 intervalli 
double midpoint(double xmin, double xmax, double f(double)){

    double h = (xmax-xmin)/5.;
    double partsum=0;
    double sum= 0;
    double ris = 0;

    for (int i=1; i<5; i++){
        partsum= f(xmin + h*i);
        if (i==1 || i==4) partsum = partsum*11;
        sum += partsum;
    }
    
    ris = (5./24.)*h*sum;
    return ris;
}
#endif
