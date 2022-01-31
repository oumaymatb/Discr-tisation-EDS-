#include <iostream>
#include "../TP_EDS/functions.hpp"
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <iomanip>
#include <functional>
#include <numeric>
#include <chrono> 

#define _USE_MATH_DEFINES
using namespace std;

int main(){

    double T = 1.;
    double sig = 0.2;
    double r = 0.05;
    double S0 = 100.0;
    int N =20;
    int M = 100000;
    int P = 6;

    vector<double> erreul(P,0.0);
    vector<double> liceul(P,0.0);
    vector<double> errmil(P,0.0);
    vector<double> licmil(P,0.0);
    vector<double> Npas(P,0.0);
    vector<double> S(M,S0);
    vector<double> Se(M,S0);
    vector<double> Sm(M,S0);
    vector<double> maxerreul(M);
    vector<double> maxerrmil(M);

    for(int i=0;i<P;i++){
        double dt = T/(double)(N);
        double B = (r-sig*sig*0.5)*dt;
        double C = r*dt;

        vector<double> init_S(M,S0);
        vector<double> init_err(M);
        S = init_S;
        Se = init_S;
        Sm = init_S;
        maxerreul = init_err;
        maxerrmil = init_err;


        for(int k=0;k<N;k++){
            vector<double> g = gaussian(0.0,1.0,M);
            vector<double> A = g;
            for_each(A.begin(),A.end(),
                [&](double& x){ x = sig*sqrt(dt)*x;  }  );

            //euler
            transform(A.begin(),A.end(),Se.begin(),Se.begin(),
                [&](double x,double y){ return y*(1+x+r*dt);  } );
            
            //milstein
            transform(A.begin(),A.end(),Sm.begin(),Sm.begin(),
                [&](double x,double y){ return y*(1 + x + 0.5*x*x + B);  } );

            //EDS
            transform(A.begin(),A.end(),S.begin(),S.begin(),
                [&](double x,double y){ return y*exp(B+x) ; } );

            vector<double> absolute_err(M);
            transform(S.begin(),S.end(),Se.begin(),absolute_err.begin(),
                [](double x, double y){ return abs(x-y); }  );
            
            transform(absolute_err.begin(),absolute_err.end(),maxerreul.begin(),maxerreul.begin(),
                [](double x, double y){ return max(x,y); }  );
            
            vector<double> absolute_mil(M);
            transform(S.begin(),S.end(),Sm.begin(),absolute_mil.begin(),
                [](double x, double y){ return abs(x-y); }  );

            transform(absolute_mil.begin(),absolute_mil.end(),maxerrmil.begin(),maxerrmil.begin(),
                [](double x, double y){  return max(x,y); } );
        }

        double someul = 0.0;
        for_each(maxerreul.begin(),maxerreul.end(),
            [&](double x){ someul += x*x; } );
        double careul = 0.0;
        for_each(maxerreul.begin(),maxerreul.end(),
            [&](double x){ careul += pow(x,4.0); } );
        erreul[i] = someul/(double)(M);
        double liceul_i = careul/(double)(M)-(someul/(double)(M))*(someul/(double)(M));
        liceul[i] = 1.96*sqrt(liceul_i/(double)(M));

        double sommil = 0.0;
        for_each(maxerrmil.begin(),maxerrmil.end(),
            [&](double x){ sommil += x*x; } );
        double carmil = 0.0;
        for_each(maxerrmil.begin(),maxerrmil.end(),
            [&](double x){ carmil += pow(x,4.0); } );
        errmil[i] = sommil/(double)(M);
        double licmil_i = carmil/(double)(M)-(sommil/(double)(M))*(sommil/(double)(M));
        licmil[i] = 1.96*sqrt(licmil_i/(double)(M));

        Npas[i] = N;

        N = 2*N;

    }

cout << "Erreur euler" << endl;
print(erreul);
cout << "LIC Euler" << endl;
print(liceul);
cout << "Erreur Milstein" << endl;
print(errmil);
cout << "LIC Milstein" << endl;
print(licmil);
cout << "Npas" << endl;
print(Npas);


}