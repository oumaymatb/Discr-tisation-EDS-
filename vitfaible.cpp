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
    double K = 110.;
    double S0 = 100.0;
    int N =1;
    int M = 100000;
    int P = 6;

    double d1 = (log(S0/K)+r*T)/(sig*sqrt(T)) + sig*sqrt(T)/(double)(2);
    double d2 = d1 - sig*sqrt(T);
    double BS = K*exp(-r*T)*NN(-d2) - S0*NN(-d1);

    vector<double> erreul(P,0.0);
    vector<double> liceul(P,0.0);
    vector<double> conterreul(P,0.0);
    vector<double> contliceul(P,0.0);
    vector<double> errmil(P,0.0);
    vector<double> licmil(P,0.0);
    vector<double> conterrmil(P,0.0);
    vector<double> contlicmil(P,0.0);
    vector<double> Npas(P,0.0);
    vector<double> S(M,S0);
    vector<double> Se(M,S0);
    vector<double> Sm(M,S0);


    for(int i=0;i<P;i++){
        double dt = T/(double)(N);
        double sigdt = sig*sqrt(dt);
        double rdt = r*dt;
        double der = rdt - sigdt*sigdt*0.5;

        vector<double> init_S(M,S0);
        S = init_S;
        Se = init_S;
        Sm = init_S;


        for(int k=0;k<N;k++){
            vector<double> g = gaussian(0.0,1.0,M);
            vector<double> A = g;
            for_each(A.begin(),A.end(),
                [&](double& x){ x = sigdt*x;  }  );

            //euler
            transform(A.begin(),A.end(),Se.begin(),Se.begin(),
                [&](double x,double y){ return y*(1+x+rdt);  } );
            
            //milstein
            transform(A.begin(),A.end(),Sm.begin(),Sm.begin(),
                [&](double x,double y){ return y*(1 + x + 0.5*x*x + der);  } );

            //EDS
            transform(A.begin(),A.end(),S.begin(),S.begin(),
                [&](double x,double y){ return y*exp(der+x) ; } );

        }
        // EDS put
        vector<double> paymc;
        for_each(S.begin(),S.end(),       
            [&](double x){  paymc.push_back( (K>x) ? (K-x):0.0) ; } );

        // Euler Put
        vector<double> payeul;
        for_each(Se.begin(),Se.end(), 
            [&](double x){ payeul.push_back( (K>x) ? (K-x):0.0 ); } );

        double puteul = accumulate(payeul.begin(),payeul.end(),0.0);
        double careul = 0.0;
        for_each(payeul.begin(),payeul.end(),
            [&](double x){ careul += x*x; } );
        
        double contputeul = 0.0;

        vector<double> diff_eul(M);
        transform(paymc.begin(),paymc.end(),payeul.begin(),diff_eul.begin(),
            [](double x, double y){ return x-y;  }  );
        
        contputeul = accumulate(diff_eul.begin(),diff_eul.end(),0.0);
        double contcareul = 0.0;
        for_each(diff_eul.begin(),diff_eul.end(),
            [&](double x){ contcareul += x*x; } );
        
        // Milstein Put
        vector<double> paymil;
        for_each(Sm.begin(),Sm.end(),
            [&](double x){ paymil.push_back( (K>x) ? (K-x):0.0 ); } );
        double putmil = accumulate(paymil.begin(),paymil.end(),0.0);
        double carmil = 0.0;
        for_each(paymil.begin(),paymil.end(),
            [&](double x){ carmil += x*x; } );

        vector<double> diff_mil(M);
        transform(paymc.begin(),paymc.end(),paymil.begin(),diff_mil.begin(),
            [](double x, double y){ return x-y;  }  );
        
        double contputmil = accumulate(diff_mil.begin(),diff_mil.end(),0.0);
        double contcarmil = 0.0;
        for_each(diff_mil.begin(),diff_mil.end(),
            [&](double x){ contcarmil += x*x; } );
        
        erreul[i] = exp(-r*T)*puteul/(double)(M) - BS;
        double liceul_i =  careul/(double)(M) - pow(puteul/(double)(M),2.0);
        liceul_i = liceul_i/(double)(M);
        liceul[i] = 1.96 * exp(-r*T)*sqrt( liceul_i);
        conterreul[i] = exp(-r*T)*contputeul/(double)(M);
        liceul_i = 0.0;
        liceul_i =  contcareul/(double)(M) - pow(contputeul/(double)(M),2.0);
        liceul_i = liceul_i/(double)(M);
        contliceul[i] = 1.96*exp(-r*T)*sqrt( (contcareul/(double)(M) - pow((contputeul/(double)(M)),2.0))/(double)(M) );

        errmil[i] = exp(-r*T)*putmil/(double)(M) - BS;
        licmil[i] = 1.96 * exp(-r*T)*sqrt( (carmil/(double)(M) - pow((putmil/(double)(M)),2.0))/(double)(M) );
        conterrmil[i] = exp(-r*T)*contputmil/(double)(M);
        contlicmil[i] = 1.96*exp(-r*T)*sqrt( (contcarmil/(double)(M) - pow((contputmil/(double)(M)),2.0))/(double)(M) );

        Npas[i] = N;
        N = N*2;
    }

    cout << "Erreur faible Euler" << endl;
    print(erreul);
    cout << "LIC Euler" << endl;
    print(liceul);
    cout << "Erreur faible Euler VA controle" << endl;
    print(conterreul);
    cout << "LIC Euler VA controle" << endl;
    print(contliceul);

    cout << "Erreur faible Milstein" << endl;
    print(errmil);
    cout << "LIC Milstein" << endl;
    print(licmil);
    cout << "Erreur faible Milstein VA controle" << endl;
    print(conterrmil);
    cout << "LIC Milstein VA controle" << endl;
    print(contlicmil);

    cout <<"Npas" <<endl;
    print(Npas);

}