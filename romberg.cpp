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
    double T = 1.0;
    double sig = 0.2;
    double r = 0.05;
    double S0 = 100.0;
    double K = 110.0;
    int N = 1;
    int R = 10;
    int M = 100000;
    int RM = R*M;
    int P = 4;

    vector<double> erreul(P,0.0);
    vector<double> liceul(P,0.0);
    vector<double> errmil(P,0.0);
    vector<double> licmil(P,0.0);
    vector<double> Npas(P,0.0);
    vector<double> S(M,S0);
    vector<double> Se(M,S0);
    vector<double> S2e(M,S0);
    vector<double> Sm(M,S0);
    vector<double> S2m(M,S0);

    for(int i = 0; i<P; i++){
        double dt = T/(double)(N);
        double sigdt=sig*sqrt(dt);
        double rdt=r*dt;
        double der=rdt-sigdt*sigdt*0.5;

        double payeul = 0.0;
        double careul = 0.0;
        double paymil = 0.0;
        double carmil = 0.0;

        for(int j = 0; j<R; j++){
            vector<double> S_init(M,S0);
            S = S_init;
            Se = S_init;
            S2e = S_init;
            Sm = S_init;
            S2m = S_init;

            for(int k = 0; k<N; k++){
                vector<double> g1 = gaussian(0.0,1.0,M);
                for_each(g1.begin(),g1.end(),
                    [](double& x){ x = x/sqrt(2); } );
                vector<double> g2 = gaussian(0.0,1.0,M);
                for_each(g2.begin(),g2.end(),
                    [](double& x){ x = x/sqrt(2); } );
                vector<double> g(M);
                transform(g1.begin(),g1.end(),g2.begin(),g.begin(),
                    [](double x, double y){ return x+y ; } );
            
                // EDS Solution
                transform(g.begin(),g.end(),S.begin(),S.begin(),
                    [&](double x,double y){ return y*( exp( sigdt*x + der ) ); } );
            
                // Euler
                transform(g.begin(),g.end(),Se.begin(),Se.begin(),
                    [&](double x,double y){ return y*(1 + sigdt*x + rdt); } );
                /*    
                transform(g1.begin(),g1.end(),S2e.begin(),S2e.begin(),
                    [&](double x, double y){ return y*(1 + sigdt*x + rdt*0.5) ; } );
                transform(g2.begin(),g2.end(),S2e.begin(),S2e.begin(),
                    [&](double x, double y){ return y*(1 + sigdt*x + rdt*0.5) ; } );*/
            
                // Milstein
            
                transform(g.begin(),g.end(),Sm.begin(),Sm.begin(),
                    [&](double x,double y){ return y*(1 + sigdt*x + pow(sigdt*x,2.0)*0.5 + der); } );
                /*
                transform(g1.begin(),g1.end(),S2m.begin(),S2m.begin(),
                    [&](double x,double y){ return y*(1 + sigdt*x + sigdt*x*sigdt*x*0.5 + rdt*0.5); } );
                transform(g2.begin(),g2.end(),S2m.begin(),S2m.begin(),
                    [&](double x,double y){ return y*(1 + sigdt*x + sigdt*x*sigdt*x*0.5 + rdt*0.5); } );*/
                
                for(int l=0;l<M;l++){
                    double x1 = g1[l];
                    double x2 = g2[l];
                    // Euler
                    S2e[l] = S2e[l]*(1+sigdt*x1+rdt*0.5)*(1+sigdt*x2+rdt*0.5);
                    // Milstein
                    S2m[l] = S2m[l]*(1+sigdt*x1+pow(sigdt*x1,2.0)*0.5+der*0.5)*(1+sigdt*x2+pow(sigdt*x2,2.0)*0.5+der*0.5);
                }
            }
            // contribution euler
            vector<double> eul;
            for_each(S2e.begin(),S2e.end(),
                [&](double x){ eul.push_back( (K>x) ? (2*(K-x)):0.0 ); } );
            vector<double> diff_eul;
            for_each(Se.begin(),Se.end(),
                [&](double x){ diff_eul.push_back( (K>x) ? (K-x):0.0 ); } );
            transform(eul.begin(),eul.end(),diff_eul.begin(),eul.begin(),
                [&](double x, double y){ return x-y ; } );
            diff_eul.clear();
            for_each(S.begin(),S.end(),
                [&](double x){ diff_eul.push_back( (K>x) ? (K-x):0.0 ); } );
            transform(eul.begin(),eul.end(),diff_eul.begin(),eul.begin(),
                [&](double x, double y){ return x-y ; } );
            payeul += (accumulate(eul.begin(),eul.end(),0.0))/(double)(M);
            double sum_squared = 0.0;
            for_each(eul.begin(),eul.end(),
                [&](double x){ sum_squared += x*x; } );
            careul = careul + sum_squared/(double)(M);

            // contribution milstein
            vector<double> mil;
            for_each(S2m.begin(),S2m.end(),
                [&](double x){ mil.push_back( (K>x) ? (2*(K-x)):0.0 ); } );
            vector<double> diff_mil;
            for_each(Sm.begin(),Sm.end(),
                [&](double x){ diff_mil.push_back( (K>x) ? (K-x):0.0 ); } );
            transform(mil.begin(),mil.end(),diff_mil.begin(),mil.begin(),
                [&](double x, double y){ return x-y ; } );
            diff_mil.clear();
            for_each(S.begin(),S.end(),
                [&](double x){ diff_mil.push_back( (K>x) ? (K-x):0.0 ); } );
            transform(mil.begin(),mil.end(),diff_mil.begin(),mil.begin(),
                [&](double x, double y){ return x-y ; } );
            paymil += (accumulate(mil.begin(),mil.end(),0.0))/(double)(M);
            sum_squared = 0.0;
            for_each(mil.begin(),mil.end(),
                [&](double x){ sum_squared += x*x; } );
            carmil = carmil + sum_squared/(double)(M);
        }

        erreul[i] = exp(-r*T)*payeul/(double)(R);
        double cal_eul = (payeul/(double)(R))*(payeul/(double)(R));
        cal_eul = (careul/(double)(R) - cal_eul);
        liceul[i] = 1.96*exp(-r*T)*sqrt(cal_eul/(double)(RM));

        errmil[i] = exp(-r*T)*paymil/(double)(R);
        double cal_mil = (paymil/(double)(R))*(paymil/(double)(R));
        cal_mil =( carmil/(double)(R) - cal_mil);
        licmil[i] = 1.96*exp(-r*T)*sqrt(cal_mil/(double)(RM));
        
        Npas[i] = N;
        N = N*2;
    }

    cout << "Romberg Euler" <<endl;
    print(erreul);
    cout << "LIC Euler" <<endl;
    print(liceul);
    cout << "Romberg Milstein" <<endl;
    print(errmil);
    cout << "LIC Milstein" <<endl;
    print(licmil);

    cout << "Npas" << endl;
    print(Npas);


}
