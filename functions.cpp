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

double NN(double x){
    return erfc(-x/sqrt(2))/2;
}

void print_head(vector<double> L){
    cout << "[";
    for(int i =0;i<10;i++){
        cout << L[i] <<", ";
    }
    cout <<"]" <<endl;
}

void print(vector<double> L){
    cout << "[";
    for(auto l:L){
        cout << l <<", ";
    }
    cout <<"]" <<endl;
}

double carre(double x){
    return x*x;
}

double quadruple(double x){
    return pow(x,4.0);
}

vector<double> gaussian(double mu, double sigma, size_t size){
    vector<double> G;
    default_random_engine generator;
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
    normal_distribution<double> distribution(mu,sigma);
    for(int i=0;i<size;i++){
        G.push_back(distribution(generator));
    }
    return G;
}