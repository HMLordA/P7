//
//  Robbins_Monro_Algo.h
//  MC_Project
//
//  Created by Nazar KOSTYUCHYK on 01/03/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK. All rights reserved.
//

#ifndef Robbins_Monro_Algo_h
#define Robbins_Monro_Algo_h

#include <stdio.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

typedef boost::numeric::ublas::vector<double> Vector;

using namespace std;

double N(double d){
    
    double res;
    
    if (d>0)
        res = 0.5 * erf(d/sqrt(2)) + 0.5; // pas besoin de sigma ici car on incorpore sigma dans note d1/2
    else
        res = 0.5 + 0.5 * erf(d/sqrt(2));
    
    return res;
    
}


double Call_Price(double S0, double T, double vol, double r, double K){
    
    double d1 = log( S0 / (K*exp(-r*T)) ) / (vol * sqrt(T)) + 0.5 * vol * sqrt(T);
    double d2 = log( S0 / (K*exp(-r*T)) ) / (vol * sqrt(T)) - 0.5 * vol * sqrt(T);
    double prix = S0 * N(d1) - K*exp(-r*T) * N(d2);
    
    return prix;
    
}

template<class T, class S, double (T::*F_Payoff)(const Vector &) const, double (T::*F_Tilda_Control)(const Vector &) const>
Vector Robbins_Monro_Algo(int M, double alpha, double gamma0, Vector theta, double c, const T& Obj, S& G){

    double S1=0.0;
    double S2=0.0;
    Vector g;
    double Sreal1=0.0;
    double Sreal2=0.0;
//    double Sth1=0.0;
//    double Sth2=0.0;
//    double th=2.1963;
    
    for (int n=0; n<M; ++n){
        
        g = G();
        
        theta = theta - (gamma0/pow(n+1,alpha)) * (2*theta - g) * (pow((Obj.*F_Payoff)(g-theta),2)/(1+pow((Obj.*F_Tilda_Control)(-theta),2 * c)));

        S1 += (Obj.*F_Payoff)(g+theta) * exp(inner_prod(-theta,g) - 0.5 * inner_prod(theta,theta));
        S2 += pow((Obj.*F_Payoff)(g+theta) * exp(inner_prod(-theta,g) - 0.5 * inner_prod(theta,theta)),2);
/*
        Sth1 += (Obj.*F_Payoff)(g+th) * exp((-th * g) - 0.5 * pow(th,2));
        Sth2 += pow((Obj.*F_Payoff)(g+th) * exp((-th * g) - 0.5 * pow(th,2)),2);
*/
        Sreal1 += (Obj.*F_Payoff)(g);
        Sreal2 += pow((Obj.*F_Payoff)(g),2);
        
    }
    
    cout <<"Le theta optimal : "<< theta << endl;
    
    cout <<"L'esperance simple : "<< exp(-Obj.r*Obj.T) * Sreal1/M << endl;
//    cout <<"Vrai Prix Call : "<< Call_Price(Obj.S0,Obj.T,Obj.vol,Obj.r,Obj.K) << endl;
    cout <<"L'esperence avec changement : "<< exp(-Obj.r*Obj.T) * S1/M << endl;
//    cout <<"L'esperence avec changement thhhh: "<< exp(-Obj.r*Obj.T) * Sth1/M << endl;
    
    cout <<"L'esperance du carré simple : "<< S2/M << endl;
    cout <<"L'esperance du carré avec changement : "<< Sreal2/M << endl;
    
    // double variance = 1.0/(M-1)*(Sreal2-M*(Sreal1/M)*(Sreal1/M));
    // cout <<"La variance : "<< variance << endl;
    
    cout <<"La variance simple : "<< Sreal2/M-pow(Sreal1/M,2) << endl;
    cout <<"La variance avec changement : "<< S2/M-pow(S1/M,2) << endl;
//    cout <<"La variance avec changement thhhh: "<< Sth2/M-pow(Sth1/M,2) << endl;
    
    return theta;

}

#endif /* Robbins_Monro_Algo_h */
