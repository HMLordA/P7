//
//  main.cpp
//  MC_Project
//
//  Created by Nazar KOSTYUCHYK on 01/03/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK. All rights reserved.
//

#include <iostream>
#include "Variables.h"
#include <cmath>
#include "Robbins_Monro_Call.h"
#include "Robbins_Monro_BestOfCall.h"
#include "Robbins_Monro_Algo.h"
#include "Robbins_Monro_Algo_Normal_Distrib.h"

using namespace std;
/*
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
*/
double StockBS(double S0, double T, double vol, double r, double G){

    return S0 * exp((r - 0.5 * pow(vol,2)) * T + vol * sqrt(T) * G );
    
}

double Payoff_call (double S0, double T, double vol, double r, double G, double K){
    
    double S = StockBS(S0,T,vol,r,G);
    
    return (S-K>0)?S-K:0;

}

// Best of Calls solution Exacte
double Payoff_maxcall (double L, double G1, double G2, double S1, double S2, double T, double vol1, double vol2, double r){
    double ST1 = S1*exp( (r-pow(vol1,2)/2) *T + vol1*G1 );
    double ST2 = S2*exp( (r-pow(vol2,2)/2) *T + vol2*G2 );
    double payoff = max(max(ST1, ST2)-L, 0.0);
    return payoff;
}

int main(int argc, const char * argv[]) {
    
	const int NB_ASSETS = 1;

    Matrix m(NB_ASSETS,NB_ASSETS,0);
    for (unsigned i = 0; i < m.size1 (); ++ i){
        m (i, i) = 1;
    }
    
    
    double S01=100.0;
    double S02=105.0;
    double K=140.0;
    double T=1.0;
    double vol1=0.2;
    double vol2=0.15;
    double r=0.05;
    
    int M=10000000;
    double alpha=0.500001;
    double gamma0=1.0;
    Vector theta(NB_ASSETS,0);
    double c=1.0;

	Vector v(NB_ASSETS,0);
    
    //Robbins_Monro_BestOfCall rmb(S01, S02, K, T, vol1, vol2, r);
	Robbins_Monro_Call rmc(S01,K,T,vol1,r);
    Gaussian_Vector G(v,m,NB_ASSETS);

    //Vector thet = Robbins_Monro_Algo<Robbins_Monro_BestOfCall, Gaussian_Vector, &Robbins_Monro_BestOfCall::Payoff_BestOfCall, &Robbins_Monro_BestOfCall::StockBS_BestOfCall>(M, alpha, gamma0, theta, c, rmb, G);
    Vector thet = Robbins_Monro_Algo<Robbins_Monro_Call, Gaussian_Vector, &Robbins_Monro_Call::Payoff_Call, &Robbins_Monro_Call::StockBS>(M, alpha, gamma0, theta, c, rmc, G);
    //Vector thet = Robbins_Monro_Algo_Normal_Distrib<Robbins_Monro_Call, Gaussian_Vector, &Robbins_Monro_Call::Payoff_Call>(M, alpha, gamma0, theta, c, rmc, G);


/*
    //Identity_Matrix m (3);
    Matrix m(3,3,0);
    for (unsigned i = 0; i < m.size1 (); ++ i){
            m (i, i) = 1;
    }
    
    cout << m << endl;
    
    Vector v(3,0);
    
    Gaussian_Vector G(v,m,3);
    
    double S01=100.0;
    double S02=105.0;
    double K=110.0;
    double T=1.0;
    double vol1=0.2;
    double vol2=0.15;
    double r=0.05;
    
    Robbins_Monro_BestOfCall(S01, S02, K, T, vol1, vol2, r);
    
    v=G();
    
    for(int i=0;i<v.size();++i) {
        std::cout << v(i) << endl;
    }
*/
    
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
/*
 double a[] = { 0.11, 0.12, 0.13,
 0.21, 0.22, 0.23 };
 
 double b[] = { 1011, 1012,
 1021, 1022,
 1031, 1032 };
 
 double c[] = { 0.00, 0.00,
 0.00, 0.00 };
 
 gsl_matrix_view A = gsl_matrix_view_array(a, 2, 3);
 gsl_matrix_view B = gsl_matrix_view_array(b, 3, 2);
 gsl_matrix_view C = gsl_matrix_view_array(c, 2, 2);
 
 // Compute C = A B
 gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
 1.0, &A.matrix, &B.matrix,0.0, &C.matrix);
 
 printf ("[ %g, %g\n", c[0], c[1]);
 printf (" %g, %g ]\n", c[2], c[3]);
 
*/
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
/*
    Matrix m(2,2,0);
    for (unsigned i = 0; i < m.size1 (); ++ i){
        m (i, i) = 1;
    }
    
    Vector v(2,0);
    
    double S01=100.0;
    double S02=105.0;
    double K=110.0;
    double T=1.0;
    double vol1=0.2;
    double vol2=0.15;
    double r=0.05;
    
    int M=100000;
    double alpha=0.500001;
    double gamma0=1.0;
    Vector theta(2,0);
    double c=1.0;
    
    Robbins_Monro_BestOfCall rmb(S01, S02, K, T, vol1, vol2, r);
    Gaussian_Vector G(v,m,2);
    
    Vector thet = Robbins_Monro_Algo<Robbins_Monro_BestOfCall, Gaussian_Vector, &Robbins_Monro_BestOfCall::Payoff_BestOfCall, &Robbins_Monro_BestOfCall::StockBS_BestOfCall>(M, alpha, gamma0, theta, c, rmb, G);
 */
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
/*
 
 Matrix m(1,1,0);
 for (unsigned i = 0; i < m.size1 (); ++ i){
 m (i, i) = 1;
 }
 
 Vector v(1,0.0);
 
 double S0=100.0;
 double K=110.0;
 double T=1.0;
 double vol=0.2;
 double r=0.05;
 
 int M=100000;
 double alpha=0.500001;
 double gamma0=1.0;
 Vector theta(1,0);
 double c=1.0;
 
 Robbins_Monro_Call rmc(S0,K,T,vol,r);
 Gaussian_Vector G(v,m,1);
 
 Vector thet = Robbins_Monro_Algo<Robbins_Monro_Call, Gaussian_Vector, &Robbins_Monro_Call::Payoff_Call, &Robbins_Monro_Call::StockBS>(M, alpha, gamma0, theta, c, rmc, G);
 
*/
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
/*
    Gaussian G(0.0,1.0);
    double S0=100.0;
    double K=140.0;
    double T=1.0;
    double vol=0.2;
    double r=0.05;
    
    double alpha=0.500001;
    int M=10000000;
    
    double S1=0.0;
    double S2=0.0;
    double theta=0.0;
    double gamma0=1.0;
    double g=0.0;
    double Sreal1=0.0;
    double Sreal2=0.0;
    double Sth1=0.0;
    double Sth2=0.0;
    double th=2.1963;
    
    for (int n=0; n<M; ++n){
        
        g = G();
        
        S1 += Payoff_call(S0,T,vol,r,g+theta,K) * exp((-theta * g) - 0.5 * pow(theta,2));
        S2 += pow(Payoff_call(S0,T,vol,r,g+theta,K) * exp((-theta * g) - 0.5 * pow(theta,2)),2);
       
        theta = theta - (gamma0/pow(n+1,alpha)) * (2*theta - g) * (pow(Payoff_call(S0,T,vol,r,g-theta,K),2)/(1+pow(StockBS(S0,T,vol,r,-theta),2)));
        
        Sth1 += Payoff_call(S0,T,vol,r,g+th,K) * exp((-th * g) - 0.5 * pow(th,2));
        Sth2 += pow(Payoff_call(S0,T,vol,r,g+th,K) * exp((-th * g) - 0.5 * pow(th,2)),2);
        
        Sreal1 += Payoff_call(S0,T,vol,r,g,K);
        Sreal2 += pow(Payoff_call(S0,T,vol,r,g,K),2);
        
    }
    
    cout <<"Le theta optimal : "<< theta << endl;
    
    cout <<"L'esperance simple : "<< exp(-r*T) * Sreal1/M << endl;
        cout <<"Vrai Prix Call : "<< Call_Price(S0,T,vol,r,K) << endl;
    cout <<"L'esperence avec changement : "<< exp(-r*T) * S1/M << endl;
    cout <<"L'esperence avec changement thhhh: "<< exp(-r*T) * Sth1/M << endl;


    cout <<"L'esperance du carré simple : "<< S2/M << endl;
    cout <<"L'esperance du carré avec changement : "<< Sreal2/M << endl;
    
    // double variance = 1.0/(M-1)*(Sreal2-M*(Sreal1/M)*(Sreal1/M));
    // cout <<"La variance : "<< variance << endl;
    
    cout <<"La variance simple : "<< Sreal2/M-pow(Sreal1/M,2) << endl;
    cout <<"La variance avec changement : "<< S2/M-pow(S1/M,2) << endl;
    cout <<"La variance avec changement thhhh: "<< Sth2/M-pow(Sth1/M,2) << endl;
*/
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    return 0;
    
}
