//
//  Robbins_Monro_Algo.h
//  MC_Project
//
//  Created by Nazar KOSTYUCHYK - JC DIETRICH on 01/03/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK - JC DIETRICH. All rights reserved.
//

#ifndef Robbins_Monro_Algo_Normal_h
#define Robbins_Monro_Algo_Normal_h

#include <stdio.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/special_functions/erf.hpp>
#include "Robbins_Monro_Algo.h"

typedef boost::numeric::ublas::vector<double> Vector;

using namespace std;

template<class T, class S, double (T::*F_Payoff)(const Vector &) const>
Vector Robbins_Monro_Algo_Normal_Distrib(int M, double alpha, double gamma0, Vector theta, double c, const T& Obj, S& G){

    double S1=0.0;
    double S2=0.0;
    Vector g;
    double Sreal1=0.0;
    double Sreal2=0.0;

    int counter = 0;
    for (int n=0; n<M; ++n){
        
        g = G();
        if (n == counter)
		{
			cout <<"Theta at "<<n<<" : " <<theta[0]<<endl;
			if (n < 100)
				counter+=1;
			else
				counter += M/100;
		}
        double b = 2;
		double lambda = 0.25;
		double payoff = (Obj.*F_Payoff)(g-theta); 
		Vector d_theta = (gamma0/pow(n+1,alpha)) * exp(-lambda/2*pow(sqrt(inner_prod(theta,theta)),b)) * (pow(payoff,2)*(2*theta-g));

		theta = theta - d_theta;

        S1 += (Obj.*F_Payoff)(g+theta) * exp(inner_prod(-theta,g) - 0.5 * inner_prod(theta,theta));
        S2 += pow((Obj.*F_Payoff)(g+theta) * exp(inner_prod(-theta,g) - 0.5 * inner_prod(theta,theta)),2);

        Sreal1 += (Obj.*F_Payoff)(g);
        Sreal2 += pow((Obj.*F_Payoff)(g),2);
        
    }
    
    cout <<"Le theta optimal : "<< theta << endl;
    
    cout <<"L'esperance simple : "<< exp(-Obj.r*Obj.T) * Sreal1/M << endl;
    cout <<"Vrai Prix Call : "<< Call_Price(Obj.S0,Obj.T,Obj.vol,Obj.r,Obj.K) << endl;
    cout <<"L'esperence avec changement : "<< exp(-Obj.r*Obj.T) * S1/M << endl;
    
    cout <<"L'esperance du carre simple : "<< S2/M << endl;
    cout <<"L'esperance du carre avec changement : "<< Sreal2/M << endl;
    
    cout <<"La variance simple : "<< Sreal2/M-pow(Sreal1/M,2) << endl;
    cout <<"La variance avec changement : "<< S2/M-pow(S1/M,2) << endl;
    
    return theta;

}

#endif /* Robbins_Monro_Algo_h */
