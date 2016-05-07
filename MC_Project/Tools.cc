//
//  Tools.cc
//  MC_git_Project
//
//  Created by Nazar KOSTYUCHYK - JC DIETRICH on 25/03/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK - JC DIETRICH. All rights reserved.
//

#include "Tools.h"

double N(double d){
    
    double res;
    
    if (d>0)
        res = 0.5 * boost::math::erf(d/sqrt(2.0)) + 0.5; // pas besoin de sigma ici car on incorpore sigma dans note d1/2
    else
        res = 0.5 + 0.5 * boost::math::erf(d/sqrt(2.0));
    
    return res;
    
}


double Call_Price(double S0, double T, double vol, double r, double K){
    
    double d1 = log( S0 / (K*exp(-r*T)) ) / (vol * sqrt(T)) + 0.5 * vol * sqrt(T);
    double d2 = log( S0 / (K*exp(-r*T)) ) / (vol * sqrt(T)) - 0.5 * vol * sqrt(T);
    double prix = S0 * N(d1) - K*exp(-r*T) * N(d2);
    
    return prix;
    
}