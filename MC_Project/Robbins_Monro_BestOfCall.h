//
//  Robbins_Monro_BestOfCall.h
//  MC_Project
//
//  Created by Nazar KOSTYUCHYK on 02/03/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK. All rights reserved.
//

#ifndef Robbins_Monro_BestOfCall_h
#define Robbins_Monro_BestOfCall_h

#include <stdio.h>
#include <cmath>
#include <boost/numeric/ublas/vector.hpp>

typedef boost::numeric::ublas::vector<double> Vector;

class Robbins_Monro_BestOfCall{
    
public:
    
    Robbins_Monro_BestOfCall();
    Robbins_Monro_BestOfCall(double S01, double S02, double K, double T, double vol1, double vol2, double r);
    
    double Payoff_BestOfCall(const Vector & G) const;
    double StockBS_BestOfCall(const Vector & G) const;
    
    double S01;
    double S02;
    double K;
    double T;
    double vol1;
    double vol2;
    double r;

};


#endif /* Robbins_Monro_BestOfCall_h */
