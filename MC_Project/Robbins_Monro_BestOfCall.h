//
//  Robbins_Monro_BestOfCall.h
//  MC_Project
//
//  Created by Nazar KOSTYUCHYK - JC DIETRICH on 02/03/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK - JC DIETRICH. All rights reserved.
//

#ifndef Robbins_Monro_BestOfCall_h
#define Robbins_Monro_BestOfCall_h

#include <stdio.h>
#include <cmath>
#include <boost/numeric/ublas/vector.hpp>
#include "Robbins_Monro_Payoff.h"

typedef boost::numeric::ublas::vector<double> Vector;

class Robbins_Monro_BestOfCall:public Robbins_Monro_Payoff{
    
public:
    
    Robbins_Monro_BestOfCall();
    Robbins_Monro_BestOfCall(double S01, double S02, double K, double my_T, double vol1, double vol2, double my_r);
    
    virtual double Payoff(const Vector & G) const override;
    virtual double StockBS(const Vector & G) const override;
    
    double S01;
    double S02;
    double K;
    double vol1;
    double vol2;

};


#endif /* Robbins_Monro_BestOfCall_h */
