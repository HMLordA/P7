//
//  Robbins_Monro_Call.h
//  MC_Project
//
//  Created by Nazar KOSTYUCHYK on 01/03/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK. All rights reserved.
//

#ifndef Robbins_Monro_Call_h
#define Robbins_Monro_Call_h

#include <stdio.h>
#include <cmath>
#include <boost/numeric/ublas/vector.hpp>
#include "Robbins_Monro_Payoff.h"

class Robbins_Monro_Call:public Robbins_Monro_Payoff{

public:
    
    Robbins_Monro_Call();
    Robbins_Monro_Call(double S0, double K, double my_T, double vol, double my_r);
    
    virtual double Payoff(const Vector & G) const override;
    virtual double StockBS(const Vector &  G) const override;
    
    double S0;
    double K;
    //double T;
    double vol;
    //double r;

};

#endif /* Robbins_Monro_Call_h */

