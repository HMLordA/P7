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

class Robbins_Monro_Call{

public:
    
    Robbins_Monro_Call();
    Robbins_Monro_Call(double S0, double K, double T, double vol, double r);
    
    double Payoff_Call(double G) const;
    double StockBS(double G) const;
    
    double S0;
    double K;
    double T;
    double vol;
    double r;

};

#endif /* Robbins_Monro_Call_h */

