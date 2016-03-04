//
//  Robbins_Monro_Call.cpp
//  MC_Project
//
//  Created by Nazar KOSTYUCHYK on 01/03/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK. All rights reserved.
//

#include "Robbins_Monro_Call.h"

Robbins_Monro_Call::Robbins_Monro_Call(){}

Robbins_Monro_Call::Robbins_Monro_Call(double S0, double K, double T, double vol, double r):S0(S0), K(K), T(T), vol(vol), r(r){}

double Robbins_Monro_Call::Payoff_Call(const Vector & G) const {
    
    double S = StockBS(G);
    
    return (S-K>0)?S-K:0;
    
}

double Robbins_Monro_Call::StockBS(const Vector & G) const {
    
    return S0 * exp((r - 0.5 * pow(vol,2)) * T + vol * sqrt(T) * G(0) );
    
}