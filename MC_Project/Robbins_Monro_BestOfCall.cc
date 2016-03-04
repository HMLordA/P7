//
//  Robbins_Monro_BestOfCall.cpp
//  MC_Project
//
//  Created by Nazar KOSTYUCHYK on 02/03/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK. All rights reserved.
//

#include "Robbins_Monro_BestOfCall.h"

Robbins_Monro_BestOfCall::Robbins_Monro_BestOfCall(){}

Robbins_Monro_BestOfCall::Robbins_Monro_BestOfCall(double S01, double S02, double K, double T,  double vol1, double vol2, double r):S01(S01), S02(S02), K(K), T(T), vol1(vol1), vol2(vol2), r(r){}


// Best of Calls solution Exacte
double Robbins_Monro_BestOfCall::Payoff_BestOfCall (const Vector & G) const {
    
    double S = StockBS_BestOfCall(G);
    
    return (S-K>0)?S-K:0;
    
}

double Robbins_Monro_BestOfCall::StockBS_BestOfCall(const Vector & G) const {
    
    double ST1 = S01*exp( (r-pow(vol1,2)/2) *T + vol1*G[0] );
    double ST2 = S02*exp( (r-pow(vol2,2)/2) *T + vol2*G[1] );
    
    return (ST1>ST2)?ST1:ST2;
    
}