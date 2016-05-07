//
//  Robbins_Monro_BestOfCall.cpp
//  MC_Project
//
//  Created by Nazar KOSTYUCHYK - JC DIETRICH on 02/03/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK - JC DIETRICH. All rights reserved.
//

#include "Robbins_Monro_BestOfCall.h"

Robbins_Monro_BestOfCall::Robbins_Monro_BestOfCall(){}

Robbins_Monro_BestOfCall::Robbins_Monro_BestOfCall(double S01, double S02, double K, double my_T,  double vol1, double vol2, double my_r):S01(S01), S02(S02), K(K), vol1(vol1), vol2(vol2){T=my_T;r=my_r;}


// Best of Calls solution Exacte
double Robbins_Monro_BestOfCall::Payoff(const Vector & G) const {
    
    double S = StockBS(G);
    
    return (S-K>0)?S-K:0;
    
}

double Robbins_Monro_BestOfCall::StockBS(const Vector & G) const {
    
    double ST1 = S01*exp( (r-pow(vol1,2)/2) *T + vol1*G[0] );
    double ST2 = S02*exp( (r-pow(vol2,2)/2) *T + vol2*G[1] );
    
    return (ST1>ST2)?ST1:ST2;
    
}