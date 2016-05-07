//
//  Robbins_Monro_Call.cpp
//  MC_Project
//
//  Created by Nazar KOSTYUCHYK - JC DIETRICH on 01/03/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK - JC DIETRICH. All rights reserved.
//

#include "Robbins_Monro_Call.h"

Robbins_Monro_Call::Robbins_Monro_Call(){}

Robbins_Monro_Call::Robbins_Monro_Call(double S0, double K, double my_T, double vol, double my_r):S0(S0), K(K), vol(vol){T=my_T;r=my_r;}

double Robbins_Monro_Call::Payoff(const Vector & G) const {
    
    double S = StockBS(G);
    
    return (S-K>0)?S-K:0; //Call
	//return (S-K>0)?1:0; //Digital
	//return (K-S>0)?K-S:0; //Put
    
}

double Robbins_Monro_Call::StockBS(const Vector & G) const {
    
    return S0 * exp((r - 0.5 * pow(vol,2)) * T + vol * sqrt(T) * G(0) );
    
}