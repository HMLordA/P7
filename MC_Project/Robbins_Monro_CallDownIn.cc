//
//  Robbins_Monro_CallDownIn.cpp
//  MC_git_Project
//
//  Created by Nazar KOSTYUCHYK on 18/03/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK. All rights reserved.
//

#include "Robbins_Monro_CallDownIn.h"

Robbins_Monro_CallDownIn::Robbins_Monro_CallDownIn(){}

Robbins_Monro_CallDownIn::Robbins_Monro_CallDownIn(double S0, double T, double vol, double r, double L, double K):S0(S0), T(T), vol(vol), r(r), L(L), K(K){}

double Robbins_Monro_CallDownIn::Payoff_Call(const diffusion & eds) const {
    
    bool activation=false;

    for(cst_iter j = eds.begin(); j != eds.end(); ++j)
    {
        
        if((*j).second<=L) {
            activation = true;
            break;}

    }
    
    if(activation){return std::max(eds.back().second-K,0.0); }
    
    else
        return 0.0;
    
}