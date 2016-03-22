//
//  Robbins_Monro_Call_EDS.cpp
//  MC_git_Project
//
//  Created by Nazar KOSTYUCHYK on 18/03/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK. All rights reserved.
//

#include "Robbins_Monro_Call_EDS.h"

Robbins_Monro_Call_EDS::Robbins_Monro_Call_EDS(){}

Robbins_Monro_Call_EDS::Robbins_Monro_Call_EDS(double S0, double T, double vol, double r, double K):S0(S0), T(T), vol(vol), r(r), K(K){}

double Robbins_Monro_Call_EDS::Payoff_Call(const diffusion & eds) const {
    
    return std::max(eds.back().second-K,0.0); 
 
}