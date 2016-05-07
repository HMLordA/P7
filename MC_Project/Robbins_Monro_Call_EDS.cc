//
//  Robbins_Monro_Call_EDS.cpp
//  MC_git_Project
//
//  Created by Nazar KOSTYUCHYK on 18/03/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK. All rights reserved.
//

#include "Robbins_Monro_Call_EDS.h"

Robbins_Monro_Call_EDS::Robbins_Monro_Call_EDS(){}

Robbins_Monro_Call_EDS::Robbins_Monro_Call_EDS(double S0, double my_T, double vol, double my_r, double K):S0(S0), vol(vol), K(K){r=my_r;T=my_T;}

double Robbins_Monro_Call_EDS::Payoff_Call(const diffusion & eds) const {
    
    return std::max(eds.back().second-K,0.0); 
 
}