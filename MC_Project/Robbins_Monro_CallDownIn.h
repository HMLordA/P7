//
//  Robbins_Monro_CallDownIn.h
//  MC_git_Project
//
//  Created by Nazar KOSTYUCHYK on 18/03/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK. All rights reserved.
//

#ifndef Robbins_Monro_CallDownIn_h
#define Robbins_Monro_CallDownIn_h

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <list>

class Robbins_Monro_CallDownIn{
    
public:
    
    typedef std::pair<double, double> state;
    typedef std::list<state> diffusion;
    typedef std::list<state>::iterator iter;
    typedef std::list<state>::const_iterator cst_iter;
    
    
    Robbins_Monro_CallDownIn();
    Robbins_Monro_CallDownIn(double S0, double T, double vol, double r, double L, double K);
    
    double Payoff_Call(const diffusion & eds) const;
    
    double S0;
    double T;
    double vol;
    double r;
    double L;
    double K;
    
};


#endif /* Robbins_Monro_CallDownIn_h */
