//
//  Robbins_Monro_CallDownIn.h
//  MC_git_Project
//
//  Created by Nazar KOSTYUCHYK - JC DIETRICH on 18/03/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK - JC DIETRICH. All rights reserved.
//

#ifndef Robbins_Monro_CallDownIn_h
#define Robbins_Monro_CallDownIn_h

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <list>
#include "Robbins_Monro_Payoff_EDS.h"

class Robbins_Monro_CallDownIn:public Robbins_Monro_Payoff_EDS{
    
public:

    typedef std::list<state>::iterator iter;
    typedef std::list<state>::const_iterator cst_iter;
    
    
    Robbins_Monro_CallDownIn();
    Robbins_Monro_CallDownIn(double S0, double my_T, double vol, double my_r, double L, double K);
    
    virtual double Payoff_Call(const diffusion & eds) const override;
    
    double S0;
    double vol;
    double L;
    double K;
    
};


#endif /* Robbins_Monro_CallDownIn_h */
