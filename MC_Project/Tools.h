//
//  Tools.h
//  MC_git_Project
//
//  Created by Nazar KOSTYUCHYK - JC DIETRICH on 25/03/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK - JC DIETRICH. All rights reserved.
//

#ifndef Tools_h
#define Tools_h

#include <stdio.h>
#include <boost/math/special_functions/erf.hpp>
#include "Theta.h"
#include <vector>
#include <list>

// Polynome Legendre
#include <boost/math/special_functions/legendre.hpp>

using namespace std;

double N(double d);

double Call_Price(double S0, double T, double vol, double r, double K);

class LegendreCarre{

public:
    
    LegendreCarre(int degre): l(degre) {}
    
    double value(double x) const{        
        return pow(boost::math::legendre_p(l, x),2);     
    }
    
private:
    int l;

};



// Integration de fonctions algo trouvé sur BOOST!!!
template<typename value_type, typename function_type>
inline value_type integral2P(const value_type a,
                             const value_type b,
                             const int l,
                             const value_type tol,
                             function_type func)
{
    unsigned n = 1U;
    
    value_type h = (b - a);
    value_type I = (func(l,a) + func(l,b)) * (h / 2);
    
    for(unsigned k = 0U; k < 8U; k++)
    {
        h /= 2;
        
        value_type sum(0);
        for(unsigned j = 1U; j <= n; j++)
        {
            sum += func(l,a + (value_type((j * 2) - 1) * h));
        }
        
        const value_type I0 = I;
        I = (I / 2) + (h * sum);
        
        const value_type ratio     = I0 / I;
        const value_type delta     = ratio - 1;
        const value_type delta_abs = ((delta < 0) ? -delta : delta);
        
        if((k > 1U) && (delta_abs < tol))
        {
            break;
        }
        
        n *= 2U;
    }
    
    return I;
}

// Integration de fonctions algo trouvé sur BOOST!!!
template<typename value_type>
inline value_type integral1P(const value_type a,
                             const value_type b,
                             const value_type tol,
                             const Theta& th,
							 unsigned int iterations = 8U)
{
    unsigned n = 1U;
    
    value_type h = (b - a);

    value_type I = (th.value(a) + th.value(b)) * (h / 2);
    
    for(unsigned k = 0U; k < iterations; k++)
    {
        h /= 2;
        
        value_type sum(0);
        for(unsigned j = 1U; j <= n; j++)
        {
            sum += th.value(a + (value_type((j * 2) - 1) * h));
        }
        
        const value_type I0 = I;
        I = (I / 2) + (h * sum);

		const value_type ratio     = I0 / I;
        const value_type delta     = ratio - 1;
        const value_type delta_abs = ((delta < 0) ? -delta : delta);
        
        if((k > 1U) && (delta_abs < tol))
        {
            break;
        }
        
        n *= 2U;
    }
    
    return I;
}

// Integration de fonctions algo trouvé sur BOOST!!!
template<typename value_type, typename func>
inline value_type integralSTO(
                             const list<pair<value_type,value_type>> proc,
                             const func & th)
{
    
    double res = 0.0;
    auto list_itr_plus_1=proc.begin();
    
    list_itr_plus_1++;
    for (auto list_itr=proc.begin();list_itr_plus_1!=proc.end();++list_itr)
    {
        res += (th.value(list_itr->first)) * (list_itr_plus_1->second-list_itr->second);
		list_itr_plus_1++;
    }
    
    return res;
}

#endif /* Tools_h */
