//
//  Theta.h
//  MC_git_Project
//
//  Created by Nazar KOSTYUCHYK on 25/03/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK. All rights reserved.
//

#ifndef Theta_h
#define Theta_h

#include <stdio.h>
#include <vector>
// Polynome Legendre
#include <boost/math/special_functions/legendre.hpp>

using namespace std;

class Theta{
    
public:
    Theta(){}
    Theta(vector<double> & thet): theta_i(thet){}
    
    virtual double value(double t) const =0;
    
    // Getter
    vector<double> getTh() const{
        
        return theta_i;
        
    }
    
    // Setter
    void setTh(vector<double> & theta){
        
        theta_i = theta;
    }
    
    void setTheta_i(int i, double theti){
        
        theta_i[i] = theti;
    }
    
    double getTheta_i(int i){
        
        return theta_i[i];
    }
    
private:
    vector<double> theta_i;
    
};

class Theta_Legendre: public Theta{
    
public:
    
    Theta_Legendre(){}
    
    Theta_Legendre(vector<double> & thet){ setTh(thet); }
    
    double value(double t) const override {
        
        double theta=0.0;
        
        for(unsigned int j = 0; j < getTh().size(); j++){
            
            theta += getTh()[j] * boost::math::legendre_p(j, 2*t-1);
        }
        return theta;
    }
    
};

class Theta_Legendre_squared: public Theta{
    
public:
    
    Theta_Legendre_squared(){}
    
    Theta_Legendre_squared(vector<double> & thet){ setTh(thet); }
    
    double value(double t) const override {
        
        double theta=0.0;
        
        for(unsigned int j = 0; j < getTh().size(); j++){
            
            theta += getTh()[j] * boost::math::legendre_p(j, 2*t-1);
        }
        return theta*theta;
    }
    
};

#endif /* Theta_h */
