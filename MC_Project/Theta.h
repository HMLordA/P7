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
    
	virtual double value(double t) const { throw 1; return 0.0; }
	
	virtual Theta* copy() {return 0;}
    
    // Getter
    vector<double> getTh() const{
        
        return theta_i;
        
    }
    
    // Setter
    virtual void setTh(vector<double> & theta){
        
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

	virtual Theta_Legendre* copy() override {return new Theta_Legendre(this->getTh());}
    
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

	virtual Theta_Legendre_squared* copy() override {return new Theta_Legendre_squared(this->getTh());}
    
    double value(double t) const override {
        
        double theta=0.0;
        
        for(unsigned int j = 0; j < getTh().size(); j++){
            
            theta += getTh()[j] * boost::math::legendre_p(j, 2*t-1);
        }
        return theta*theta;
    }
    
};


double Phi_n_k(double n, double k, double t);

class Theta_Haar: public Theta{
    
public:
    
    Theta_Haar():n(0.0){}
    
    Theta_Haar(vector<double> & thet){ setTh(thet); }

	virtual void setTh(vector<double> & theta) override {
		Theta::setTh(theta);
		//n = int(log((double)theta.size())/log(2.0)); 
		n = int(log(double(theta.size()-1)+1.0)/log(2.0))-1;
	}

	virtual Theta_Haar* copy() override {return new Theta_Haar(this->getTh());}
    
    double value(double t) const override {
        
        double theta=0.0;
		//double N = getTh().size();
		theta += getTh()[0]; //constant first
		int counter = 1;
		for (int current_n=0;current_n<=n;current_n++)
		{
			double j = pow(2.0,current_n)-1; 
			for(int current_j=0; current_j<=j; current_j++){ 

				theta += getTh()[counter] * Phi_n_k(current_n,current_j,t);
				counter++;
			}
		}
        return theta;
    }

private:
	double n;
    
};

class Theta_Haar_squared: public Theta{
    
public:
    
    Theta_Haar_squared():n(0.0){}
    
    Theta_Haar_squared(vector<double> & thet){ setTh(thet); }

	virtual Theta_Haar_squared* copy() override {return new Theta_Haar_squared(this->getTh());}

	virtual void setTh(vector<double> & theta) override {
		Theta::setTh(theta);
		//n = int(log((double)theta.size())/log(2.0)); 
		n = int(log(double(theta.size()-1)+1.0)/log(2.0))-1;
	}
    
    double value(double t) const override {
        
        double theta=0.0;
		double N = getTh().size();
		theta += getTh()[0];
		int counter = 1;
		for (int current_n=0;current_n<=n;current_n++)
		{
			double j = pow(2.0,current_n)-1; 
			for(int current_j=0; current_j<=j; current_j++){   
				//for(unsigned int j = 0; j < getTh().size(); j++){
				theta += getTh()[counter] * Phi_n_k(current_n,current_j,t);
				counter++;
			}
		}
        return theta*theta;
    }

private:
	double n;
    
};

class HaarCarre{

public:
    
    HaarCarre(int my_n,int my_k): n(my_n),k(my_k) {}
    
    double value(double x) const{

		double phi = 1;
		if (n >= 0) //non constant polynom
			phi = Phi_n_k(n,k,x);
		else 
			phi = phi*(x>=0)*(x<=1);
        return phi/* * phi*/;       
    }
    
private:
    int n;
	int k;

};

#endif /* Theta_h */
