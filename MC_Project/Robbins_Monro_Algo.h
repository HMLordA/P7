//
//  Robbins_Monro_Algo.h
//  MC_Project
//
//  Created by Nazar KOSTYUCHYK on 01/03/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK. All rights reserved.
//

#ifndef Robbins_Monro_Algo_h
#define Robbins_Monro_Algo_h

#include <stdio.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <vector>

// Polynome Legendre
#include <boost/math/special_functions/legendre.hpp>


typedef boost::numeric::ublas::vector<double> Vector;

using namespace std;

double N(double d){
    
    double res;
    
    if (d>0)
        res = 0.5 * boost::math::erf(d/sqrt(2.0)) + 0.5; // pas besoin de sigma ici car on incorpore sigma dans note d1/2
    else
        res = 0.5 + 0.5 * boost::math::erf(d/sqrt(2.0));
    
    return res;
    
}


double Call_Price(double S0, double T, double vol, double r, double K){
    
    double d1 = log( S0 / (K*exp(-r*T)) ) / (vol * sqrt(T)) + 0.5 * vol * sqrt(T);
    double d2 = log( S0 / (K*exp(-r*T)) ) / (vol * sqrt(T)) - 0.5 * vol * sqrt(T);
    double prix = S0 * N(d1) - K*exp(-r*T) * N(d2);
    
    return prix;
    
}

// Integration de fonctions algo trouvé sur BOOST!!!
template<typename value_type, typename function_type>
inline value_type integral(const value_type a,
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

class Theta{
    
public:
    Theta(){}
    Theta(vector<double> & thet): theta_i(thet){}
    
    virtual double value(double t)=0;

    // Getter
    vector<double> getTh(){
        
        return theta_i;
    
    }
    
    // Setter
    void setTh(vector<double> & theta){
    
        theta_i = theta;
    }
    
    void setTheta_i(int i, double theti){
        
        theta_i[i] = theti;
    }
    
private:
    vector<double> theta_i;

};

class Theta_Legendre: public Theta{

public:
    
    Theta_Legendre(){}
    
    Theta_Legendre(vector<double> & thet){ setTh(thet); }
    
    double value(double t) override {
        
        double theta=0.0;
        
        for(int j = 0; j < getTh().size(); j++){
            
            theta += getTh()[j] * boost::math::legendre_p(j, 2*t-1);
        }
        return theta;
    }
    
};

template<class T, class S, double (T::*F_Payoff)(const Vector &) const, double (T::*F_Tilda_Control)(const Vector &) const>
Vector Robbins_Monro_Algo(int M, double alpha, double gamma0, Vector theta, double c, const T& Obj, S& G){

    double S1=0.0;
    double S2=0.0;
    Vector g;
    double Sreal1=0.0;
    double Sreal2=0.0;
//    double Sth1=0.0;
//    double Sth2=0.0;
//    double th=2.1963;
    int counter = 0;
    for (int n=0; n<M; ++n){
        
        g = G();
        if (n == counter)
		{
			cout <<"Theta at "<<n<<" : " <<theta[0]<<endl;
			if (n < 100)
				counter+=1;
			else
				counter += M/100;
		}
		Vector d_theta = ((gamma0/pow(n+1,alpha))) * (2*theta - g) * (pow((Obj.*F_Payoff)(g-theta),2)/(1+pow((Obj.*F_Tilda_Control)(-theta),2 * c)));
        theta = theta - d_theta;

        S1 += (Obj.*F_Payoff)(g+theta) * exp(inner_prod(-theta,g) - 0.5 * inner_prod(theta,theta));
        S2 += pow((Obj.*F_Payoff)(g+theta) * exp(inner_prod(-theta,g) - 0.5 * inner_prod(theta,theta)),2);
/*
        Sth1 += (Obj.*F_Payoff)(g+th) * exp((-th * g) - 0.5 * pow(th,2));
        Sth2 += pow((Obj.*F_Payoff)(g+th) * exp((-th * g) - 0.5 * pow(th,2)),2);
*/
        Sreal1 += (Obj.*F_Payoff)(g);
        Sreal2 += pow((Obj.*F_Payoff)(g),2);
        
    }
    
    cout <<"Le theta optimal : "<< theta << endl;
    
    cout <<"L'esperance simple : "<< exp(-Obj.r*Obj.T) * Sreal1/M << endl;
//    cout <<"Vrai Prix Call : "<< Call_Price(Obj.S0,Obj.T,Obj.vol,Obj.r,Obj.K) << endl;
    cout <<"L'esperence avec changement : "<< exp(-Obj.r*Obj.T) * S1/M << endl;
//    cout <<"L'esperence avec changement thhhh: "<< exp(-Obj.r*Obj.T) * Sth1/M << endl;
    
    cout <<"L'esperance du carre simple : "<< S2/M << endl;
    cout <<"L'esperance du carre avec changement : "<< Sreal2/M << endl;
    
    // double variance = 1.0/(M-1)*(Sreal2-M*(Sreal1/M)*(Sreal1/M));
    // cout <<"La variance : "<< variance << endl;
    
    cout <<"La variance simple : "<< Sreal2/M-pow(Sreal1/M,2) << endl;
    cout <<"La variance avec changement : "<< S2/M-pow(S1/M,2) << endl;
//    cout <<"La variance avec changement thhhh: "<< Sth2/M-pow(Sth1/M,2) << endl;
    
    return theta;

}

// joue le rôle de notre polynome
double legendreCarre(int l, double x){
    
    return pow(boost::math::legendre_p(l, x),2);
    
}

template<class T, class S, class SS, class U, double (T::*F_Payoff)(const std::list<std::pair<double,double>> &) const>
void Robbins_Monro_SDE_Algo(int M, double alpha, double gamma0, Theta_Legendre& theta, double c, const T& Obj, S& EDS1, SS& EDS2, U& G){
    
    double S1=0.0;
    double S2=0.0;
    double Sreal1=0.0;
    double Sreal2=0.0;
    double g;
    vector<U> gaussian;
    
    for(int j=0; j<theta.getTh().size(); j++){
        
        gaussian.push_back(U(0.0,integral(0.0, 1.0, j, 0.01,legendreCarre)));
        
    }
    
    for (int n=0; n<M; ++n){
        
        EDS1.setNewTheta(theta);
        
        EDS1();
        EDS2();
        
        // Mise a jour des thetas
        for(int i=0; i<theta.getTh().size(); i++){
            
            g = gaussian[i]();
            
            theta.setTheta_i( i, theta.getTh()[i] - (gamma0/(pow(n+1,alpha)+0.0001)) * ( pow((Obj.*F_Payoff)(EDS1.current()),2)*( 2*theta.getTh()[i] - g ) ) );
            
        }
        //_______________________
        
        S1 += (Obj.*F_Payoff)(EDS1.current());
        S2 += pow((Obj.*F_Payoff)(EDS1.current()),2);
        
        Sreal1 += (Obj.*F_Payoff)(EDS2.current());
        Sreal2 += pow((Obj.*F_Payoff)(EDS2.current()),2);
        
    }
    
    cout <<"L'esperance simple : "<< exp(-Obj.r*Obj.T) * Sreal1/M << endl;
    cout <<"L'esperence avec changement : "<< exp(-Obj.r*Obj.T) * S1/M << endl;
    
    cout <<"L'esperance du carre simple : "<< Sreal2/M << endl;
    cout <<"L'esperance du carre avec changement : "<< S2/M << endl;
    
    cout <<"La variance simple : "<< Sreal2/M-pow(Sreal1/M,2) << endl;
    cout <<"La variance avec changement : "<< S2/M-pow(S1/M,2) << endl;
    
}

#endif /* Robbins_Monro_Algo_h */
