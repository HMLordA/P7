//
//  Processus.hpp
//  EK_MC
//
//  Created by Nazar KOSTYUCHYK on 19/02/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK. All rights reserved.
//

#ifndef Processus_h
#define Processus_h

#include <stdio.h>
#include <ostream>
#include <vector>
#include <list>
#include <math.h>
#include "Variables.h"

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___PROCESSUS___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// Classe générique on l'implemente dans ce fichier .h
template <typename T>
class Processus{
    
public:
    // state: notre valeur du processus liée a la date t : first correspond a ti second a valeur en ti
    // result_type: liste des etats du processus (liste des ti et valeurs en ti)
    // iter:
    // cst_iter:
    typedef std::pair<double, T> state;
    typedef std::list<state> result_type;
    typedef typename result_type::iterator iter;
    typedef typename result_type::const_iterator cst_iter;
    
    // constructeur qui joue le role de constructeur par default à la fois
    Processus(int size = 0) : value(size) {}
    
    // surcharge de l'operator () qui va generer le processus complet à chaque appel
    virtual result_type operator()() = 0;
    
    // getteur qui nous retourne notre processus généré
    result_type current() const { return value; }
    
    // affichage du processus qui sera a définir en dehors de la classe
    template <typename S>
    friend std::ostream& operator<<(std::ostream &o, const Processus<S> &p);
    
protected:
    // result_type: liste des etats du processus (liste des ti et valeurs en ti)
    result_type value;
};

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___PROCESSUS___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___BROWNIAN___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// Soit notre classe qui generera le processus: Mouvement Brownien
class Brownian : public Processus<double> {
    
public:
    // Implementation dans le .cc
    Brownian(int n, double T=1);
    
    // Implementation dans le .cc
    result_type operator()();
    
    // Implementation dans le .cc
    result_type affine();
    
    //
    friend class Black_scholes;
    
protected:
    int n;
    double h, T;
    Gaussian G;
};

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___BROWNIAN___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___PONT_BROWNIAN___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class Brownian_bridge : public Brownian {
    
public:
    
    // Implementation dans le .cc
    Brownian_bridge(int n, double T = 1, double B_T = 0);
    
    // Implementation dans le .cc
    result_type operator()();
    
private:
    // le point d'arrivé de notre pont en T (ZERO par default)
    double B_T;
};

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___PONT_BROWNIAN___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___BROWNIAN_BIAISE___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class Brownian_biased : public Brownian {
    
public:
    
    // Implementation dans le .cc
    Brownian_biased(int n, Var_alea<double> & X, double T = 1);
    
    // Implementation dans le .cc
    result_type operator()();

private:
    //
    Var_alea<double> & X;
};

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___BROWNIAN_BIAISE___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___BLACK & SCHOLES___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class Black_scholes : public Brownian {
    
public:
    
    Black_scholes(int n, double x0, double r, double s, double T=1);
    
    result_type operator()();
    
private:
    
    class Fun_bs : public std::unary_function<state, state> {
        
    public:
        
        Fun_bs(double x0, double r, double s);
        
        state operator()(const state & x);
        
    private:
        
        double x0, s, mu;
        
    } bs;
    
};

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___BLACK & SCHOLES___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___B & S DRIFT(t)___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
template<class S>
class BS_Drift_t: public Brownian {
    
public:
    
    BS_Drift_t(int n, double x0, double r, double s, S& theta, double T=1): Brownian(n, T), bs(x0, r, s, theta) {}
    
    result_type operator()(){
        
        Brownian::operator()();
        std::transform(value.begin(), value.end(), value.begin(), bs);
        return value;
    
    }
    
    void setNewTheta(S& thet){
    
        bs.setTheta(thet);
    }
    
private:
    
    class Fun_bs : public std::unary_function<state, state> {
        
    public:
        
        Fun_bs(double x0, double r, double s, S& thet):x0(x0), s(s), mu(r-0.5*s*s), theta(thet) {}
        
        state operator()(const state & x){
            
            return state(x.first, x0*exp((mu-theta.value(x.first))*x.first + s*x.second));
        
        }
        
        void setTheta(S& thet){
            
            theta = thet;
        }
        
    private:
        
        double x0, s, mu;
        S theta;
        
    } bs;

public:
    
    Fun_bs getBs(){
        
        return bs;
        
    }
    
};

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___B & S DRIFT(t)___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___ORNSTEIN & UHLENBECK___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class Ornstein_uhlenbeck : public Processus<double>{
    
public:
    
    Ornstein_uhlenbeck(int n, double x0, double lambda, double mu, double s, double T = 1);
    
    result_type operator()();
    
private:
    
    double x0, h, T, ret;
    Gaussian G;
};

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___ORNSTEIN & UHLENBECK___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___POISSON___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class Ppoisson : public Processus<unsigned> {
    
public:
    
    Ppoisson(double lambda, double T = 1.);
    
    result_type operator()();
    
    double getT();

protected:
    
    double T;
    Expo E;

};

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___POISSON___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___POISSON BIS___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
/*
class Ppoisson_bis : public Processus<unsigned>{
    
public:
    
    Ppoisson_bis(double lambda, double T = 1): U(0,1), P(lambda*T), T(T) {};
    
    result_type operator()() {
        
        value.clear();
        value.push_back(state(0,0));
        
        std::vector<double> vect_Uk(P());
        std::generate(vect_Uk.begin(), vect_Uk.end(), U);
        std::sort(vect_Uk.begin(), vect_Uk.end());
        
        for (int k = 0; k < vect_Uk.size(); k++) {
            
            value.push_back(state(T*vect_Uk[k], k+1));
            
        }
        
        value.push_back(state(T, P.current()));
        return value;
    };
    
protected:
    
    double T;
    Uniform U;
    Ppoisson P;
};
*/
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___POISSON BIS___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___POISSON COMPOSE___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class Ppoisson_compose : public Processus<double>{
    
public:
    
    Ppoisson_compose(double lambda, Var_alea<double> & X, double T = 1);
    
    result_type operator()();
    
protected:
    
    Ppoisson N;
    Var_alea<double> & X;
    
};

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___POISSON COMPOSE___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// Fonction qui s'occupe de l'affichage de nos processus
template <typename T>
std::ostream& operator<<(std::ostream &o, const Processus<T> &p) {
    
    
    typename Processus<T>::cst_iter i;
    for(i = p.value.begin(); i != p.value.end(); ++i)
        o << (*i).first << "\t" << (*i).second << std::endl; return o;

}

#endif /* Processus_h */
