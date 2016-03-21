//
//  Processus.cpp
//  EK_MC
//
//  Created by Nazar KOSTYUCHYK on 19/02/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK. All rights reserved.
//

#include "Processus.h"

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___BROWNIAN___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// constructeur va construire le processus de taille 2^n + 1, et de pas de temps T/2^n
Brownian::Brownian(int n, double T): Processus<double>(pow(2.0,n)+1), n(n), T(T), h(T/pow(2., n)), G(0,sqrt(h)) {}

// genere une nouvelle liste présentant notre Mouvement Brownien
// procède par incrementation
Brownian::result_type Brownian::operator()() {
    
    value.clear();
    state val_k(0,0);
    value.push_back(val_k);
    
    do {
        val_k.first += h;
        val_k.second += G();
        value.push_back(val_k);
    } while (val_k.first < T);
    
    return value;

}

// cette fonction multiplie par 2 le nb de valeurs du processus et affine donc le pas de subdivision le
// rendant plus petit.
// toutes les autres caracteristiques sont changés : n, h, G, value mais pas T
Brownian::result_type Brownian::affine() {
    
    n++;
    h *= 0.5;
    G = Gaussian(0, sqrt(0.5*h));
    
    iter precedent = value.begin(),
    current = ++value.begin();
    
    while (current != value.end()) {
        value.insert(current, state(0.5*((*current).first+(*precedent).first),
                                    0.5*((*current).second+(*precedent).second)+G()));
        precedent = current;
        current++;
    }
    return value;
}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___BROWNIAN___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___PONT_BROWNIAN___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// constructeur va construire le processus de taille 2^n + 1, et de pas de temps T/2^n
Brownian_bridge::Brownian_bridge(int n, double T, double B_T) : Brownian(n, T), B_T(B_T) {}

// crée notre pont brownien en fixant la condition initiale et condition finale puis utilise
// la fonction affine hérité du brownien pour affiner le pas de discretisation h.
// notre nombre de valeurs est 2^n+1
Brownian_bridge::result_type Brownian_bridge::operator()() {
    
    value.clear();
    value.push_back(state(0,0));
    value.push_back(state(T, B_T));
    
    int n_tmp = n;
    n = 0;
    h = T;
    
    for (int j = 0; j < n_tmp; j++) {
        affine();
    }
    
    n = n_tmp;
    return value;
}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___PONT_BROWNIAN___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___BROWNIAN_BIAISE___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Brownian_biased::Brownian_biased(int n, Var_alea<double> & X, double T) : Brownian(n, T), X(X) {}

Brownian_biased::result_type Brownian_biased::operator()() {
    
    value.clear();
    value.push_back(state(0,0));
    value.push_back(state(T, X()));
    
    int n_tmp = n;
    n = 0;
    h = T;
    
    for (int j = 0; j < n_tmp; j++) {
        affine();
    }
    
    n = n_tmp;
    return value;
}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___BROWNIAN_BIAISE___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___BLACK & SCHOLES___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Black_scholes::Black_scholes(int n, double x0, double r, double s, double T) : Brownian(n, T), bs(x0, r, s) {}

Black_scholes::result_type Black_scholes::operator()() {
    
    Brownian::operator()();
    std::transform(value.begin(), value.end(), value.begin(), bs);
    return value;
    
}

Black_scholes::Fun_bs::Fun_bs(double x0, double r, double s): x0(x0), s(s), mu(r-0.5*s*s) {}

Black_scholes::state Black_scholes::Fun_bs::operator()(const state & x) {
    
    return state(x.first, x0*exp(mu*x.first + s*x.second));
    
}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___BLACK & SCHOLES___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___ORNSTEIN & UHLENBECK___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Ornstein_uhlenbeck::Ornstein_uhlenbeck(int n, double x0, double lambda, double mu, double s, double T) : T(T), h(pow(2.,-n)), x0(x0), ret(exp(-lambda*h)), G(mu*(1-ret), s*s*(1-exp(-2*lambda*h))/(2*lambda)) {}

Ornstein_uhlenbeck::result_type Ornstein_uhlenbeck::operator()() {
    
    value.clear();
    state val_k(0, x0);
    value.push_back(val_k);
    
    do {
        val_k.first += h;
        val_k.second = val_k.second*ret + G();
        value.push_back(val_k);
    } while (val_k.first < T);
    
    return value;
}
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___ORNSTEIN & UHLENBECK___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___POISSON___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Ppoisson::Ppoisson(double lambda, double T) : E(lambda), T(T) {}

// On simule l'instant du premier saut, et on regarde si cela ne dépasse pas T, puis un deuxieme unstant plus le premier, puis un troisieme... et ainsi de suite jusqu'a ce que notre somme des instants de sauts ne dépasse pas T.
Ppoisson::result_type Ppoisson::operator()() {
    
    value.clear();
    state val_k(0, 0);
    value.push_back(val_k);
    val_k.first += E();
    
    while (val_k.first < T) {
        val_k.second += 1;
        value.push_back(val_k);
        val_k.first += E();
    };
    
    val_k.first = T;
    value.push_back(val_k);
    return value;

}

double Ppoisson::getT(){

    return T;
    
}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___POISSON___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___POISSON COMPOSE___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Ppoisson_compose::Ppoisson_compose(double lambda, Var_alea<double> & X, double T): N(lambda, T), X(X) {}

Ppoisson_compose::result_type Ppoisson_compose::operator()() {
    
    value.clear();
    state val_k = state(0,0);
    value.push_back(val_k);
    N();
    // probleme au niveau de définitions des lists et itherators dans la STL
    // faire attention et ne pas passer en direct
    std::list<std::pair<double,unsigned>> Nvalues = N.current();
    
    for (std::list<std::pair<double,unsigned>>::iterator i = Nvalues.begin(); i!=Nvalues.end(); ++i) {
    
        val_k.first = (*i).first;
        val_k.second += X();
        value.push_back(val_k);
    
    }
    
    value.push_back(state(N.getT(), val_k.second));
    return value;
    
}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___POISSON COMPOSE___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
