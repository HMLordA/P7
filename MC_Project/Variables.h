//
//  Variables.h
//  EK_MC
//
//  Created by Nazar KOSTYUCHYK on 18/02/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK. All rights reserved.
//

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/cholesky.hpp>
#include <boost/numeric/ublas/io.hpp>

/*
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
*/

typedef boost::mt19937 generator;
typedef boost::normal_distribution<> normal_dist;
typedef boost::random::uniform_real_distribution<> uniform_dist;
typedef boost::random::exponential_distribution<> exponentiel_dist;
typedef boost::random::bernoulli_distribution<> bernoulli_dist;
typedef boost::variate_generator< generator&,normal_dist > normal_rv;

typedef boost::numeric::ublas::matrix<double> Matrix;
typedef boost::numeric::ublas::identity_matrix<double> Identity_Matrix;
typedef boost::numeric::ublas::vector<double> Vector;

#ifndef Variables_h
#define Variables_h

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___VARIABLE ALEA___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

template <class T>
class Var_alea
{
    
public:
    typedef T result_type;
    Var_alea() : value(0) {}
    Var_alea(T value) : value(value) {}
   // virtual ~Var_alea() {};
    virtual T operator()() = 0;
    T current() const { return value; }
    
protected:
    T value;
    generator gen;
};

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___VARIABLE ALEA___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___UNIFORM___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// Variable aléatoire uniforme
class Uniform : public Var_alea <double> {
    
public:
    // Constructeur qui peut jouer le role du constructeur par default
    Uniform(double a = 0.0, double b = 1.0);
    
    // Surcharge de l'opérateur: Simulation d'une variable uniforme
    double operator()();
    
private:
    double a, size;
};

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___UNIFORM___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___GAUSSIAN___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// Variable aléatoire Gausienne
class Gaussian : public Var_alea <double> {
    
public:
    // Constructeur qui peut jouer le role du constructeur par default
    Gaussian(double mu = 0.0, double sigma = 1.0);
    
    // Surcharge de l'opérateur: Simulation d'une variable gaussienne
    double operator()();
    
private:
    double mu, sigma;
};

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___GAUSSIAN___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___GAUSSIAN_VECTOR___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// Variable aléatoire Gausienne
class Gaussian_Vector : public Var_alea <Vector> {
    
public:
    // Constructeur qui peut jouer le role du constructeur par default
    Gaussian_Vector(Vector & mu, Matrix & sigma, int size = 0);
    
    // Surcharge de l'opérateur: Simulation d'une variable gaussienne
    Vector operator()();
    
private:
    Gaussian G;
    int size;
    Vector mu;
    Matrix sigma;
    Matrix sigmaSqrt;
};

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___GAUSSIAN_VECTOR___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
/* OLD Version
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___GAUSSIAN_VECTOR___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// Variable aléatoire Gausienne
class Gaussian_Vector : public Var_alea <std::vector<double>> {
    
public:
    // Constructeur qui peut jouer le role du constructeur par default
    Gaussian_Vector(gsl_vector * mu = nullptr, gsl_matrix * sigma = nullptr, int size = 0);
    
    // Surcharge de l'opérateur: Simulation d'une variable gaussienne
    std::vector<double> operator()();
    
private:
    int size;
    gsl_vector * mu;
    gsl_matrix * sigma;
    gsl_matrix * sigmaSqrt;
};

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___GAUSSIAN_VECTOR___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
*/
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___EXPO___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// Variable aléatoire Exponentielle
class Expo : public Var_alea <double> {
    
public:
    // Constructeur qui peut jouer le role du constructeur par default
    Expo(double lambda = 1.0);
    
    // Surcharge de l'opérateur: Simulation d'une variable exponentielle
    double operator()();
    
private:
    double lambda;
};

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___EXPO___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___BERNOUILLI___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// Variable aléatoire Bernoulli
class Bernoulli : public Var_alea <unsigned> {
    
public:
    // Constructeur qui peut jouer le role du constructeur par default
    Bernoulli(double p = 0.5);
    
    // Surcharge de l'opérateur: Simulation d'une variable bernoulli
    unsigned operator()();
    
private:
    double p;
};

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___BERNOUILLI___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#endif /* Variables_h */


