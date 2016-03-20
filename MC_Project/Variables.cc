//
//  variables.cpp
//  EK_MC
//
//  Created by Nazar KOSTYUCHYK on 18/02/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK. All rights reserved.
//

#include "Variables.h"
#include <ctime>
#include <iostream>

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___UNIFORM___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Uniform::Uniform(double a, double b) : a(a), size(b-a) {
    
	unsigned int seed = static_cast<unsigned int>(std::time(0));
    gen.seed(seed);
	std::cout <<"Seed : "<<seed<<std::endl;

}

// Surcharge de l'opérateur simulation simulation d'une variable uniforme
double Uniform::operator()() {
    //gen.seed(static_cast<unsigned int>(std::time(0)));
    uniform_dist u(a,a+size);
    value = u(gen);
    return value;
}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___UNIFORM___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___GAUSSIAN___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Gaussian::Gaussian(double mu, double sigma) : mu(mu), sigma(sigma) {
	
	unsigned int seed = static_cast<unsigned int>(std::time(0));
    gen.seed(seed);
    //gen.seed(0);
	std::cout <<"Seed : "<<seed<<std::endl;

}

// Surcharge de l'opérateur simulation simulation d'une variable gaussienne
double Gaussian::operator()() {
    
    //gen.seed(static_cast<unsigned int>(std::time(0)));
    normal_dist n(mu,sigma);
    value = n(gen);
    return value;
    
}

void Gaussian::setSigma(double sig){
    
    sigma=sig;
    
}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___GAUSSIAN___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___GAUSSIAN_VECTOR___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// Variable aléatoire Gausienne

// Constructeur qui peut jouer le role du constructeur par default
Gaussian_Vector::Gaussian_Vector(Vector & mu, Matrix & sigma, int size): mu(mu), sigma(sigma),size(size),G(0.0,1.0){
    
    //G = Gaussian(0.0,1.0);
    sigmaSqrt = sigma;
    cholesky_decompose<Matrix,Matrix>(sigma,sigmaSqrt);
    
}

// Surcharge de l'opérateur: Simulation d'une variable gaussienne
Vector Gaussian_Vector::operator()(){
    
    value.clear();
    
    Vector v(size);
    
    for (int i=0; i<size; i++){
        v(i)=G();
    }
    
    value = prod(sigmaSqrt, v);
    
    return value;
    
}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___GAUSSIAN_VECTOR___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
/* OLD VERSION
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___GAUSSIAN_VECTOR___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// Variable aléatoire Gausienne

    // Constructeur qui peut jouer le role du constructeur par default
Gaussian_Vector::Gaussian_Vector(gsl_vector * mu, gsl_matrix * sigma, int size): mu(mu), sigma(sigma),size(size){

    sigmaSqrt = gsl_matrix_calloc(size, size);
    
    gsl_matrix_memcpy(sigmaSqrt, sigma); // This function copies the elements of the matrix m into the matrix A. The two matrices must have the same size
    ///////ERROR POSSIBLE
    gsl_linalg_cholesky_decomp(sigmaSqrt); // Effectue la decomposition de cholesky
    
    // set to zero the upper diagonal of A
    for (int i =0; i<size; i++){
        for (int j=i+1; j<size; j++)
            gsl_matrix_set(sigmaSqrt, i, j, 0);
    }
    
}

    // Surcharge de l'opérateur: Simulation d'une variable gaussienne
std::vector<double> Gaussian_Vector::operator()(){

    value.clear();
    Gaussian G(0.0,1.0);
    double g;
    gsl_vector * F = gsl_vector_calloc(size);
    gsl_vector * x = gsl_vector_alloc(size);
    gsl_vector_memcpy(x,mu);
    
    for (int i=0; i<size; i++){
        g = G();
        gsl_vector_set(F, i, g);
    }
    
    // int gsl_blas_strmv (CBLAS UPLO t Uplo, CBLAS TRANSPOSE t TransA, CBLAS DIAG t Diag, const gsl matrix float * A, gsl vector float * x)
    // gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A, G, 0, X); // produit matriciel Compute C = A B on ecrit gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, &A.matrix, &B.matrix, 0.0, &C.matrix);
    
    gsl_blas_dgemv(CblasNoTrans, 1.0, sigmaSqrt, F, 1.0, x);
   
    for (int i=0; i<size; i++){
        
        value.push_back(gsl_vector_get(x, i));
        
    }
    
    gsl_vector_free(F);
    gsl_vector_free(x);
    
    return value;

}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___GAUSSIAN_VECTOR___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
*/
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___EXPO___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Expo::Expo(double lambda) : lambda(lambda) {
    
    gen.seed(static_cast<unsigned int>(std::time(0)));
    
}

// Surcharge de l'opérateur simulation simulation d'une variable exponentielle
double Expo::operator()() {
    
    //gen.seed(static_cast<unsigned int>(std::time(0)));
    exponentiel_dist ex(lambda);
    value = ex(gen);
    return value;
    
}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___EXPO___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___BERNOUILLI___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Bernoulli::Bernoulli(double p) : p(p) {
    
    gen.seed(static_cast<unsigned int>(std::time(0)));
    
}

// Surcharge de l'opérateur simulation simulation d'une variable exponentielle
unsigned Bernoulli::operator()() {
    
    //gen.seed(static_cast<unsigned int>(std::time(0)));
    bernoulli_dist bern(p);
    value = bern(gen);
    return value;
    
}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@___BERNOUILLI___@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@