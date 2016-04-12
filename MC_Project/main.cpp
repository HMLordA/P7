//
//  main.cpp
//  MC_Project
//
//  Created by Nazar KOSTYUCHYK on 01/03/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK. All rights reserved.
//

#include <iostream>
#include "Variables.h"
#include <cmath>
#include "Robbins_Monro_Call.h"
#include "Robbins_Monro_BestOfCall.h"
#include "Robbins_Monro_CallDownIn.h"
#include "Robbins_Monro_Call_EDS.h"
#include "Robbins_Monro_Algo.h"
#include "Robbins_Monro_Algo_Normal_Distrib.h"
#include "Processus.h"

using namespace std;


double StockBS(double S0, double T, double vol, double r, double G){

	return S0 * exp((r - 0.5 * pow(vol,2)) * T + vol * sqrt(T) * G );

}

double Payoff_call (double S0, double T, double vol, double r, double G, double K){

	double S = StockBS(S0,T,vol,r,G);

	return (S-K>0)?S-K:0;

}

// Best of Calls solution Exacte
double Payoff_maxcall (double L, double G1, double G2, double S1, double S2, double T, double vol1, double vol2, double r){
	double ST1 = S1*exp( (r-pow(vol1,2)/2) *T + vol1*G1 );
	double ST2 = S2*exp( (r-pow(vol2,2)/2) *T + vol2*G2 );
	double payoff = max(max(ST1, ST2)-L, 0.0);
	return payoff;
}

enum Payoff_Choice{
	EURO_CALL,
	EURO_BEST_OF,
	EDS_CALL_LEGENDRE,
	EDS_CALL_HAAR,
	EDS_BARRIER_LEGENDRE,
	EDS_BARRIER_HAAR
};

int main(int argc, const char * argv[]) {


	//Values by default
	double S=100.0; double K=115.0; double L=60.0; double T=1.0; double vol=0.7; double r=0.04; double alpha=0.5001; double gamma0=0.0001;
	double S2 = 105.0; double vol2 = 0.3;
	//Choice of parameters depending of the payoff
	
	double c=1.0;
	int M=1000000;
	int n = 5; const int NB_ASSETS = 2; const int DIM_THETA = 4;
	Payoff_Choice my_payoff = EURO_BEST_OF; // EURO_CALL // EURO_BEST_OF // EDS_CALL // EDS_BARRIER_LEGENDRE // EDS_BARRIER_HAAR 
	
	if (my_payoff== EURO_CALL){
		//S=100.0; K=115.0; L=60.0; T=1.0; vol=0.3; r=0.04; alpha=0.5001; gamma0=0.00001; //OTM
		//S=100.0; K=100.0; L=60.0; T=1.0; vol=0.3; r=0.04; alpha=0.5001; gamma0=0.00001; //ATM
		S=100.0; K=85.0; L=60.0; T=1.0; vol=0.3; r=0.04; alpha=0.7501; gamma0=0.00001; //ITM

	}else if (my_payoff== EURO_BEST_OF){
		//S=100.0; S2 = 95.0; K=115.0; L=60.0; T=1.0; vol=0.3; vol2=0.4; r=0.04; alpha=0.7501; gamma0=0.0001; //OTM
		//S=100.0; S2 = 95.0; K=100.0; L=60.0; T=1.0; vol=0.3; vol2=0.4; r=0.04; alpha=0.7501; gamma0=0.0001; //ATM
		S=100.0; S2 = 95.0; K=85.0; L=60.0; T=1.0; vol=0.3; vol2=0.4; r=0.04; alpha=0.7501; gamma0=0.0001; //ITM

	}else if (my_payoff== EDS_CALL_LEGENDRE){
		S=100.0; K=115.0; L=60.0; T=1.0; vol=0.7; r=0.04; alpha=0.5001; gamma0=0.0001;
	}else if (my_payoff== EDS_CALL_HAAR){
		S=100.0; K=115.0; L=60.0; T=1.0; vol=0.7; r=0.04; alpha=0.5001; gamma0=0.0001;
	}else if (my_payoff== EDS_BARRIER_LEGENDRE){
		S=100.0; K=115.0; L=60.0; T=1.0; vol=0.7; r=0.04; alpha=0.5001; gamma0=0.0001;
	}else if (my_payoff== EDS_BARRIER_HAAR){
		S=100.0; K=115.0; L=60.0; T=1.0; vol=0.7; r=0.04; alpha=0.5001; gamma0=0.0001;
	}

	vector<double> theta;

	//PAYOFF
	Robbins_Monro_CallDownIn rmbCID(S, T, vol, r, L, K);
	Gaussian G(0.0,1.0);
	Black_scholes BS(n, S, r, vol, T);
	BS();


	//EDS
	if (my_payoff== EDS_CALL_LEGENDRE ||my_payoff== EDS_CALL_HAAR ||my_payoff== EDS_BARRIER_LEGENDRE ||my_payoff== EDS_BARRIER_HAAR){
		for (int i=0; i<DIM_THETA; i++){
			theta.push_back(0.0);
		}
		if (my_payoff== EDS_CALL_LEGENDRE ||my_payoff== EDS_BARRIER_LEGENDRE)
		{
			Theta_Legendre* thetaL = new Theta_Legendre(theta);
			BS_Drift_t<Theta> BS_Drift(n, S, r, vol, thetaL, T);
			BS_Drift();
			Robbins_Monro_SDE_Algo<Robbins_Monro_CallDownIn, BS_Drift_t<Theta>, Black_scholes, Gaussian, &Robbins_Monro_CallDownIn::Payoff_Call>(M, alpha, gamma0, thetaL, c, rmbCID, BS_Drift, BS, G);

			delete thetaL;
		}
		//HAAR
		if (my_payoff== EDS_CALL_HAAR ||my_payoff== EDS_BARRIER_HAAR)
		{
			for (int i=0; i<std::pow(2.0,DIM_THETA+1.0); i++){
				theta.push_back(0.0);
			}
			Theta_Haar* thetaL = new Theta_Haar(theta);
			BS_Drift_t<Theta> BS_Drift(n, S, r, vol, thetaL, T);
			BS_Drift();
			Robbins_Monro_SDE_Algo<Robbins_Monro_CallDownIn, BS_Drift_t<Theta>, Black_scholes, Gaussian, &Robbins_Monro_CallDownIn::Payoff_Call>(M, alpha, gamma0, thetaL, c, rmbCID, BS_Drift, BS, G);
			delete thetaL;
		}
	
		//EUROPEAN
	} else if (my_payoff==EURO_CALL||my_payoff==EURO_BEST_OF)
	{

		//Old call
		double gamma0 = 1.0;
		Matrix m(NB_ASSETS,NB_ASSETS,0);
		for (unsigned i = 0; i < m.size1 (); ++ i){
			m (i, i) = 1;
		}
		Vector thetaOld(NB_ASSETS,0.0);
		double cOld=1.0;

		Vector v(NB_ASSETS,0);

		Robbins_Monro_Payoff* rmc = 0;
		if (my_payoff==EURO_CALL)
			rmc = new Robbins_Monro_Call(S, K, T, vol, r);
		else if (my_payoff==EURO_BEST_OF)
			rmc = new Robbins_Monro_BestOfCall(S,S2, K, T, vol, vol2, r); //TO CHANGE
		Gaussian_Vector GOld(v,m,NB_ASSETS);


		//if (my_payoff==EURO_CALL)
		Vector thet = Robbins_Monro_Algo<Robbins_Monro_Payoff, Gaussian_Vector, &Robbins_Monro_Payoff::Payoff, &Robbins_Monro_Payoff::StockBS>(M, alpha, gamma0, thetaOld, cOld, *rmc, GOld);
		/*else if (my_payoff==EURO_BEST_OF)
			Vector thet = Robbins_Monro_Algo<Robbins_Monro_Payoff, Gaussian_Vector, &Robbins_Monro_Payoff::Payoff, &Robbins_Monro_Payoff::StockBS>(M, alpha, gamma0, thetaOld, c, *rmc, GOld);
		//Vector thet = Robbins_Monro_Algo_Normal_Distrib<Robbins_Monro_Call, Gaussian_Vector, &Robbins_Monro_Call::Payoff_Call>(M, alpha, gamma0, theta, c, rmc, G);*/
		delete rmc;
	}

	return 0;
}
