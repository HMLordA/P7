//
//  main.cpp
//  MC_Project
//
//  Created by Nazar KOSTYUCHYK - JC DIETRICH on 01/03/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK - JC DIETRICH. All rights reserved.
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
	int M=100000;
	int n = 5; const int NB_ASSETS = 1; const int DIM_THETA = 2;
	Payoff_Choice my_payoff = EDS_BARRIER_HAAR; // EURO_CALL // EURO_BEST_OF // EDS_CALL_HAAR // EDS_CALL_LEGENDRE // EDS_BARRIER_LEGENDRE // EDS_BARRIER_HAAR 
	
	if (my_payoff== EURO_CALL){
		//S=100.0; K=130.0; T=1.0; vol=0.3; r=0.04; alpha=0.5001; gamma0=1.0; //OTM
		//S=100.0; K=105.0; T=1.0; vol=0.3; r=0.04; alpha=0.5001; gamma0=1.0; //ATM
		S=100.0; K=70.0; T=1.0; vol=0.3; r=0.04; alpha=0.5001; gamma0=1.0; //ITM

	}else if (my_payoff== EURO_BEST_OF){
		//S=100.0; S2 = 95.0; K=130.0;  T=1.0; vol=0.3; vol2=0.4; r=0.04; alpha=0.5001; gamma0=1.0; //OTM
		//S=100.0; S2 = 95.0; K=105.0; T=1.0; vol=0.3; vol2=0.4; r=0.04; alpha=0.7501; gamma0=1.0; //ATM
		S=100.0; S2 = 95.0; K=70.0; T=1.0; vol=0.3; vol2=0.4; r=0.04; alpha=0.7501; gamma0=1.0; //ITM

	}else if (my_payoff== EDS_CALL_LEGENDRE){
		S=100.0; K=130.0; L=60.0; T=1.0; vol=0.3; r=0.04; alpha=0.5001; gamma0=0.0001; n=0;
		//S=100.0; K=115.0; L=60.0; T=1.0; vol=0.7; r=0.04; alpha=0.5001; gamma0=0.0001; n=5;
	}else if (my_payoff== EDS_CALL_HAAR){
		S=100.0; K=115.0; L=60.0; T=1.0; vol=0.7; r=0.04; alpha=0.5001; gamma0=0.0001;
	}else if (my_payoff== EDS_BARRIER_LEGENDRE){
		//S=100.0; K=105.0; L=90.0; T=1.0; vol=0.5; r=0.04; alpha=0.7501; gamma0=0.0001; n=5; //Scen 1
		S=100.0; K=120.0; L=60.0; T=1.0; vol=0.5; r=0.04; alpha=0.5001; gamma0=0.001; n=5; //Scen 2
	}else if (my_payoff== EDS_BARRIER_HAAR){
		//S=100.0; K=105.0; L=90.0; T=1.0; vol=0.5; r=0.04; alpha=0.5001; gamma0=0.0001; n=5; //Scen 1
		S=100.0; K=120.0; L=60.0; T=1.0; vol=0.5; r=0.04; alpha=0.5001; gamma0=0.001; n=5; //Scen 2
	}

	vector<double> theta;

	//PAYOFF
	Robbins_Monro_Payoff_EDS* rmbCID = 0;
	if (my_payoff== EDS_CALL_LEGENDRE||my_payoff== EDS_CALL_HAAR)
		rmbCID=new Robbins_Monro_Call_EDS(S, T, vol, r, K);
	else if (my_payoff== EDS_BARRIER_LEGENDRE||my_payoff== EDS_BARRIER_HAAR)
		rmbCID=new Robbins_Monro_CallDownIn(S, T, vol, r, L, K);
	Gaussian G(0.0,1.0);
	Black_scholes BS(n, S, r, vol, T);
	BS();


	//EDS
	if (my_payoff== EDS_CALL_LEGENDRE ||my_payoff== EDS_CALL_HAAR ||my_payoff== EDS_BARRIER_LEGENDRE ||my_payoff== EDS_BARRIER_HAAR){

		if (my_payoff== EDS_CALL_LEGENDRE ||my_payoff== EDS_BARRIER_LEGENDRE)
		{		
			for (int i=0; i<DIM_THETA; i++){
				theta.push_back(0.0);
			}
			Theta_Legendre* thetaL = new Theta_Legendre(theta);
			BS_Drift_t<Theta> BS_Drift(n, S, r, vol, thetaL, T);
			BS_Drift();
			Robbins_Monro_SDE_Algo<Robbins_Monro_Payoff_EDS, BS_Drift_t<Theta>, Black_scholes, Gaussian, &Robbins_Monro_Payoff_EDS::Payoff_Call>(M, alpha, gamma0, thetaL, c, *rmbCID, BS_Drift, BS, G);
			delete thetaL;
		}
		//HAAR
		if (my_payoff== EDS_CALL_HAAR ||my_payoff== EDS_BARRIER_HAAR)
		{
			for (int i=0; i<std::pow(2.0,DIM_THETA+1.0); i++){
				theta.push_back(0.0);
			}
			Theta_Haar* thetaL = new Theta_Haar(theta);
			// Creation of the diffusion object
			BS_Drift_t<Theta> BS_Drift(n, S, r, vol, thetaL, T);
			// First simulation of the diffusion
			BS_Drift(); 
			// Run of the RM algo
			Robbins_Monro_SDE_Algo<Robbins_Monro_Payoff_EDS, BS_Drift_t<Theta>, Black_scholes, Gaussian, &Robbins_Monro_Payoff_EDS::Payoff_Call>(M, alpha, gamma0, thetaL, c, *rmbCID, BS_Drift, BS, G);
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
		
		Vector thet = Robbins_Monro_Algo<Robbins_Monro_Payoff, Gaussian_Vector, &Robbins_Monro_Payoff::Payoff, &Robbins_Monro_Payoff::StockBS>(M, alpha, gamma0, thetaOld, cOld, *rmc, GOld);
		delete rmc;
	}
	if (rmbCID!=0)
		delete rmbCID;

	return 0;
}
