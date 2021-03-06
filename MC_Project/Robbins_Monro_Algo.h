//
//  Robbins_Monro_Algo.h
//  MC_Project
//
//  Created by Nazar KOSTYUCHYK - JC DIETRICH on 01/03/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK - JC DIETRICH. All rights reserved.
//

#ifndef Robbins_Monro_Algo_h
#define Robbins_Monro_Algo_h

#include <stdio.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <vector>
#include <list>
#include "Tools.h"
#include "Theta.h"
#include <iostream>
#include <fstream>

#include "Variables.h"

// Polynome Legendre
#include <boost/math/special_functions/legendre.hpp>

typedef boost::numeric::ublas::vector<double> Vector;


using namespace std;

//RM Algo for European products
template<class T, class S, double (T::*F_Payoff)(const Vector &) const, double (T::*F_Tilda_Control)(const Vector &) const>
Vector Robbins_Monro_Algo(int M, double alpha, double gamma0, Vector theta, double c, const T& Obj, S& G){

    double S1=0.0;
    double S2=0.0;
    Vector g;
    double Sreal1=0.0;
    double Sreal2=0.0;
	double Sreal1Drift=0.0;
    double Sreal2Drift=0.0;
	double thetaDrift = 1.6;

    int counter = 0;
    for (int n=0; n<M; ++n){
        
        g = G();
        if (n == counter)
		{
			cout <<"Theta at "<<n<<" : ";
			for (unsigned int ll=0;ll<theta.size();ll++)
				cout<<theta[ll]<<",";
			cout<<endl;
			if (n < 100)
				counter+=1;
			else
				counter += M/100;
		}

		//Update of the multidim theta
		Vector d_theta = ((gamma0/pow(n+1,alpha))) * (2*theta - g) * (pow((Obj.*F_Payoff)(g-theta),2)/(1+pow((Obj.*F_Tilda_Control)(-theta),2 * c)));
        theta = theta - d_theta;

        S1 += (Obj.*F_Payoff)(g+theta) * exp(inner_prod(-theta,g) - 0.5 * inner_prod(theta,theta));
        S2 += pow((Obj.*F_Payoff)(g+theta) * exp(inner_prod(-theta,g) - 0.5 * inner_prod(theta,theta)),2);

        Sreal1 += (Obj.*F_Payoff)(g);
        Sreal2 += pow((Obj.*F_Payoff)(g),2);

		Vector v_thetaDrift = Vector(g.size(),thetaDrift);
		Sreal1Drift += (Obj.*F_Payoff)(g+v_thetaDrift)*exp(inner_prod(-v_thetaDrift,g) - 0.5 * inner_prod(v_thetaDrift,v_thetaDrift));
        Sreal2Drift += pow((Obj.*F_Payoff)(g+v_thetaDrift)*exp(inner_prod(-v_thetaDrift,g) - 0.5 * inner_prod(v_thetaDrift,v_thetaDrift)),2);
        
    }
    
    cout <<"Le theta optimal : ";
	for (unsigned int ll=0;ll<theta.size();ll++)
		cout<<theta[ll]<<",";
	cout<<endl;
    
	double meanSimple = exp(-Obj.r*Obj.T) * Sreal1/M;
	double meanDrift = exp(-Obj.r*Obj.T) * Sreal1Drift/M;
	double meanChanged = exp(-Obj.r*Obj.T) * S1/M;
	double meanSqSimple = exp(-2*Obj.r*Obj.T)*Sreal2/M;
	double meanSqDrift = exp(-2*Obj.r*Obj.T)*Sreal2Drift/M;
	double varSimple = meanSqSimple-exp(-2*Obj.r*Obj.T)*pow(Sreal1/M,2);
	double varDrift = meanSqDrift-exp(-2*Obj.r*Obj.T)*pow(Sreal1Drift/M,2);
	double meanSqChanged = exp(-2*Obj.r*Obj.T)*S2/M;
	double varChanged = meanSqChanged-exp(-2*Obj.r*Obj.T)*pow(S1/M,2);  
    

	cout <<"L'esperance simple : "<< meanSimple << endl;
	cout <<"L'esperance simple avec drift "<<thetaDrift<<" : "<< meanDrift << endl;
    cout <<"L'esperence avec changement : "<< meanChanged << endl;

	cout <<"L'esperance du carre simple : "<< meanSqSimple << endl;
	cout <<"L'esperance du carre avec drift : "<< meanSqDrift << endl;
    cout <<"L'esperance du carre avec changement : "<< meanSqChanged << endl;
    
    cout <<"La variance simple : "<< varSimple << endl;
	cout <<"La variance avec drift : "<< varDrift << endl;
    cout <<"La variance avec changement : "<< varChanged << endl;

	cout<<"Intervalle de confiance pour le prix simple:";
    cout<<"["<<meanSimple-1.96*sqrt(varSimple/M)<<";"<<meanSimple+1.96*sqrt(varSimple/M)<<"]"<<endl;
	cout<<"Intervalle de confiance pour le prix avec drift:";
    cout<<"["<<meanDrift-1.96*sqrt(varDrift/M)<<";"<<meanDrift+1.96*sqrt(varDrift/M)<<"]"<<endl;
	cout<<"Intervalle de confiance pour le prix avec changement:";
    cout<<"["<<meanChanged-1.96*sqrt(varChanged/M)<<";"<<meanChanged+1.96*sqrt(varChanged/M)<<"]"<<endl;
    
    return theta;

}

//RM Algo for products priced by EDS
template<class T, class S, class SS, class U, double (T::*F_Payoff)(const std::list<std::pair<double,double>> &) const>
void Robbins_Monro_SDE_Algo(int M, double alpha, double gamma0, Theta* mytheta, double c, const T& Obj, S& EDS1, SS& EDS2, U& G){
    
    double S1=0.0;
    double S2=0.0;
    double Sreal1=0.0;
    double Sreal2=0.0;

	double JJ = 0.0;
	U mynorm = U();
	Theta* plusTheta = mytheta->copy();
	//Theta_Legendre plusTheta(mytheta);

	//LEGENDRE : uncomment for Legendre run
	/*vector<LegendreCarre> LC;
	for (unsigned int j=0;j<mytheta->getTh().size();j++)
	{     
        LC.push_back(LegendreCarre(j));    
    }
	vector<double> th = mytheta->getTh();
	Theta_Legendre_squared* th_sq= new Theta_Legendre_squared(th);
	*/

	//HAAR
	vector<HaarCarre> LC;
	int n = int(log(double(mytheta->getTh().size())+1.0)/log(2.0))-1;
	LC.push_back(HaarCarre(-1,-1));
    for (int current_n=0;current_n<=n;current_n++)
	{
		double j = pow(2.0,current_n)-1; 
		for(int current_j=0; current_j<=j; current_j++){       
        //LC.push_back(LegendreCarre(j));
        LC.push_back(HaarCarre(current_n,current_j));
		}     
    }
 	vector<double> th = mytheta->getTh();
	Theta_Haar_squared* th_sq = new Theta_Haar_squared(th); 
	
	
	int counter = 0;

    for (int n=0; n<M; ++n){
        
        //EDS1.setNewTheta(mytheta);
        
		//simulation of the diffusions used
        EDS1();  
        EDS2();
        
		//print of the first steps
		if (n == counter)
		{
			cout <<"Theta at "<<n<<" : " ;
			for (unsigned int ll=0;ll<mytheta->getTh().size();ll++)
				cout<<mytheta->getTheta_i(ll)<<",";
			cout<<endl;
			/*if (n < 100)
				counter+=1;
			else*/
				counter += M/100;
		}
		
		//double sum_theta = integral1P<double>(0.0, 1, 0.01, th_sq);
		vector<double> thet = mytheta->getTh();
		double sum_theta = 0.0;
        for (auto th =thet.begin();th!=thet.end();th++){
            sum_theta += pow(*th,2.0);
        }

        // Thetas update
        for(unsigned int i=0; i<mytheta->getTh().size(); i++){
            
            JJ = integralSTO( EDS1.current(), LC[i]);
            
            mytheta->setTheta_i( i, mytheta->getTh()[i] - (gamma0/(pow(n+1,alpha)+0*0.0001)) * 1/(1+sum_theta)/*exp(sum_theta)*/*( pow((Obj.*F_Payoff)(EDS1.current_BS_Drift()),2)*( 2*mytheta->getTh()[i] - JJ ) ) );
            th_sq->setTheta_i(i,mytheta->getTheta_i(i));
            plusTheta->setTheta_i(i,-mytheta->getTheta_i(i));
            
        }
		
        //Rediffusion with the same brownian but updated thetas
        EDS1.setNewTheta(plusTheta);
        EDS1.rediffusion();
        
        double gg = integralSTO(EDS1.current(), *mytheta);     
		sum_theta = integral1P<double>(0.0, 1, 0.01, *th_sq);

        /*for (auto th =thet.begin();th!=thet.end();th++){
            sum_theta += pow(*th,2.0);
        }*/

		S1 += (Obj.*F_Payoff)(EDS1.current_BS_Drift()) * exp (-(gg + sum_theta/2));
        S2 += pow((Obj.*F_Payoff)(EDS1.current_BS_Drift()) * exp (-(gg + sum_theta/2)),2);
        
        Sreal1 += (Obj.*F_Payoff)(EDS2.current());
        Sreal2 += pow((Obj.*F_Payoff)(EDS2.current()),2);
        
    }

	double meanSimple = exp(-Obj.r*Obj.T) * Sreal1/M;
	double meanChanged = exp(-Obj.r*Obj.T) * S1/M;
	double meanSqSimple = exp(-2*Obj.r*Obj.T)*Sreal2/M;
	double varSimple = meanSqSimple-exp(-2*Obj.r*Obj.T)*pow(Sreal1/M,2);
	double meanSqChanged = exp(-2*Obj.r*Obj.T)*S2/M;
	double varChanged = meanSqChanged-exp(-2*Obj.r*Obj.T)*pow(S1/M,2);  
	
	//Uncomment to write the thetas Haar in a file 
	/*cout<<"Theta graphe:"<<endl;
	ofstream myfile;
	myfile.open ("E:\\Documents\\Work\\M2MO\\MC\\Project\\MC_Project\\Resultats\\Haar_res\\result.csv");
	for(double z=0.01;z<=1.0;z+=0.01)
	{
		cout<<"\t"<<z<<"\t"<<mytheta->value(z)<<endl;
	}
	for(double z=0.0;z<=1.01;)
	{
		myfile<<z<<";"<<mytheta->value(z)<<endl;
		z+=0.05;
	}
	myfile.close();
	*/
	cout <<"L'esperance simple : "<< meanSimple << endl;
    cout <<"L'esperence avec changement : "<< meanChanged << endl;
    cout <<"L'esperance du carre simple : "<< meanSqSimple << endl;
    cout <<"L'esperance du carre avec changement : "<< meanSqChanged << endl;
    
    cout <<"La variance simple : "<< varSimple << endl;
    cout <<"La variance avec changement : "<< varChanged << endl;
	cout<<"Intervalle de confiance pour le prix simple:";
    cout<<"["<<meanSimple-1.96*sqrt(varSimple/M)<<";"<<meanSimple+1.96*sqrt(varSimple/M)<<"]"<<endl;
	cout<<"Intervalle de confiance pour le prix avec changement:";
    cout<<"["<<meanChanged-1.96*sqrt(varChanged/M)<<";"<<meanChanged+1.96*sqrt(varChanged/M)<<"]"<<endl;
	
	delete plusTheta;
	delete th_sq;
    
}

#endif /* Robbins_Monro_Algo_h */
