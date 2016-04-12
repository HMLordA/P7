#include "random-singleton.h"
#include "random-helper.h"
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <numeric>
#include <functional>

#include <iterator>

using namespace std;

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
//#include <boost/numeric/ublas/io.hpp>

typedef boost::numeric::ublas::vector<double> Vector;
typedef boost::numeric::ublas::matrix<double> Matrix;
typedef boost::numeric::ublas::identity_matrix<double> Identity_Matrix;

#include <ctime>

static const double T = 1.0; 
static const double S01 = 100.0; static const double S02 = 95.0;
static const double sigma1 = 0.3; static const double sigma2 = 0.4;
static const double r = 0.04;

//static const double L = 115.0; //OTM
//static const double L = 100.0; //ATM
static const double L = 85.0; //ITM

static const double rho = 0.0; //0.3;
//const double L = 60.0;

static const int DIM = 2;
static const double epsilon = 0.000001;
static const double M = 1000000;
static const double N = 10000;


double Payoff(const Vector& G)
{
	double W1 = sqrt(T)*G[0];


	//Assets prices at maturity
	double ST1 = S01*exp((r-sigma1*sigma1/2)*T+sigma1*W1);

	double W3 = sqrt(T)*G[1];
	//Generation of the correlated brownian
	double W2 = rho*W1+sqrt(1-rho*rho)*W3;
	double ST2 = S02*exp((r-sigma2*sigma2/2)*T+sigma2*W2);
	
	return max(max(ST1,ST2)-L,0.0);
	//return max(ST1-L,0.0);
}

//Gradient calculation
Vector nablacomputation(const Vector& theta, const Matrix& G, int M)
{
	Vector res = Vector(G.size2(),0.0);
	for (int i = 0; i < M; i++)
	{
		const Vector& Gi = row(G,i);
		double payoff = Payoff(Gi);
		double theta_Gi = inner_prod(theta,Gi);
		double sq_norm_theta = inner_prod(theta,theta);

		//vector<double> interm_res = vector<double>(2,0.0);
		double shift = payoff*payoff*exp(-theta_Gi+sq_norm_theta/2);
		res += (theta-Gi)*shift;
	}
	res /= M;
	return res;
}

//Hessian calculation
Matrix nabla2computation(const Vector& theta, const Matrix& G, int M)
{
	Matrix res = Matrix(G.size2(),G.size2(),0.0);
	Matrix id = Identity_Matrix(G.size2());
	for (int i = 0; i < M; i++)
	{
		const Vector& Gi = row(G,i);
		double payoff = Payoff(Gi);
		double theta_Gi = inner_prod(theta,Gi);
		double sq_norm_theta = inner_prod(theta,theta);

		//vector<double> interm_res = vector<double>(2,0.0);
		double shift = payoff*payoff*exp(-theta_Gi+sq_norm_theta/2);
		Vector theta_m_Gi = theta - Gi;
		Matrix mat_theta_m_Gi = outer_prod(theta_m_Gi,theta_m_Gi);
		res += (mat_theta_m_Gi + id)*shift;
		/*res[0][0] = (1 + theta_m_Gi[0]*theta_m_Gi[0])*shift;
		res[0][1] = (theta_m_Gi[0]*theta_m_Gi[1])*shift;
		res[1][0] = (theta_m_Gi[0]*theta_m_Gi[0])*shift;
		res[1][1] = (1 + theta_m_Gi[1]*theta_m_Gi[1])*shift;*/
		/*double res = 0.0;
		for (int i = 0; i < M; i++)
		{
		double Gi1 = G[0][i];
		double Gi2 = G[1][i];
		double payoff = Payoff(Gi1,Gi2);
		res += (1 + (theta - Gi)*(theta - Gi))*payoff*payoff*exp(-theta*Gi+theta*theta/2);
		}*/
	}
	res /= M;
	return res;
}

Vector NewtonRaphson1d(const Matrix& G, double epsilon, int M)
{
	// G est un vecteur contenant l'échantillon
	// epsilon est la précision, M le nombre de simulation
	Vector theta=Vector(G.size2(),0.0);
	//theta = 0; //par exemple
	Vector nablav = nablacomputation(theta, G, M); //calcul du gradient de v
	Matrix nabla2v = nabla2computation(theta, G, M); //calcul de la hessienne

	double norm = sqrt(inner_prod(nablav,nablav));
	while(norm > epsilon)
	{
		cout<<"Norm:"<<norm<<endl;
		cout<<endl<<"Theta : ";
		for (int k=0;k<G.size2();k++)
			cout<<theta[k]<<",";
		cout<<endl;

		Matrix inv_nabla2v = Matrix(G.size2(),G.size2(),0.0);
		InvertMatrix(nabla2v,inv_nabla2v); 
		theta -= prod(inv_nabla2v,nablav);

		/*double determ = nabla2v[0][0]*nabla2v[1][1]-nabla2v[1][0]*nabla2v[0][1];
		theta[0] -= 1/determ*(nabla2v[1][1]*nablav[0]-nabla2v[0][1]*nablav[1]);
		theta[1] -= 1/determ*(-nabla2v[1][0]*nablav[0]+nabla2v[0][0]*nablav[1]);*/
		nablav = nablacomputation(theta, G, M);
		nabla2v = nabla2computation(theta, G, M); 
		norm = sqrt(inner_prod(nablav,nablav));
	}
	cout<<"Norm:"<<norm<<endl;
	cout<<endl<<"Theta : ";
	for (int k=0;k<G.size2();k++)
		cout<<theta[k]<<",";
	cout<<endl;
	return theta;
}

int main(int argc, const char** argv)
{
	Random::Randomize(time(0));
	cout<<setprecision(5);

	//3.2 Prix standard de l'option BestOf
	{
		double res_sum = 0.0;
		double res_sum_sq = 0.0;
		for (int i = 0; i < M; i++)
		{
			//double payoff = Payoff(NormaleBoxMuller(0,1),S0,r,sigma,T,L);
			Vector G (DIM,0.0);
			for (int k=0;k<DIM;++k)
				G(k)=Random::Gaussian();
			double payoff = Payoff(G);
			res_sum += payoff;
			res_sum_sq += payoff * payoff;
		}
		double mean = exp(-r*T)*res_sum/M;
		double var = 1.0/(M-1)*(exp(-2*r*T)*res_sum_sq - M*mean*mean);

		cout<<"MC Digital on "<<M<< " scenarii"<<endl;
		cout<<"Standard Mean:"<<mean<<endl;
		cout<<"Standard Variance:"<<var<<endl;
		cout<<"95perc-confidence interval for the price:";
		cout<<"["<<mean-1.96*sqrt(var/M)<<";"<<mean+1.96*sqrt(var/M)<<"]"<<endl;
	}
	//Estimation of the optimal theta
	//Generation of the gaussian vector
	Matrix v_gaussian = Matrix((int)N,DIM, 0.0);

	for (int i=0;i<N;i++)
	{
		for (int j=0;j<DIM;j++)
		{
			v_gaussian(i,j) = Random::Gaussian();
		}
	}
	Vector optimal_theta = NewtonRaphson1d(v_gaussian,epsilon,N);
	
	//Vector optimal_theta = Vector(DIM,1.20629);
	cout<<endl<<"Optimal theta : ";
	for (int k=0;k<DIM;k++)
		cout<<optimal_theta[k]<<",";
	cout<<endl;


	//3 Estimator calculation + variance + confidence interval
	{
		double res_sum = 0.0;
		double res_sum_sq = 0.0;
		for (int i = 0; i < M; i++)
		{
			//double payoff = Payoff(NormaleBoxMuller(0,1),S0,r,sigma,T,L);
			Vector Gi = Vector(DIM,0.0);
			for (int k=0;k<DIM;k++)
				Gi[k] = Random::Gaussian();
			//Vector Gi_p_theta = Vector(2,0.0);
			Vector Gi_p_theta = Gi+optimal_theta;
			double theta_Gi = inner_prod(optimal_theta,Gi);
			double sq_norm_theta = inner_prod(optimal_theta,optimal_theta);
			double payoff = Payoff(Gi_p_theta)*exp(-theta_Gi-sq_norm_theta/2);
			res_sum += payoff;
			res_sum_sq += payoff * payoff;
		}
		double mean = exp(-r*T)*res_sum/M;
		double var = 1.0/(M-1)*(exp(-2*r*T)*res_sum_sq - M*mean*mean);

		cout<<"MC Digital on "<<M<< " scenarii"<<endl;
		cout<<"Optimized Mean:"<<mean<<endl;
		cout<<"Optimized Variance:"<<var<<endl;
		cout<<"95perc-confidence interval for the price:";
		cout<<"["<<mean-1.96*sqrt(var/M)<<";"<<mean+1.96*sqrt(var/M)<<"]"<<endl;
	}
}
