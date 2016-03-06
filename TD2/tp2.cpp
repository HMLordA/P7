#include "random-singleton.h"
#include "random-helper.h"
#include <cstdio>
#include <iostream>
#include <iomanip>

static const double T = 10.0;
static const double S0 = 100.0;
static const double sigma = 0.2;
static const double r = 0.05;
static const double L = 110.0;

double Payoff(double G)
{

    //const double L = 60.0;
    
    double ST = S0*exp((r-sigma*sigma/2)*T+sigma*sqrt(T)*G);
    return (ST > L)?ST-L:0.0;
}

//Gradient calculation
double nablacomputation(double theta, double* G, int M)
{
    double res = 0.0;
    for (int i = 0; i < M; i++)
    {
        double Gi = G[i];
        double payoff = Payoff(Gi);
        res += (theta - Gi)*payoff*payoff*exp(-theta*Gi+theta*theta/2);
    }

    return res/M;
}

//Hessian calculation
double nabla2computation(double theta, double* G, int M)
{
    double res = 0.0;
    for (int i = 0; i < M; i++)
    {
        double Gi = G[i];
        double payoff = Payoff(Gi);
        res += (1 + (theta - Gi)*(theta - Gi))*payoff*payoff*exp(-theta*Gi+theta*theta/2);
    }

    return res/M;
}

double NewtonRaphson1d(double *G, double epsilon, int M)
{
    // G est un vecteur contenant l'échantillon
    // epsilon est la précision, M le nombre de simulation
    double theta, nablav, nabla2v;
    theta = 0; //par exemple
    nablav = nablacomputation(theta, G, M); //calcul du gradient de v
    nabla2v = nabla2computation(theta, G, M); //calcul de la hessienne
    while(abs(nablav) > epsilon)
    {
        theta -= nablav/nabla2v;
        nablav = nablacomputation(theta, G, M);
       nabla2v = nabla2computation(theta, G, M); 
    }
    return theta;
}

int main(int argc, const char** argv)
{
    Random::Randomize(0);
    cout<<setprecision(5);
    const double M = atof(argv[1]);

    //1 Prix standard de l'option
    {
        double res_sum = 0.0;
        double res_sum_sq = 0.0;
        for (int i = 0; i < M; i++)
        {
            //double payoff = Payoff(NormaleBoxMuller(0,1),S0,r,sigma,T,L);
            double payoff = Payoff(Random::Gaussian());
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
    //2 Estimation of the optimal theta
    //Generation of the gaussian vector
    double* v_gaussian = new double[(int)M];
    double epsilon = 0.000001;
    for (int i=0;i<M;i++)
    {
        v_gaussian[i] = Random::Gaussian();
    }
    double optimal_theta = NewtonRaphson1d(v_gaussian,epsilon,M);

    //3 Estimator calculation + variance + confidence interval
    {
        double res_sum = 0.0;
        double res_sum_sq = 0.0;
        for (int i = 0; i < M; i++)
        {
            //double payoff = Payoff(NormaleBoxMuller(0,1),S0,r,sigma,T,L);
            double Gi = Random::Gaussian();
            double payoff = Payoff(Gi+optimal_theta)*exp(-optimal_theta*Gi-optimal_theta*optimal_theta/2);
            res_sum += payoff;
            res_sum_sq += payoff * payoff;
        }
        double mean = exp(-r*T)*res_sum/M;
        double var = 1.0/(M-1)*(exp(-2*r*T)*res_sum_sq - M*mean*mean);

        cout<<endl<<"Optimal theta : "<<optimal_theta<<endl;
        cout<<"MC Digital on "<<M<< " scenarii"<<endl;
        cout<<"Optimized Mean:"<<mean<<endl;
        cout<<"Optimized Variance:"<<var<<endl;
        cout<<"95perc-confidence interval for the price:";
        cout<<"["<<mean-1.96*sqrt(var/M)<<";"<<mean+1.96*sqrt(var/M)<<"]"<<endl;
    }

    //3.2 Best of - Prix standard
    /*{
        double res_sum = 0.0;
        double res_sum_sq = 0.0;
        for (int i = 0; i < M; i++)
        {
            //double payoff = Payoff(NormaleBoxMuller(0,1),S0,r,sigma,T,L);
            double payoff = PayoffBestOf(Random::Gaussian(),Random::Gaussian());
            res_sum += payoff;
            res_sum_sq += payoff * payoff;
        }
        double mean = res_sum/M;
        double var = 1.0/(M-1)*(res_sum_sq - M*mean*mean);

        cout<<"MC Digital on "<<M<< " scenarii"<<endl;
        cout<<"Standard Mean:"<<mean<<endl;
        cout<<"Standard Variance:"<<var<<endl;
        cout<<"95perc-confidence interval for the price:";
        cout<<"["<<mean-1.96*sqrt(var/M)<<";"<<mean+1.96*sqrt(var/M)<<"]"<<endl;
    }*/
    //3.2 Estimator calculation + variance + confidence interval
    /*{
        double res_sum = 0.0;
        double res_sum_sq = 0.0;
        for (int i = 0; i < M; i++)
        {
            //double payoff = Payoff(NormaleBoxMuller(0,1),S0,r,sigma,T,L);
            double Gi = Random::Gaussian();
            double payoff = Payoff(Gi+optimal_theta)*exp(-optimal_theta*Gi-optimal_theta*optimal_theta/2);
            res_sum += payoff;
            res_sum_sq += payoff * payoff;
        }
        double mean = res_sum/M;
        double var = 1.0/(M-1)*(res_sum_sq - M*mean*mean);

        cout<<endl<<"Optimal theta : "<<optimal_theta<<endl;
        cout<<"MC Digital on "<<M<< " scenarii"<<endl;
        cout<<"Optimized Mean:"<<mean<<endl;
        cout<<"Optimized Variance:"<<var<<endl;
        cout<<"95perc-confidence interval for the price:";
        cout<<"["<<mean-1.96*sqrt(var/M)<<";"<<mean+1.96*sqrt(var/M)<<"]"<<endl;
    }*/
return 0;
}
