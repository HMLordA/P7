#include "random-helper.h"
#include "random-singleton.h"
#include <cmath>
#include <vector>
#include <map>
#include <ctime>
#include <exception>
#include <algorithm>
#include <cstdlib>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


using namespace std;

static const double e1 = exp(1.0);

int Bernoulli(double p)
{
    double rand = Random::Uniform();
    return (rand>p)?1:0;
}

int Binomiale(double p, int n)
{
    int res = 0;
    for (int i = 0;i<n;++i)
        res += Bernoulli(p);

    return res;
}

double Exponentielle(double lambda)
{
    return -log(Random::Uniform())/lambda;
}

double Cauchy(double x0, double lambda)
{
    return x0 + lambda * tan(M_PI*Random::Uniform()-0.5); 
}

double Normale(double mean, double stddev)
{
    double c0 = 2.515517;
    double c1 = 0.802853;
    double c2 = 0.010328;
    double d1 = 1.432788;
    double d2 = 0.189269;
    double d3 = 0.001308;
    double rand = Random::Uniform();
    int sign = -1;
    if (rand > 0.5)
    {
        sign = 1.0;
        rand = 1.0 - rand;
    }
    double t = sqrt(-2.0 * log(rand));
    double res = sign * (t-((c2*t+c1)*t+c0)/(1.0+t*(d1+t*(d2+d3*t))));
    return mean + res * stddev;
}

double NormaleBoxMuller(double mean, double stddev)
{
    double R = sqrt(Exponentielle(0.5));
    double theta = 2 * M_PI * Random::Uniform(); 
    return mean + R * cos(theta) * stddev;
}

double G_inv(double a, double z)
{
    //z is supposed between 0 and 1
    
    //Defensive programming
    if (z<=0)
        return -10000000; //return a small value
    if (z>=1)
        return 10000000; //return a big value
    
    if (z < (e1/(a+e1)))
    {
        return pow(z*(a+e1)/e1,1.0/a);
    } else {
        return -log((1-z)*(a+e1)/(a*e1));
    }
}

double q(double a, double x)
{
    if (x<=0)
        return 0.;
    if (x<1)
        return exp(-x);
    else
        return pow(x,a-1);
}

double GammaByRejection(double a) //JCD : does not work for now
{
    while (true)
    {
        double Y = G_inv(a,Random::Uniform());
        if (Random::Uniform()<=q(a,Y))
            return Y;
    }
}

void Brownian(vector<double>& traj, double num_steps, double step_size)
{
    traj[0] = 0.0;
    double vol_one_step = sqrt(step_size);
    for (int i=1;i<num_steps;++i)
    {
        traj[i] = traj[i-1] + vol_one_step * Normale(0.0,1.0);
    }
}

double Brownian_Square_Variation(int N, double T)
{
    vector<double> traj = vector<double>(N,0.0);
    Brownian(traj,N, T/N);
    double res = 0.0;
    for (int i=0;i<N-1;++i)
    {
        double increase = traj[i+1]-traj[i];
        res += increase*increase;
    }
    return res;
}
