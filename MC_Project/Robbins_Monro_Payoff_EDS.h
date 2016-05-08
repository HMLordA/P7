#ifndef Robbins_Monro_Payoff_EDS_h
#define Robbins_Monro_Payoff_EDS_h

#include <boost/numeric/ublas/vector.hpp>
typedef boost::numeric::ublas::vector<double> Vector;

class Robbins_Monro_Payoff_EDS {

public:
	typedef std::pair<double, double> state;
	typedef std::list<state> diffusion;
	virtual double Payoff_Call(const diffusion & eds) const = 0;

	double T;
	double r;
};

#endif