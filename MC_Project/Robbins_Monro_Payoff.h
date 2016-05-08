#ifndef Robbins_Monro_Payoff_h
#define Robbins_Monro_Payoff_h

#include <boost/numeric/ublas/vector.hpp>
typedef boost::numeric::ublas::vector<double> Vector;

class Robbins_Monro_Payoff {
public:
	virtual double Payoff(const Vector & G) const = 0;
	virtual double StockBS(const boost::numeric::ublas::vector<double> & G) const = 0;

	double T;
	double r;
};

#endif