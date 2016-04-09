//
//  Theta.c
//  MC_git_Project
//
//  Created by Nazar KOSTYUCHYK on 25/03/2016.
//  Copyright © 2016 Nazar KOSTYUCHYK. All rights reserved.
//

#include "Theta.h"

static double* pow_array = 0;
static double* half_pow_array = 0;

double Phi_n_k(double n, double k, double t)
{
	if (pow_array == 0)
	{
		double nb_power_2 = DIM_THETA+1;
		pow_array = new double[(int)nb_power_2];
		half_pow_array = new double[(int)nb_power_2];
		pow_array[0] = 1.0;
		half_pow_array[0] = 1.0;
		for (int i=1;i<(int)nb_power_2;++i)
		{
			pow_array[i] = 2.0*pow_array[i-1];
			half_pow_array[i] = sqrt(2.0)*half_pow_array[i-1];
		}
	}
	double phi = 0.0;
	//double t_tilda = pow_array[(int)k]*t-n;
	double t_tilda = pow_array[(int)n]*t-k;
	if (t_tilda>=0 && t_tilda<0.5)
		phi = 1.0;
	else if (t_tilda>=0.5 && t_tilda<1)
		phi = -1.0;
	
	//if (phi !=0.0) 
		return half_pow_array[(int)n]*phi;
	//else
	//	return 0.0;
}
