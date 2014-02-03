#ifndef _BINOBJ_H_
#define _BINOBJ_H_


/*#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstdlib> 
#include <sstream>*/
#include <vector>

using namespace std;


//#include "iphas_obj.h"
//#include "iso_obj.h"
class iphas_obj;


class bin_obj2
{
	public:
	double mean_A, d_mean, sigma, d_sigma, size, error_measure, sum, diff, d_diff, median_A, median_sigma;
	bin_obj2();

	friend vector <bin_obj2> MCMC_bin(vector<iphas_obj>& aa, double r_i0_low, double r_i0_high, double r_min,double i_min,double ha_min,double r_max,double i_max, double ha_max);
	friend vector <bin_obj2> MCMC_bin_guts(vector<iphas_obj>& included, vector<bin_obj2> &backup_A_mean);
};
#endif
