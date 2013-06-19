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
#include <algorithm>

using namespace std;


//#include "iphas_obj.h"
//#include "iso_obj.h"
class iphas_obj;


class bin_obj
{

//	friend vector <bin_obj> MCMC_bin(vector<iphas_obj>& aa, double r_i0_low, double r_i0_high, double r_min,double i_min,double ha_min,double r_max,double i_max, double ha_max);
//	friend vector <bin_obj> MCMC_bin_guts(vector<iphas_obj>& included, vector<bin_obj> &backup_A_mean);

	public:
		double mean_A, sigma;
		bin_obj();

		void initial_add(iphas_obj * star);
	
		void try_add(iphas_obj * star);
		void try_remove(iphas_obj * star);
		void accept(void);		

//	private:

		vector<iphas_obj *>::iterator it_found;
		
		int last_n, test_n;


};

class bin_obj2
{

//	friend vector <bin_obj> MCMC_bin(vector<iphas_obj>& aa, double r_i0_low, double r_i0_high, double r_min,double i_min,double ha_min,double r_max,double i_max, double ha_max);
//	friend vector <bin_obj> MCMC_bin_guts(vector<iphas_obj>& included, vector<bin_obj> &backup_A_mean);


	public:
		double mean_A, d_mean, sigma, d_sigma, size, error_measure, sum, diff, d_diff;
		bin_obj2();

	private:
		vector <iphas_obj *> last_star_list;
		vector <iphas_obj *> next_star_list;


};
#endif
