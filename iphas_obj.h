#ifndef IPHASOBJ_H_
#define IPHASOBJ_H_



/*
#include <fstream>
#include <algorithm>
#include <cstdlib> 
#include <sstream>*/
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
using namespace std;

#ifndef PI
#define PI 3.14159265358979323846
#endif


//#endif
//#include "iso_obj.h"
//#include "bin_obj.h"
#include "helper.h"
//class bin_obj2;
//class iphas_sub_obj;
//class iso_obj;

// This is a class to hold the photometry
// and estimated properties of observed objects.
// One of these is created for each object in 
// the IPHAS catalogue read in.
class iphas_obj
{	

	friend class bin_obj2;
	friend class sl_obj;
	public:
		iphas_obj(double r_input, double i_input, double ha_input, double d_r_input, double d_i_input, double d_ha_input, double l_input, double b_input);
      		iphas_obj(double d_r_input, double d_i_input, double d_ha_input);
		iphas_obj(double r_input, double i_input, double ha_input, double d_r_input, double d_i_input, double d_ha_input, double l_input, double b_input, double real_dist_in, double real_A_in, double real_Mi_in, double real_logAge_in, double real_feh_in);
		void red_dist(vector<bin_obj2> &A_mean, vector<iso_obj> &isochrones, double l, double b);
		void red_dist(vector<iso_obj> &isochrones, double l, double b);
		void print();
		void set_mag_weight(double mag_weight_input);
		void dist_redMCMC(vector<iso_obj> &isochrones, vector<iso_obj> &guess_set, double l, double b);
		void dist_redMCMC(vector<iso_obj> &isochrones, vector<iso_obj> &guess_set, double l, double b, vector<bin_obj2> &A_mean);
		double prob_eval(iso_obj test_iso, double test_A, double test_dist_mod, vector<vector <double> > &A_mean);
		void new_func(vector<iso_obj> &isochrones, double l, double b, vector<bin_obj2> &A_mean, int flag);
		void initial_guess(vector<iso_obj> &isochrones, vector<iso_obj> &guess_set, vector<vector <double> > &A_mean);
		void star_try1(vector<iso_obj> &isochrones, double &l, double &b, vector<vector <double> > &A_mean);
		void mean_intervals(void);
		double get_A_prob(iso_obj test_iso, double test_A, double test_dist_mod, vector<vector <double> > &A_mean);
	//
	   
      // Returns abs mag or colour of an object given its r-i colour and luminosity class
      
     

     

   private:
   // The basic photometry of each object and its calculated parameters - with uncertainties
	double r, i, ha, d_r, d_i, d_ha, l, b;	
	double r_i0, A, dist, d_r_i0, d_A, d_dist, mag_weight;
	double Mi, logAge, feh, d_Mi, d_logAge, d_feh;
	double logT, logg;
	int distbin;
	double A_min, A_max, A_min_r, A_max_r, A_min_i, A_max_i, A_min_ha, A_max_ha, A_prob;

	vector <iso_obj> iso_obj_chain;
	vector <double> dist_mod_chain;
	vector <double> A_chain;
	vector <double> prob_chain, A_prob_chain;
	int no_accept, no_accept2;

	double mean_prob, mean_A_prob;

	iso_obj last_iso;
	double last_dist_mod, last_A, last_logT, last_logg, last_dist;
	double last_prob, last_A_prob;
	double last_rmag, last_ri;

	double best_prob, best_A_prob;
	iso_obj best_iso;
	double best_A, best_dist_mod;
	double best_it;
	
	double logT_sd, logg_sd, feh_sd, A_sd, dist_mod_sd;
	double ri_sd, rmag_sd;


	double real_dist, real_A, real_Mi, real_logAge, real_feh;

	vector <double> rx_chain, ix_chain, hax_chain;
	double rx, ix, hax;

	vector<double> acl_calc(void);
	
   // BRIGHT LIMITS (default)
 //  double r_min, i_min, ha_min;				
   // FAINT LIMITS (default)
  // double r_max, i_max, ha_max;

		
   //friends!

	friend vector <bin_obj2>  dist_redMCMC(vector<iphas_obj> &stars, vector<iso_obj> &isochrones, vector<iso_obj> &guess_set, double l, double b, vector <bin_obj2> backup_A_mean, double ri_min, double ri_max);
	friend void output_write(string filename, vector<bin_obj2> A_mean, vector<iphas_obj> colours);
	friend double real_prob(vector<iphas_obj> &stars, vector<iso_obj> &isochrones, vector<iso_obj> &guess_set, double l, double b, vector <bin_obj2> backup_A_mean, double ri_min, double ri_max);
	//friend bool bin_obj::obj_comparison(iphas_obj object1, iphas_obj object2);
//
};

#endif
