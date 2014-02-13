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

#ifndef LN_TEN
#define LN_TEN 2.302585093
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
		iphas_obj(float r_input, float i_input, float ha_input, float d_r_input, float d_i_input, float d_ha_input, float l_input, float b_input);
      		iphas_obj(float d_r_input, float d_i_input, float d_ha_input);
		iphas_obj(float r_input, float i_input, float ha_input, float d_r_input, float d_i_input, float d_ha_input, float l_input, float b_input, float real_dist_in, float real_A_in, float real_Mi_in, float real_logAge_in, float real_feh_in);
		iphas_obj(float P1_input, float P2_input, float P3_input, float d_P1_input, float d_P2_input, float d_P3_input, float l_input, float b_input, string source);
		void red_dist(vector<bin_obj2> &A_mean, vector<iso_obj> &isochrones, float l, float b);
		void red_dist(vector<iso_obj> &isochrones, float l, float b);
		void print();
		void set_mag_weight(float mag_weight_input);
		void dist_redMCMC(vector<iso_obj> &isochrones, vector<iso_obj> &guess_set, float l, float b);
		void dist_redMCMC(vector<iso_obj> &isochrones, vector<iso_obj> &guess_set, float l, float b, vector<bin_obj2> &A_mean);
		float likelihood_eval(iso_obj test_iso, float test_A, float test_dist_mod, vector<vector <float> > &A_mean);
		void new_func(vector<iso_obj> &isochrones, float l, float b, vector<bin_obj2> &A_mean, int flag);
		void initial_guess(vector<iso_obj> &isochrones, vector<iso_obj> &guess_set, vector<vector <float> > &A_mean);
		void star_try1(vector<iso_obj> &isochrones, float &l, float &b, vector<vector <float> > &A_mean);
		void mean_intervals(void);
		float get_A_prob(iso_obj test_iso, float test_A, float test_dist_mod, vector<vector <float> > &A_mean);
		void adaptive_proposal_update(int batch_len);
	//
	   
      // Returns abs mag or colour of an object given its r-i colour and luminosity class
      
     

     

   private:
   // The basic photometry of each object and its calculated parameters - with uncertainties
	float r, i, ha, d_r, d_i, d_ha;
	float J, H, K, d_J, d_H, d_K;
	float l, b;
	int pixel;	
	float r_i0, A, dist, d_r_i0, d_A, d_dist, mag_weight;
	float Mi, logAge, feh, d_Mi, d_logAge, d_feh;
	float logT, logg;
	int distbin;
	float A_min, A_max, A_min_r, A_max_r, A_min_i, A_max_i, A_min_ha, A_max_ha, A_prob;
	float A_min_J, A_max_J, A_min_H, A_max_H, A_min_K, A_max_K;
	float update_prop;

	vector <iso_obj> iso_obj_chain;
	vector <float> dist_mod_chain;
	vector <float> A_chain;
	vector <float> prob_chain, A_prob_chain;
	float no_accept, no_accept2, batch_accept;

	float mean_prob, mean_A_prob;

	iso_obj last_iso;
	float last_dist_mod, last_A, last_logT, last_logg, last_dist;
	float last_prob, last_A_prob;
	float last_rmag, last_ri;

	float best_prob, best_A_prob;
	iso_obj best_iso;
	float best_A, best_dist_mod;
	float best_it;
	
	float logT_sd, logg_sd, feh_sd, A_sd, dist_mod_sd;
	float ri_sd, rmag_sd;


	float real_dist, real_A, real_Mi, real_logAge, real_feh;

	vector <float> rx_chain, ix_chain, hax_chain;
	float rx, ix, hax;

	vector<float> acl_calc(void);
	
   // BRIGHT LIMITS (default)
 //  float r_min, i_min, ha_min;				
   // FAINT LIMITS (default)
  // float r_max, i_max, ha_max;

		
   //friends!

	friend vector <bin_obj2>  dist_redMCMC(vector<iphas_obj> &stars, vector<iso_obj> &isochrones, vector<iso_obj> &guess_set, float l, float b, vector <bin_obj2> backup_A_mean, float ri_min, float ri_max);
	friend void output_write(string filename, vector<bin_obj2> A_mean, vector<iphas_obj> colours);
	friend float real_prob(vector<iphas_obj> &stars, vector<iso_obj> &isochrones, vector<iso_obj> &guess_set, float l, float b, vector <bin_obj2> backup_A_mean, float ri_min, float ri_max);
	//friend bool bin_obj::obj_comparison(iphas_obj object1, iphas_obj object2);
//
};

float log_prior(float test_dist_mod, float test_feh, float l, float b);

#endif
