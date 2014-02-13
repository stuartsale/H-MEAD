
//#include "helper.cpp"
//#include "bin_obj.h"
#include "LF.h"
#include "cat_read.h"
using namespace std;

#ifndef PI
#define PI 3.14159265358979323846
#endif

extern string hmd_dir;

class sl_obj
{

	public:
		sl_obj(string filename, float l_in, float b_in, string datatype, float r_max_in, float i_max_in, float ha_max_in);
		sl_obj(void);
		void initial_guess(vector<iso_obj> &isochrones, vector<iso_obj> &guess_set, vector <LF> &LFs);
		void dist_redMCMC(vector<iso_obj> &isochrones, vector <LF> &LFs);
		void update(vector<iso_obj> &isochrones,  vector <LF> &LFs);
		void mean_intervals(void);
		void output_write(void);
		void neighbour_set(sl_obj * neighbour);
		void acl_calc(void);
		float it_num;
		float r_max, i_max, ha_max;

		float l, b;

	private:
		vector <iphas_obj> star_cat;



		float r_min, i_min, ha_min;	// iphas min & maxes


		float J_min, H_min, K_min;	// 2MASS min & maxes
		float J_max, H_max, K_max;

		string rootname;

		vector < vector<float> > proposal_sd ;//(150, vector <float> (2));
		vector < vector <float> > previous_rel; //(150, vector <float> (4));
		vector < vector <float> > internal_rel; //(150, vector <float> (2));
		vector < vector <float> > previous_internal_rel; //(150, vector <float> (2));
		vector < vector <float> > first_internal_rel; //(150, vector <float> (2));
		vector < vector <float> > first_rel; //(150, vector <float> (2));

		float global_previous_prob;
		float previous_hyperprior_prob, current_hyperprior_prob;
		float last_part_prior2, test_part_prior2;
		float last_part_prior1, test_part_prior1;
		float previous_norm_prob, current_norm_prob;
		float global_current_prob, global_transition_prob;
		float previous_xsl_prob, current_xsl_prob;

		// hold chains too
		vector <vector <vector <float> > > global_A_chain;

		vector<bin_obj2> A_mean;
		vector<bin_obj2> backup_A_mean;	

		vector < vector<float> > rho_final;		

		float sigma_fac, accepted;
// Set up

		int without_change;
		int thin;

		float sigma2_LN, mu_LN;
		vector <float> proposed_probs;
		size_t rel_length;



		sl_obj * neighbour_sl;
		vector <vector <float> > recv_neighbour_rel;

};
