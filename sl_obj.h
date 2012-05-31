
//#include "helper.cpp"
//#include "bin_obj.h"
#include "iphas_obj.h"
#include "cat_read.h"
using namespace std;

#ifndef PI
#define PI 3.14159265358979323846
#endif

class sl_obj
{

	public:
		sl_obj(string filename, double l_in, double b_in);
		sl_obj(void);
		void initial_guess(vector<iso_obj> &isochrones, vector<iso_obj> &guess_set);
		void dist_redMCMC(vector<iso_obj> &isochrones);
		void update(vector<iso_obj> &isochrones);
		void mean_intervals(void);
		void output_write(void);
		void neighbour_set(sl_obj * neighbour);
		float it_num;


	private:
		vector <iphas_obj> star_cat;

		double l, b;

		double r_min, i_min, ha_min;	// iphas min & maxes
		double r_max, i_max, ha_max;

		string rootname;

		vector < vector<double> > proposal_sd ;//(150, vector <double> (2));
		vector < vector <double> > previous_rel; //(150, vector <double> (4));
		vector < vector <double> > internal_rel; //(150, vector <double> (2));
		vector < vector <double> > previous_internal_rel; //(150, vector <double> (2));
		vector < vector <double> > first_internal_rel; //(150, vector <double> (2));

		double global_previous_prob;
		double previous_hyperprior_prob, current_hyperprior_prob;
		double global_current_prob, global_transition_prob;
		double previous_xsl_prob, current_xsl_prob;

		// hold chains too
		vector <vector <vector <double> > > global_A_chain;

		vector<bin_obj2> A_mean;
		vector<bin_obj2> backup_A_mean;				

		double sigma_fac, accepted;
// Set up

		int without_change;
		int thin;

		double sigma2_LN, mu_LN;
		vector <double> proposed_probs;
		double rel_length;



		sl_obj * neighbour_sl;

};
