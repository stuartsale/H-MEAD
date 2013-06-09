
//#include "helper.cpp"
//#include "bin_obj.h"
#include "LF.h"
#include "cat_read.h"
#include <Eigen/SparseCore>
#include <Eigen/Dense>


using namespace std;

#ifndef PI
#define PI 3.14159265358979323846
#endif

#ifndef fBm_s
#define fBm_s 0.555
#endif

#ifndef fBm_x2
#define fBm_x2 0.622
#endif

class sl_obj
{

	public:
		sl_obj(string filename, float l_in, float b_in, string datatype);
		sl_obj(void);
		void initial_guess(vector<iso_obj> &isochrones, vector<iso_obj> &guess_set, vector <LF> &LFs);
		void dist_redMCMC(vector<iso_obj> &isochrones, vector <LF> &LFs);
		void update(vector<iso_obj> &isochrones,  vector <LF> &LFs);
		void hyperprior_update(void);
		void mean_intervals(void);
		void output_write(void);
		void neighbour_set(sl_obj * neighbour);
		void acl_calc(void);
		float it_num;


	private:
		bool move_on;
		vector <iphas_obj> star_cat;

		float l, b;

		float r_min, i_min, ha_min;	// iphas min & maxes
		float r_max, i_max, ha_max;

		float J_min, H_min, K_min;	// 2MASS min & maxes
		float J_max, H_max, K_max;

		string rootname;

		vector < vector <float> > proposal_sd ;//(150, vector <float> (2));
		vector < vector <float> > previous_rel; //(150, vector <float> (4));
		vector < vector <float> > new_rel; //(150, vector <float> (4));
		vector < vector <float> > internal_rel; //(150, vector <float> (2));
		vector < vector <float> > previous_internal_rel; //(150, vector <float> (2));
		vector < vector <float> > hyperprior_internal_rel; //(150, vector <float> (2));
		vector < vector <float> > hyperprior_rel; //(150, vector <float> (2));

		vector < vector <float> > gen_internal_rel(vector < vector <float> > old_rel, int rel_length);
		vector < vector <float> > internal_to_external(vector < vector <float> > int_rel, int rel_length);
		float hyperprior_prob_get(vector < vector <float> > internal_rel);

		float global_previous_prob;
		float previous_hyperprior_prob, current_hyperprior_prob;
		float previous_norm_prob, current_norm_prob;
		float global_current_prob, global_transition_prob;
		float previous_xsl_prob, current_xsl_prob;

		// hold chains too
		vector <vector <vector <float> > > global_A_chain;

		vector<bin_obj2> A_mean;
		vector<float> backup_A_mean;
		vector < vector <float> > rho_mean;

		float sigma_fac, accepted;
// Set up

		int without_change;
		int thin;

		float sigma2_LN, mu_LN;
		vector <float> proposed_probs;
		float rel_length;



		sl_obj * neighbour_sl;
		vector <vector <float> > recv_neighbour_rel;

	// Covariance matrix

		Eigen::SparseMatrix<float> Cov_Mat; 

		Eigen::Matrix<float, 150, 1> Mean_vec;

		void define_cov_mat(void);


	// Wider disc params

		float previous_s_R, previous_s_z, previous_A_0;
		vector <float> s_R_chain, s_z_chain, A_0_chain;
		float s_R_mean, s_z_mean, A_0_mean;



};
