
//#include "helper.cpp"
//#include "bin_obj.h"
#include "LF.h"
#include "cat_read.h"
#include <Eigen/Sparse>
#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>




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
		sl_obj(string filename, float l_in, float b_in, string datatype, float s_R, float s_z);
		sl_obj(void);
		void initial_guess(vector<iso_obj> &isochrones, vector<iso_obj> &guess_set, vector <LF> &LFs, float s_R, float s_z, float A_0);
		void dist_redMCMC(vector<iso_obj> &isochrones, vector <LF> &LFs);
		void update(vector<iso_obj> &isochrones,  vector <LF> &LFs);
		void mean_intervals(void);
		void output_write(float s_R, float s_z);
		void neighbour_set(sl_obj * neighbour);
		void acl_calc(void);

		float it_num;


	private:
		bool move_on;
		float threshold;
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
		vector < vector <float> > trial_rel; //(150, vector <float> (2));

		vector < vector <float> > gen_internal_rel(vector < vector <float> > old_rel, int rel_length);
		vector < vector <float> > mvn_gen_internal_rel(void);
		vector < vector <float> > mvn_gen_internal_rel_no_neighbour(void);
		vector <float> mvn_gen_internal_rel_from_z_dash(const Eigen::Matrix<float, 150, 1> & z_dash);
		vector < vector <float> > internal_to_external(vector < vector <float> > int_rel, int rel_length);

		float get_rho_last_prob(void);
		float get_rho_test_prob(void);

		float hyperprior_prob_get(vector < vector <float> > internal_rel);
		void rho_to_A(void);
		void initial_rho_to_A(void);

		float global_previous_prob;
		float previous_hyperprior_prob, current_hyperprior_prob;
		float previous_norm_prob, current_norm_prob;
		float global_current_prob, global_transition_prob;
		float previous_xsl_prob, current_xsl_prob;

		double theta, theta_min, theta_max;

		// hold chains too
		vector <vector <vector <float> > > global_A_chain;

		vector<bin_obj> running_A_mean;
		vector<float> backup_A_mean;

		float sigma_fac, accepted;
// Set up

		int without_change;
		int thin;

		float sigma2_LN, mu_LN;
		vector <float> proposed_probs;
		float rel_length;



		vector <sl_obj *> neighbour_slsl;
		vector <vector <float> > recv_neighbour_rel;

	// Covariance matrix

		Eigen::SparseMatrix<float> Cov_Mat, Cov_Mat_Inv;
		Eigen::SparseMatrix<float> cond_Mat, cond_mu_Mat, rho_Mat;

		Eigen::SparseMatrix<float> chol_L, chol_L_Inv, chol_L_cond; 

		Eigen::Matrix<float, 150, 1> last_m_vec, test_m_vec;
		Eigen::Matrix<float, 150, 1> last_ln_vec, test_ln_vec;
		Eigen::Matrix<float, 150, 1> last_s_vec, test_s_vec;
		Eigen::Matrix<float, 150, 1> last_s_inv, test_s_inv;

		Eigen::Matrix<float, 150, 1> last_z_dash, test_z_dash;

		Eigen::Matrix<float, 150, 1> get_last_z_dash(void);

		void define_cov_mat(void);


	// Wider disc params

		void make_new_test_m_vec(float s_R, float s_z, float A_0);

	friend void hyperprior_update_all(vector <LF> &LFs);
	friend void neighbour_find(vector<sl_obj> &sl_list);

};
