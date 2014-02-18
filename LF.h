#include "iphas_obj.h"
#include <sstream>

class LF		// A class for holding LFs
{
	public:
	LF(string filename, float feh);
	float LF_prob2(vector < vector <float> > new_rel);
	void set_prior_lf(float l_in, float b_in, float dl, float db);
	float feh;
	float metal_prob;
	void precompute_Aminmax(void);
	float alpha, beta;

	private:
	float l,b;
	vector <vector <float> > LF_vec;
	void metal_prob_set(void);
	vector <float> prior_lf;
	float norm_lf;
	vector < vector <float> > A_min;
	vector < vector <float> > A_max;

	float gauss_table[801];

};
