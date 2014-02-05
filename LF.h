#include "iphas_obj.h"
#include <sstream>

class LF		// A class for holding LFs
{
	public:
	LF(string filename);
	float LF_prob(vector < vector <float> > new_rel);
	float LF_prob2(vector < vector <float> > new_rel);
	float feh;
	float metal_prob;

	private:
	vector <vector <float> > LF_vec;
	void metal_prob_set(void);
	vector <float> prior_lf;
	float norm_lf;
};
