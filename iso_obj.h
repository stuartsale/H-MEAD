/*#include <algorithm>
#include <cstdlib>*/
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <vector>
using namespace std;



//#ifndef PI
//#define PI 3.14159265358979323846
//#endif

//#ifndef ISOOBJ_H_
//#define IS0OBJ_H_
#pragma once
//#include "iphas_obj.h"
//#include "bin_obj.h"
//#include "helper.h"

//class iphas_obj;
class bin_obj2;

class iso_obj		// This is a class for holding isochrone points
{
	public:
	iso_obj(double feh_in, double Mi_in, double age_in, double logT_in, double logg_in, double r0_in, double i0_in, double ha0_in, double Jac_in);
	iso_obj(void);
	double IMF();
	double redline(double r_i1);
	
//	private:
	double Mi, logAge, logT, logg, feh, Jac;
	float r0, i0, ha0;
	float J0, H0, K0;

	double a,b,c, l,m,n, u,v,w, u_i,v_i,w_i, u_ha,v_ha,w_ha;	// terms governing response to extinction
	double u_J,v_J,w_J, u_H,v_H,w_H, u_K,v_K,w_K;

	
	
//	friend	void iphas_obj::red_dist(vector<bin_obj2> &A_mean, vector<iso_obj> &isochrones, double l, double b);
};

class LF		// A class for holding LFs
{
	public:
	LF(string filename);
	float LF_prob(vector < vector <float> > new_rel);
	float metal_prob;

	private:
	vector <vector <float> > LF_vec;
	void metal_prob_set(void);
};

extern vector <vector <vector <double> > > lookup_table;
vector <vector <vector <double> > > lookup_creator(void);
double int_lookup(double A_max, double A_mean, double sd);
struct params_struct {double A_max; double A_mean; double sigma;};

//#endif
