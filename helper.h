#ifndef HELPER_H_
#define HELPER_H_


/*#include <iostream>
#include <string>
#include <algorithm>
#include <cstdlib> 
#include <sstream>*/
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include </home/stuart/work/work-chile/distance_red/MEAD-files/newran03/newran.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


using namespace std;
//#ifndef PI
//#define PI 3.14159265358979323846
//#endif

//GLOBALS (checked)

//#include "bin_obj.h"
//#include "iphas_obj.h"
#include "iso_obj.h"
//class bin_obj;

//
// HELPER FUNCTIONS
//


//Mag giver
//double mag_giver(vector< vector <calibration_obj> > &cdata, double r_i1, int lumclass, string band);


// Solves a quadratic eqn
double quadratic(double a, double b, double c, int sign);	

iso_obj iso_get(double targ_feh, double targ_Mi, double targ_logAge, vector<iso_obj> &isochrones);
iso_obj iso_get_Tg(double targ_feh, double targ_logT, double targ_logg, vector<iso_obj> &isochrones);
double max_age(double targ_Mi, vector<iso_obj> &isochrones);

// Setting up random number generators

extern Uniform U;
extern Normal Z;
extern VariLogNormal VLN;

extern gsl_rng* rng_handle;

// CDF of a normal distribution
double cdf_normal (double x, double mu, double sigma);	
double inv_cdf_normal_fast(double p, double mu, double sigma);
double cdf_normal_fast(double x, double mu, double sigma);
double cdf_normal_smallx(double x, double mu, double sigma);


// set default MIN vals
extern double r_min;
extern double i_min;
extern double ha_min;				
// set default MAX vals
extern double r_max;
extern double i_max;
extern double ha_max;

#endif
