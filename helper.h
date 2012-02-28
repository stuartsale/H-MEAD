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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>



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

double in_sample(double r_in, double i_in, double ha_in, double dr_in, double di_in, double dha_in);
double in_sample2(double A_max, double mu, double sigma);

// Setting up random number generators

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

struct A_params {double A_max; double A; double sigma;};
double A_integral_func (double *mean, size_t dim, void *params);

#endif
