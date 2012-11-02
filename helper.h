#ifndef HELPER_H_
#define HELPER_H_


/*#include <iostream>*/
#include <string>
#include <algorithm>
/*#include <cstdlib> 
#include <sstream>*/
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>



using namespace std;
#ifndef PI
#define PI 3.14159265358979323846
#endif

//GLOBALS (checked)

#include "bin_obj.h"
//#include "iphas_obj.h"
#include "iso_obj.h"
#include <mpi.h>
//class bin_obj;

//
// HELPER FUNCTIONS
//

// Setting up random number generators

extern gsl_rng* rng_handle;

// Solves a quadratic eqn
double quadratic(double a, double b, double c, int sign);	

iso_obj iso_get(double targ_feh, double targ_Mi, double targ_logAge, vector<iso_obj> &isochrones);
iso_obj iso_get_Tg(double targ_feh, double targ_logT, double targ_logg, vector<iso_obj> &isochrones);
double max_age(double targ_Mi, vector<iso_obj> &isochrones);


string getStdoutFromCommand(string cmd);

void output_write(string filename, vector<bin_obj2> A_mean, vector<iphas_obj> colours);
string stringify(double x);
double StrToDbl(string s);
vector<float> backup_A_mean_find(double l_gal, double b_gal, float s_R, float s_z, bool Sch_get)	;
vector<float> density_find(double l_gal, double b_gal, float s_R, float s_z, float res)	;

vector <vector <string> > config_read(string filename);

// set default MIN vals
extern double r_min;
extern double i_min;
extern double ha_min;
extern double J_min;
extern double H_min;
extern double K_min;			
// set default MAX vals
extern double r_max;
extern double i_max;
extern double ha_max;
extern double J_max;
extern double H_max;
extern double K_max;

extern vector <vector <vector <double> > > lookup_table;
vector <vector <vector <double> > > lookup_creator(void);
double int_lookup(double A_max, double A_mean, double sd);
struct params_struct {double A_max; double A_mean; double sigma;};

#endif
