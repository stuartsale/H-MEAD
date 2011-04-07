/*#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdlib> 
#include <sstream>*/
#include <cmath>
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
	iso_obj(double feh_in, double Mi_in, double age_in, double Teff_in, double logg_in, double r0_in, double i0_in, double ha0_in, double Jac_in);
	iso_obj(void);
	double IMF();
	double redline(double r_i1);
	
//	private:
	double Mi, logAge, Teff, logg, feh, r0, i0, ha0, Jac;
	double a,b,c, l,m,n, u,v,w, u_i,v_i,w_i, u_ha,v_ha,w_ha;	// terms governing response to extinction

	
	
//	friend	void iphas_obj::red_dist(vector<bin_obj2> &A_mean, vector<iso_obj> &isochrones, double l, double b);
};

//#endif
