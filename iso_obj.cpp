#include "iso_obj.h"

// constructor
iso_obj::iso_obj(double feh_in, double Mi_in, double age_in, double logT_in, double logg_in, double r0_in, double i0_in, double ha0_in, double Jac_in)
{
	Mi=Mi_in;
	logAge=age_in;
	logT=logT_in;
	logg=logg_in;
	feh=feh_in;

	r0=r0_in;
	i0=i0_in;
	ha0=ha0_in;

	Jac=Jac_in;
	
   
//R=3.1

	u=(0.00043*(r0-i0)-0.002887);						// 2nd order on r 				
	v=(-0.02976*(r0-i0)+0.83778);							// 1st order on r

	u_i=(0.000106*(r0-i0)-0.001506) ;
	v_i=(-0.01380*(r0-i0)+0.59831);

	u_ha=(0.0000042*(r0-i0)-0.0000161);
	v_ha=(-0.0000063*(r0-i0)+0.761855);

}

iso_obj::iso_obj(void)
{
	Mi=1;
	logAge=9;
	logT=3.75;
	logg=4.5;
	feh=0;

	r0=0;
	i0=0;
	ha0=0;

	Jac=0;
}

double iso_obj::IMF(void)
{
	return pow(Mi, -2.7);
}

double iso_obj::redline(double r_i1)
{
	double A_int;
	A_int=0.;// (-(v-v_i)+sqrt(pow(v-v_i,2)-4*(u-u_i)*((w-w_i)+(r0-i0)-r_i1)))/(2*(u-u_i));  //(u-u_i, v-v_i, (w-w_i)+(r0-i0)-r_i1, +1);
	A_int=(-(r0-i0)+r_i1)/(v-v_i);
	return (u-u_ha)*pow(A_int,2) + (v-v_ha)*A_int + (r0-ha0);
}	
