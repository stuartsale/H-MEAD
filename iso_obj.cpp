#include "iso_obj.h"

// constructor
iso_obj::iso_obj(double feh_in, double Mi_in, double age_in, double logT_in, double logg_in, double r0_in, double i0_in, double ha0_in, double J0_in, double H0_in, double K0_in, double Jac_in)
{
	Mi=Mi_in;
	logAge=age_in;
	logT=logT_in;
	logg=logg_in;
	feh=feh_in;

	r0=r0_in;
	i0=i0_in;
	ha0=ha0_in;

	J0=J0_in;
	H0=H0_in;
	K0=K0_in;

	Jac=Jac_in;
	
//R=3.1

	u=(0.0007075*(r0-i0)-0.00504174) /1.08573;						// 2nd order on r 				
	v=(-0.0396596*(r0-i0)+1.1106105) /1.08573;							// 1st order on r
	w=(0.001023008*(r0-i0)+0.0008454) /1.08573;							// 0th order on r	

	u_i=(0.00004671*(r0-i0)-0.0025414) /1.08573;
	v_i=(-0.0169457*(r0-i0)+0.7929863) /1.08573;
	w_i=(0.0011835*(r0-i0)-0.0037645) /1.08573;

	u_ha=(0.00007665*(r0-i0)-0.00007222) /1.08573;
	v_ha=(-0.00054454*(r0-i0)+1.003924) /1.08573;
	w_ha=(-0.0004797*(r0-i0)+0.00028745) /1.08573;

	u_J=0;
	v_J=0.276;
	w_J=0;

	u_H=0;
	v_H=0.176;
	w_H=0;

	u_K=0;
	v_K=0.112;
	w_K=0;

}

iso_obj::iso_obj(double feh_in, double Mi_in, double age_in, double logT_in, double logg_in, double P1_in, double P2_in, double P3_in, double Jac_in, string source)
{
	Mi=Mi_in;
	logAge=age_in;
	logT=logT_in;
	logg=logg_in;
	feh=feh_in;

	if (source=="IPHAS")
	{
		r0=P1_in;
		i0=P2_in;
		ha0=P3_in;
	}
	else if (source=="2MASS")
	{
		J0=P1_in;
		H0=P2_in;
		K0=P3_in;
	}

	Jac=Jac_in;
	
//R=3.1

	u=(0.0007075*(r0-i0)-0.00504174) /1.08573;							// 2nd order on r 				
	v=(-0.0396596*(r0-i0)+1.1106105) /1.08573;							// 1st order on r
	w=(0.001023008*(r0-i0)+0.0008454) /1.08573;							// 0th order on r	

	u_i=(0.00004671*(r0-i0)-0.0025414) /1.08573;
	v_i=(-0.0169457*(r0-i0)+0.7929863) /1.08573;
	w_i=(0.0011835*(r0-i0)-0.0037645) /1.08573;

	u_ha=(0.00007665*(r0-i0)-0.00007222) /1.08573;
	v_ha=(-0.00054454*(r0-i0)+1.003924) /1.08573;
	w_ha=(-0.0004797*(r0-i0)+0.00028745) /1.08573;

	u_J=0;
	v_J=0.276;
	w_J=0;

	u_H=0;
	v_H=0.176;
	w_H=0;

	u_K=0;
	v_K=0.112;
	w_K=0;


/*//R=2.9

	u=(0.000825*(r0-i0)-0.005582) /1.08573;								// 2nd order on r 				
	v=(-0.04209*(r0-i0)+1.11201874) /1.08573;							// 1st order on r
	w=(0.00115325*(r0-i0)+0.0010913316) /1.08573;							// 0th order on r	

	u_i=(-0.0000659*(r0-i0)-0.002717) /1.08573;
	v_i=(-0.017621*(r0-i0)+0.77944) /1.08573;
	w_i=(0.001279*(r0-i0)-0.0003786) /1.08573;

	u_ha=(0.00007835*(r0-i0)-0.00007595) /1.08573;
	v_ha=(-0.0005294*(r0-i0)+0.99945) /1.08573;
	w_ha=(-0.0005279*(r0-i0)+0.0003166) /1.08573; */



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
	return (u-u_ha)*pow(A_int,2) + (v-v_ha)*A_int + (w-w_ha)+(r0-ha0);
}	
