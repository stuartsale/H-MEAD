#include "iphas_obj.h"


iphas_obj::iphas_obj(double r_input, double i_input, double ha_input, double d_r_input,double d_i_input, double d_ha_input, double l_input, double b_input)
{

	r=r_input;	
	i=i_input;
	ha=ha_input;
	d_r=d_r_input;
	d_i=d_i_input;
	d_ha=d_ha_input;
	l=l_input;
	b=b_input;

	feh_sd=0.008;
	logT_sd=0.003;
	logg_sd=0.008;
	A_sd=sqrt(d_r*d_r+d_i*d_i)/4;
	dist_mod_sd=0.5*d_r;

	ri_sd=min(sqrt(d_r*d_r+d_i*d_i), 0.04);
	rmag_sd=d_r;

	no_accept=0;
}

iphas_obj::iphas_obj(double r_input, double i_input, double ha_input, double d_r_input,double d_i_input, double d_ha_input, double l_input, double b_input, double real_dist_in, double real_A_in, double real_Mi_in, double real_logAge_in, double real_feh_in)
{

	r=r_input;	
	i=i_input;
	ha=ha_input;
	d_r=d_r_input;
	d_i=d_i_input;
	d_ha=d_ha_input;
	l=l_input;
	b=b_input;

	feh_sd=0.008;
	logT_sd=0.003;
	logg_sd=0.008;
	A_sd=sqrt(d_r*d_r+d_i*d_i)/4;
	dist_mod_sd=0.5*d_r;

	ri_sd=min(sqrt(d_r*d_r+d_i*d_i), 0.04);
	rmag_sd=d_r;

	no_accept=0;
	no_accept2=0;

	real_dist=real_dist_in;
	real_A=real_A_in;
	real_Mi=real_Mi_in;
	real_logAge=real_logAge_in;
	real_feh=real_feh_in;
}

iphas_obj::iphas_obj(double d_r_input, double d_i_input, double d_ha_input)
{
	d_r=d_r_input;
	d_i=d_i_input;
	d_ha=d_ha_input;
}



void iphas_obj::set_mag_weight(double mag_weight_input)
{
   mag_weight=mag_weight_input;
}

double iphas_obj::prob_eval(iso_obj test_iso, double test_A, double test_dist_mod, vector<vector <double> > &A_mean)
{	
// Find probability of this parameter set - distance
	//cout << test_iso.Mi << " " << test_iso.logAge << " " << test_A << " " << test_dist_mod << " " ;
	// Find P(S|y,x,sigma_y)
		double current_prob1=0;
		//current_prob1-=log(1+exp(3*((test_iso.u*pow(test_A,2)+test_iso.v*test_A+test_iso.w)+test_dist_mod+test_iso.r0-r_max+0.5)));
		//current_prob1-=log(1+exp(3*((test_iso.u_i*pow(test_A,2)+test_iso.v_i*test_A+test_iso.w_i)+test_dist_mod+test_iso.i0-i_max+0.5)));
		//current_prob1-=log(1+exp(3*((test_iso.u_ha*pow(test_A,2)+test_iso.v_ha*test_A+test_iso.w_ha)+test_dist_mod+test_iso.ha0-ha_max+0.5)));

		//current_prob1+=log(in_sample((test_iso.u*pow(test_A,2)+test_iso.v*test_A+test_iso.w)+test_dist_mod+test_iso.r0, (test_iso.u_i*pow(test_A,2)+test_iso.v_i*test_A+test_iso.w_i)+test_dist_mod+test_iso.i0, (test_iso.u_ha*pow(test_A,2)+test_iso.v_ha*test_A+test_iso.w_ha)+test_dist_mod+test_iso.ha0, d_r, d_i, d_ha));
		//cout << r << 

	// Find p(y|x,sigma_y) 
		current_prob1+=-pow(r-(test_iso.u*pow(test_A,2)+test_iso.v*test_A+test_iso.w)-test_dist_mod-test_iso.r0,2)/(2*d_r*d_r) ;
	//	current_prob1+=+1.0857*((test_iso.u*pow(test_A,2)+test_iso.v*test_A+test_iso.w)+test_dist_mod+test_iso.r0) - pow((1-pow(10,-0.4*(r-(test_iso.u*pow(test_A,2)+test_iso.v*test_A+test_iso.w)-test_dist_mod-test_iso.r0))/d_r)*1.0857,2)/2;
//	cout<< current_prob1 << " " << d_r << " "  << r <<" " <<(test_iso.u*pow(test_A,2)+test_iso.v*test_A+test_iso.w)+test_dist_mod+test_iso.r0 << " " << test_dist_mod+test_iso.r0  << endl;
		current_prob1+=	-pow(i-(test_iso.u_i*pow(test_A,2)+test_iso.v_i*test_A+test_iso.w_i)-test_dist_mod-test_iso.i0,2)/(2*d_i*d_i) ;
	//	current_prob1+=+1.0857*((test_iso.u_i*pow(test_A,2)+test_iso.v_i*test_A+test_iso.w_i)+test_dist_mod+test_iso.i0) - pow((1-pow(10,-0.4*(i-(test_iso.u_i*pow(test_A,2)+test_iso.v_i*test_A+test_iso.w_i)-test_dist_mod-test_iso.i0))/d_i)*1.0857,2)/2;
	//cout<< current_prob1 << " " << d_i << " " << i << " " << (test_iso.u_i*pow(test_A,2)+test_iso.v_i*test_A+test_iso.w_i)+test_dist_mod+test_iso.i0 << " ";
		current_prob1+=	-pow(ha-(test_iso.u_ha*pow(test_A,2)+test_iso.v_ha*test_A+test_iso.w_ha)-test_dist_mod-test_iso.ha0,2)/(2*d_ha*d_ha);
	//	current_prob1+=+1.0857*((test_iso.u_ha*pow(test_A,2)+test_iso.v_ha*test_A+test_iso.w_ha)+test_dist_mod+test_iso.ha0) - pow((1-pow(10,-0.4*(ha-(test_iso.u_ha*pow(test_A,2)+test_iso.v_ha*test_A+test_iso.w_ha)-test_dist_mod-test_iso.ha0))/d_ha)*1.0857,2)/2;
	//cout<< current_prob1 << " "  << d_ha << " " << ha << " " << (test_iso.u_ha*pow(test_A,2)+test_iso.v_ha*test_A+test_iso.w_ha)+test_dist_mod+test_iso.ha0 << " ";
	// Find p(x) - includes A(dist), IMF, SFH, disc density & metallicity profiles - remember dist^2 term

		double test_dist=pow(10,test_dist_mod/5+1);


		// IMF - Scalo type?
		current_prob1+=log(test_iso.IMF());
	//	current_prob1-=test_iso.logAge;



		double /*cosb=1, sinb=0, cosl=-1,*/ R_gal;
		R_gal=sqrt(test_dist*test_dist*cos(b)*cos(b) + 64000000-16000*test_dist*cos(l)*cos(b));

		// density profile
	//	current_prob1+=-R_gal/2500 -test_dist*sinb/200;
		if (R_gal<130000){current_prob1+=-R_gal/3000 -test_dist*sin(b)/200;}
		else {current_prob1+=-R_gal/1200 +6.5 -test_dist*sin(b)/200;} 		// -6.5=13000/3000 - 13000/1200

		// metallicity profile
		current_prob1+=-pow(test_iso.feh-(R_gal-8000.)*0.00007,2)/(2*0.5);
	//	current_prob1+=-pow(test_iso.feh,2)/(2*pow(0.05,2));
		



	//	// dist^2 term 
		current_prob1+=2*log(test_dist);

		current_prob1+=log(test_dist/(2*test_A*(test_iso.u-test_iso.u_i)+test_iso.v-test_iso.v_i));
		current_prob1+=log(test_iso.Jac);
//*/

	// Also the the contribution to p(x) from the distance-reddening relationship

		A_max_r=test_A+(r_max-r)/test_iso.v;
		A_min_r=test_A-(r-r_min)/test_iso.v;
		A_max_i=test_A+(i_max-i)/test_iso.v_i;
		A_min_i=test_A-(i-i_min)/test_iso.v_i;
		A_max_ha=test_A+(ha_max-ha)/test_iso.v_ha;
		A_min_ha=test_A-(ha-ha_min)/test_iso.v_ha;

		A_max=min(min(A_max_r,A_max_i),A_max_ha);
		A_min=max(max(A_min_r,A_min_i),A_min_ha);



		if (floor(test_dist*1.04/100)<A_mean.size())
		{
			A_prob=0;
			current_prob1+=-log(A_mean[floor(test_dist*1.04/100)][3]*test_A) - pow(log(test_A)-A_mean[floor(test_dist*1.04/100)][2],2)/(2*pow(A_mean[floor(test_dist*1.04/100)][3],2));

	// Correction to prior to account for incompletness due to mag limits

		//	A_prob=-log(in_sample2(A_max, A_mean[floor(test_dist/100)][2],A_mean[floor(test_dist/100)][3]));
	
		/*	if (A_min>0){A_prob=-log(cdf_normal_fast(log(A_max),A_mean[floor(test_dist/100)][2],A_mean[floor(test_dist/100)][3])-cdf_normal_fast(log(A_min),A_mean[floor(test_dist/100)][2],A_mean[floor(test_dist/100)][3]));}
			else {A_prob=-log(cdf_normal_fast(log(A_max),A_mean[floor(test_dist/100)][2],A_mean[floor(test_dist/100)][3]));}
			if (isinf(A_prob))
			{
				A_prob=-log(cdf_normal_smallx(log(A_max),A_mean[floor(test_dist/100)][2],A_mean[floor(test_dist/100)][3]));
				//cout << A_max << " " << A_min << " " << cdf_normal_smallx(log(A_max),A_mean[floor(test_dist/100)][2],A_mean[floor(test_dist/100)][3]) << " " << A_mean[floor(test_dist/100)][0] << endl;
			}*/
	
		//	current_prob1+=A_prob;
		}
		else
		{
			A_prob=0;
			current_prob1+=-log(A_mean[A_mean.size()-1][3]*test_A) - pow(log(test_A)-A_mean[A_mean.size()-1][2],2)/(2*pow(A_mean[A_mean.size()-1][3],2));

	// Correction to prior to account for incompletness due to mag limits

		//	A_prob=-log(in_sample2(A_max, A_mean[A_mean.size()-1][2],A_mean[A_mean.size()-1][3]));
	
		/*	if (A_min>0){A_prob=-log(cdf_normal_fast(log(A_max),A_mean[A_mean.size()-1][2],A_mean[A_mean.size()-1][3])-cdf_normal_fast(log(A_min),A_mean[A_mean.size()-1][2],A_mean[A_mean.size()-1][3]));}
			else {A_prob=-log(cdf_normal_fast(log(A_max),A_mean[A_mean.size()-1][2],A_mean[A_mean.size()-1][3]));}
			if (isinf(A_prob))
			{
				A_prob=-log(cdf_normal_smallx(log(A_max),A_mean[A_mean.size()-1][2],A_mean[A_mean.size()-1][3]));
				//cout << A_max << " " << A_min << " " << cdf_normal_smallx(log(A_max),A_mean[floor(test_dist/100)][2],A_mean[floor(test_dist/100)][3]) << " " << A_mean[floor(test_dist/100)][0] << endl;
			}*/

		//	current_prob1+=A_prob;
		}

	if (current_prob1!=current_prob1){/*cout<< test_dist/100<< " " << " " << A_prob<< " " << current_prob1 << "" "" << A_max <<endl;*/ current_prob1=-1E6;}

	return current_prob1;
} 

double iphas_obj::get_A_prob(iso_obj test_iso, double test_A, double test_dist_mod, vector<vector <double> > &A_mean)
{
	// Find P(S|y,x,sigma_y)
		double current_prob1=0;
/*		current_prob1-=log(1+exp(3*((test_iso.u*pow(test_A,2)+test_iso.v*test_A+test_iso.w)+test_dist_mod+test_iso.r0-r_max+0.5)));
		current_prob1-=log(1+exp(3*((test_iso.u_i*pow(test_A,2)+test_iso.v_i*test_A+test_iso.w_i)+test_dist_mod+test_iso.i0-i_max+0.5)));
		current_prob1-=log(1+exp(3*((test_iso.u_ha*pow(test_A,2)+test_iso.v_ha*test_A+test_iso.w_ha)+test_dist_mod+test_iso.ha0-ha_max+0.5)));
		//cout << r << */

	// Find p(y|x,sigma_y) 
		current_prob1+=-pow(r-(test_iso.u*pow(test_A,2)+test_iso.v*test_A+test_iso.w)-test_dist_mod-test_iso.r0,2)/(2*d_r*d_r) ;
//	cout<< current_prob1 << " " << d_r << " "  << r <<" " <<(test_iso.u*pow(test_A,2)+test_iso.v*test_A+test_iso.w)+test_dist_mod+test_iso.r0 << " " << test_dist_mod+test_iso.r0  << endl;
	//	current_prob1+=	- pow(i-(test_iso.u_i*pow(test_A,2)+test_iso.v_i*test_A+test_iso.w_i)-test_dist_mod-test_iso.i0,2)/(2*d_i*d_i) ;
	//cout<< current_prob1 << " " << d_i << " " << i << " " << (test_iso.u_i*pow(test_A,2)+test_iso.v_i*test_A+test_iso.w_i)+test_dist_mod+test_iso.i0 << " ";
	//	current_prob1+=	- pow(ha-(test_iso.u_ha*pow(test_A,2)+test_iso.v_ha*test_A+test_iso.w_ha)-test_dist_mod-test_iso.ha0,2)/(2*d_ha*d_ha);
	//cout<< current_prob1 << " "  << d_ha << " " << ha << " " << (test_iso.u_ha*pow(test_A,2)+test_iso.v_ha*test_A+test_iso.w_ha)+test_dist_mod+test_iso.ha0 << " ";
	// Find p(x) - includes A(dist), IMF, SFH, disc density & metallicity profiles - remember dist^2 term

//*/		// IMF - Scalo type?
//		current_prob1+=log(test_iso.IMF());
		//current_prob1+=test_iso.logAge;

		double test_dist=pow(10,test_dist_mod/5+1);

/*		double cosb=1, sinb=0, cosl=-1, R_gal;
		R_gal=sqrt(test_dist*test_dist*cosb*cosb + 64000000-16000*test_dist*cosl*cosb);

		// density profile
		current_prob1+=-R_gal/2500 -test_dist*sinb/200;
	//	if (R_gal<130000){current_prob1+=-R_gal/3000 -test_dist*sinb/200;}
	//	else {current_prob1+=-R_gal/1200 +6.5 -test_dist*sinb/200;} 		// -6.5=13000/3000 - 13000/1200

		// metallicity profile
	//	current_prob1+=-pow(test_iso.feh-(R_gal-8000.)*0.00007,2)/(2*0.5);
	//	current_prob1+=-pow(test_iso.feh,2)/(2*pow(0.05,2));

	//	// dist^2 term */

	//	current_prob1+=2*log(test_dist);

	//	current_prob1+=log(test_dist/(2*test_A*(test_iso.u-test_iso.u_i)+test_iso.v-test_iso.v_i));
	//	current_prob1+=log(test_iso.Jac);

//current_prob1+=in_sample((test_iso.u*pow(test_A,2)+test_iso.v*test_A+test_iso.w)+test_dist_mod+test_iso.r0, (test_iso.u_i*pow(test_A,2)+test_iso.v_i*test_A+test_iso.w_i)+test_dist_mod+test_iso.i0, (test_iso.u_ha*pow(test_A,2)+test_iso.v_ha*test_A+test_iso.w_ha)+test_dist_mod+test_iso.ha0, d_r, d_i, d_ha);

		if (floor(test_dist/100)<A_mean.size())
		{
			A_prob=0;
	//		current_prob1+=-log(A_mean[floor(test_dist/100)][3]*test_A) - pow(log(test_A)-A_mean[floor(test_dist/100)][2],2)/(2*pow(A_mean[floor(test_dist/100)][3],2));

	// Correction to prior to account for incompletness due to mag limits

		//	A_prob=-log(in_sample2(A_max, A_mean[floor(test_dist/100)][2],A_mean[floor(test_dist/100)][3]));
	
		/*	if (A_min>0){A_prob=-log(cdf_normal_fast(log(A_max),A_mean[floor(test_dist/100)][2],A_mean[floor(test_dist/100)][3])-cdf_normal_fast(log(A_min),A_mean[floor(test_dist/100)][2],A_mean[floor(test_dist/100)][3]));}
			else {A_prob=-log(cdf_normal_fast(log(A_max),A_mean[floor(test_dist/100)][2],A_mean[floor(test_dist/100)][3]));}
			if (isinf(A_prob))
			{
				A_prob=-log(cdf_normal_smallx(log(A_max),A_mean[floor(test_dist/100)][2],A_mean[floor(test_dist/100)][3]));
				//cout << A_max << " " << A_min << " " << cdf_normal_smallx(log(A_max),A_mean[floor(test_dist/100)][2],A_mean[floor(test_dist/100)][3]) << " " << A_mean[floor(test_dist/100)][0] << endl;
			}*/
	
		//	current_prob1+=A_prob;


		}
		else
		{
			A_prob=0;
	//		current_prob1+=-log(A_mean[A_mean.size()-1][3]*test_A) - pow(log(test_A)-A_mean[A_mean.size()-1][2],2)/(2*pow(A_mean[A_mean.size()-1][3],2));

	// Correction to prior to account for incompletness due to mag limits

		//	A_prob=-log(in_sample2(A_max, A_mean[A_mean.size()-1][2],A_mean[A_mean.size()-1][3]));
	
		/*	if (A_min>0){A_prob=-log(cdf_normal_fast(log(A_max),A_mean[A_mean.size()-1][2],A_mean[A_mean.size()-1][3])-cdf_normal_fast(log(A_min),A_mean[A_mean.size()-1][2],A_mean[A_mean.size()-1][3]));}
			else {A_prob=-log(cdf_normal_fast(log(A_max),A_mean[A_mean.size()-1][2],A_mean[A_mean.size()-1][3]));}
			if (isinf(A_prob))
			{
				A_prob=-log(cdf_normal_smallx(log(A_max),A_mean[A_mean.size()-1][2],A_mean[A_mean.size()-1][3]));
				//cout << A_max << " " << A_min << " " << cdf_normal_smallx(log(A_max),A_mean[floor(test_dist/100)][2],A_mean[floor(test_dist/100)][3]) << " " << A_mean[floor(test_dist/100)][0] << endl;
			}*/

		//	current_prob1+=A_prob;*/
		}


	return current_prob1;
}

void iphas_obj::initial_guess(vector<iso_obj> &isochrones, vector<iso_obj> &guess_set, vector<vector <double> > &A_mean)
{
	for (int it=0; it<guess_set.size();it++)
	{
		if (r-ha>guess_set[it].redline(r-i))
		{
			last_logT=(((r-ha)-guess_set[it].redline(r-i))*guess_set[it-1].logT + (guess_set[it-1].redline(r-i)-(r-ha))*guess_set[it].logT)/(guess_set[it-1].redline(r-i)-guess_set[it].redline(r-i));
			last_logg=(((r-ha)-guess_set[it].redline(r-i))*guess_set[it-1].logg + (guess_set[it-1].redline(r-i)-(r-ha))*guess_set[it].logg)/(guess_set[it-1].redline(r-i)-guess_set[it].redline(r-i));
			last_iso=iso_get_Tg(0.,last_logT, last_logg, isochrones);
			break;
		}
	}
	last_A=quadratic(last_iso.u-last_iso.u_i, last_iso.v-last_iso.v_i, (last_iso.w-last_iso.w_i)+(last_iso.r0-last_iso.i0)-(r-i), +1);
	last_dist_mod=r-(last_iso.u*pow(last_A,2)+last_iso.v*last_A+last_iso.w)-last_iso.r0;
//	last_dist=pow(10,last_dist_mod/5+1);

//	last_iso=iso_get(0.,real_Mi, 7.08, isochrones);
//	last_A=real_A;
//	last_dist_mod=5*log10(real_dist/10);

	last_rmag=r;
	last_ri=(last_iso.r0-last_iso.i0)+(last_iso.u-last_iso.u_i)*pow(last_A,2) + (last_iso.v-last_iso.v_i)*last_A + (last_iso.w-last_iso.w_i);

	if (last_A<0){last_A=0.02;}

// Find prob

	last_prob=prob_eval(last_iso, last_A, last_dist_mod, A_mean);

	best_prob=last_prob;
	best_iso=last_iso;
	best_A=last_A;
	best_dist_mod=last_dist_mod;
	best_it=0;
	
}




void iphas_obj::star_try1(vector<iso_obj> &isochrones, double &l, double &b, vector<vector <double> > &A_mean)
{
	iso_obj test_iso;
	double test_dist_mod, test_A, test_dist;
	double test_feh, test_logT, test_logg;
	double test_rmag, test_ri;

	double current_prob, transition_prob=0;
	double sigma2_LN, mu_LN;

//-----------------------------------------------------------------------------------------------------------
// First isochrone posn.

	test_feh=last_iso.feh;//+gsl_ran_gaussian_ziggurat(rng_handle,feh_sd);//Z.Next()*feh_sd;
	test_logT=last_iso.logT+gsl_ran_gaussian_ziggurat(rng_handle,logT_sd);//Z.Next()*logT_sd;
	test_logg=last_iso.logg+gsl_ran_gaussian_ziggurat(rng_handle,logg_sd);//Z.Next()*logg_sd;
	try {test_iso=iso_get_Tg(test_feh, test_logT, test_logg, isochrones);}
	catch (int e){no_accept++; return;}
//	test_A=VLN.Next(last_A, A_sd);//A_chain.back()+Z.Next()*A_sd;//
//	test_dist_mod=last_dist_mod + Z.Next()*dist_mod_sd;

	test_rmag=last_rmag+gsl_ran_gaussian_ziggurat(rng_handle,rmag_sd);//Z.Next()*d_r/2;
	test_ri=last_ri+gsl_ran_gaussian_ziggurat(rng_handle,ri_sd);//Z.Next()*(d_r*d_r+d_i*d_i)/2;

	test_A=quadratic(test_iso.u-test_iso.u_i, test_iso.v-test_iso.v_i, (test_iso.w-test_iso.w_i)+(test_iso.r0-test_iso.i0)-test_ri, +1);
	test_dist_mod=test_rmag-(test_iso.u*pow(test_A,2)+test_iso.v*test_A+test_iso.w)-test_iso.r0;


// Find probability of this parameter set -isochrone
//	current_prob=prob_eval(test_iso, last_A, last_dist_mod, A_mean);
//	current_prob=prob_eval(test_iso, test_A, last_dist_mod, A_mean);
	current_prob=prob_eval(test_iso, test_A, test_dist_mod, A_mean);

// Metropolis-Hastings algorithm step

	transition_prob=0;

	// From new to old
	// Mi
//	sigma2_LN=log(1+pow(Mi_sd/test_Mi,2));
//	mu_LN=log(test_Mi)-sigma2_LN/2;		
//	transition_prob+=-log(sigma2_LN)/2-log(last_iso.Mi) - pow(log(last_iso.Mi)-mu_LN,2)/(2*sigma2_LN);
	// A
//	sigma2_LN=log(1+pow(A_sd/test_A,2));
//	mu_LN=log(test_A)-sigma2_LN/2;		
//	transition_prob+=-log(sigma2_LN)/2-log(last_A) - pow(log(last_A)-mu_LN,2)/(2*sigma2_LN);

	// From old to new
	// Mi
//	sigma2_LN=log(1+pow(Mi_sd/last_iso.Mi,2));
//	mu_LN=log(last_iso.Mi)-sigma2_LN/2;		
//	transition_prob-=-log(sigma2_LN)/2-log(test_Mi) - pow(log(test_Mi)-mu_LN,2)/(2*sigma2_LN);
	// A
//	sigma2_LN=log(1+pow(A_sd/last_A,2));
//	mu_LN=log(last_A)-sigma2_LN/2;		
//	transition_prob-=-log(sigma2_LN)/2-log(test_A) - pow(log(test_A)-mu_LN,2)/(2*sigma2_LN);


	if (current_prob-last_prob+transition_prob>0 )		// New parameter set better => Accept
	{
		last_iso=test_iso;
		last_A=test_A;
		last_dist_mod=test_dist_mod;

		last_rmag=test_rmag;
		last_ri=test_ri;

		last_prob=current_prob;
		no_accept++;
		no_accept2++;

		/*if (no_accept/400.==floor(no_accept/400.))
		{
			iso_obj_chain.push_back(last_iso);
		//	dist_mod_chain.push_back(last_dist_mod);
			dist_mod_chain.push_back(test_dist_mod);
		//	A_chain.push_back(last_A);
			A_chain.push_back(test_A);
		}*/
	
		//without_change=0;

		if (current_prob>best_prob)
		{
			best_prob=current_prob;
			best_iso=test_iso;
		//	best_A=last_A;
			best_A=test_A;
		//	best_dist_mod=last_dist_mod;
			best_dist_mod=test_dist_mod;
			best_it=A_chain.size();
		}
	}

	else if (exp(current_prob-last_prob+transition_prob)>U.Next() )	// New set worse => accept with P=P(new)/P(old)
	{
		last_iso=test_iso;
		last_A=test_A;
		last_dist_mod=test_dist_mod;

		last_rmag=test_rmag;
		last_ri=test_ri;

		last_prob=current_prob;

		no_accept++;
		no_accept2++;

	/*	if (no_accept/400.==floor(no_accept/400.))
		{
			iso_obj_chain.push_back(last_iso);
		//	dist_mod_chain.push_back(last_dist_mod);
			dist_mod_chain.push_back(test_dist_mod);
		//	A_chain.push_back(last_A);
			A_chain.push_back(test_A);
		}*/

	}
	else 
	{
		//without_change++;
		/* if (without_change>100){cout << "fail: " << without_change << " "  << current_prob << " " << previous_prob << " " << transition_prob << endl;}*/
		no_accept++;
	}

	if (no_accept/100.==floor(no_accept/100.))
	{
		iso_obj_chain.push_back(last_iso);
		dist_mod_chain.push_back(last_dist_mod);
		A_chain.push_back(last_A);
		prob_chain.push_back(last_prob);
		A_prob_chain.push_back(get_A_prob(last_iso, last_A, last_dist_mod, A_mean));

		rx_chain.push_back((last_iso.u*pow(last_A,2)+last_iso.v*last_A+last_iso.w)+last_dist_mod+last_iso.r0);
		ix_chain.push_back((last_iso.u_i*pow(last_A,2)+last_iso.v_i*last_A+last_iso.w_i)+last_dist_mod+last_iso.i0);
		hax_chain.push_back((last_iso.u_ha*pow(last_A,2)+last_iso.v_ha*last_A+last_iso.w_ha)+last_dist_mod+last_iso.ha0);
	}


}

void iphas_obj::mean_intervals(void)
{
// measure estimated means and ~credible intervals for each param. - then output

//	cout << "out: " << dist_mod_1/A_chain.size() << " " << A_1/A_chain.size() << " " << feh_1/A_chain.size() << " " << Mi_1/A_chain.size() << " " << A_chain.size() << endl;

	double A_sum_in=0, A_sum2_in=0, d_sum_in=0, d_sum2_in=0, r_i0_sum_in=0, r_i0_sum2_in=0;
	double Mi_sum=0, Mi_sum2=0, logAge_sum=0, logAge_sum2=0, feh_sum=0, feh_sum2=0;
	double logT_sum=0, logT_sum2=0, logg_sum=0, logg_sum2=0;
	double prob_sum=0, A_prob_sum=0;
	double rx_sum=0, ix_sum=0, hax_sum=0;
	for (int n=floor(0.5*A_chain.size()); n<A_chain.size(); n++)
	{
		A_sum_in+=A_chain[n];
		A_sum2_in+=pow(A_chain[n],2);
		d_sum_in+=pow(10,dist_mod_chain[n]/5 +1);
		d_sum2_in+=pow(pow(10,dist_mod_chain[n]/5 +1),2);
		r_i0_sum_in+=iso_obj_chain[n].r0 - iso_obj_chain[n].i0;
		r_i0_sum2_in+=pow(iso_obj_chain[n].r0 - iso_obj_chain[n].i0,2);

		Mi_sum+=iso_obj_chain[n].Mi;
		Mi_sum2+=pow(iso_obj_chain[n].Mi,2);
		logAge_sum+=iso_obj_chain[n].logAge;
		logAge_sum2+=pow(iso_obj_chain[n].logAge,2);
		feh_sum+=iso_obj_chain[n].feh;
		feh_sum2+=pow(iso_obj_chain[n].feh,2);

		logT_sum+=iso_obj_chain[n].logT;
		logT_sum2+=pow(iso_obj_chain[n].logT,2);
		logg_sum+=iso_obj_chain[n].logg;
		logg_sum2+=pow(iso_obj_chain[n].logg,2);

		prob_sum+=prob_chain[n];
		A_prob_sum+=A_prob_chain[A_prob_chain.size()-1];

		rx_sum+=rx_chain[A_prob_chain.size()-1];
		ix_sum+=ix_chain[A_prob_chain.size()-1];
		hax_sum+=hax_chain[A_prob_chain.size()-1];
	}

//	A =best_A;//A_chain.back();
	A= A_sum_in/ceil(0.5*A_chain.size());
//	dist = pow(10,best_dist_mod/5+1);// d_sum_in/ceil(0.5*dist_mod_chain.size());
	dist = d_sum_in/ceil(0.5*dist_mod_chain.size());
//	r_i0 =best_iso.r0-best_iso.i0;//iso_obj_chain.back().r0 - iso_obj_chain.back().i0;//
	r_i0 =  r_i0_sum_in/ceil(0.5*iso_obj_chain.size());

	Mi=Mi_sum/ceil(0.5*iso_obj_chain.size());
	logAge=logAge_sum/ceil(0.5*iso_obj_chain.size());
	feh=feh_sum/ceil(0.5*iso_obj_chain.size());

	logT=last_iso.logT;//logT_sum/ceil(0.5*iso_obj_chain.size());
	logg=last_iso.logg;//logg_sum/ceil(0.5*iso_obj_chain.size());

	d_A=sqrt(A_sum2_in/ceil(0.5*A_chain.size()) - pow(A_sum_in/ceil(0.5*A_chain.size()),2));
	d_dist=1.*no_accept2/(1.*no_accept);//(last_iso.u*pow(last_A,2)+last_iso.v*last_A+last_iso.w)+last_dist_mod+last_iso.r0;//sqrt(d_sum2_in/ceil(0.5*A_chain.size()) - pow(d_sum_in/ceil(0.5*dist_mod_chain.size()),2));
//	d_dist=best_iso.r0;
	d_r_i0 =sqrt(r_i0_sum2_in/ceil(0.5*A_chain.size()) - pow(r_i0_sum_in/ceil(0.5*A_chain.size()),2));
//	d_r_i0=best_iso.logAge-max_age(best_iso.Mi, isochrones);

	d_Mi=sqrt(Mi_sum2/ceil(0.5*iso_obj_chain.size())-pow(Mi,2));
	d_logAge=sqrt(logAge_sum2/ceil(0.5*iso_obj_chain.size())-pow(logAge,2));
	d_feh=sqrt(feh_sum2/ceil(0.5*iso_obj_chain.size())-pow(feh,2));
	
	distbin=floor(dist/100);

	mean_prob=prob_sum/ceil(0.5*prob_chain.size());
	mean_A_prob=A_prob_sum/ceil(0.5*A_prob_chain.size());

	rx=rx_sum/ceil(0.5*rx_chain.size());
	ix=ix_sum/ceil(0.5*ix_chain.size());
	hax=hax_sum/ceil(0.5*hax_chain.size());
}


