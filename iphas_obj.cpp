#include "iphas_obj.h"


iphas_obj::iphas_obj(float r_input, float i_input, float ha_input, float d_r_input,float d_i_input, float d_ha_input, float l_input, float b_input)
{

	r=r_input;	
	i=i_input;
	ha=ha_input;
	J=-99;
	H=-99;
	K=-99;
	d_r=d_r_input;
	d_i=d_i_input;
	d_ha=d_ha_input;
	l=l_input * PI/180.;
	b=b_input* PI/180.;

	feh_sd=0.008;
	logT_sd=0.003;
	logg_sd=0.008;
	A_sd=sqrt(d_r*d_r+d_i*d_i)/4;
	dist_mod_sd=0.5*d_r;

	ri_sd=min(sqrt(d_r*d_r+d_i*d_i), 0.04f);
	rmag_sd=d_r;

	no_accept=0;
	cluster_weight=1;

}

iphas_obj::iphas_obj(float r_input, float i_input, float ha_input, float d_r_input,float d_i_input, float d_ha_input, float l_input, float b_input, float real_dist_in, float real_A_in, float real_Mi_in, float real_logAge_in, float real_feh_in)
{

	r=r_input;	
	i=i_input;
	ha=ha_input;
	J=-99;
	H=-99;
	K=-99;
	d_r=d_r_input;
	d_i=d_i_input;
	d_ha=d_ha_input;
	l=180* PI/180.;//l_input;
	b=0* PI/180.;//b_input;

	feh_sd=0.008;
	logT_sd=0.003*4.;
	logg_sd=0.008*4.;
	A_sd=sqrt(d_r*d_r+d_i*d_i)/4;
	dist_mod_sd=0.5*d_r;

	ri_sd=min(sqrt(d_r*d_r+d_i*d_i), 0.04f);
	rmag_sd=d_r;

	no_accept=0;
	no_accept2=0;

	real_dist=real_dist_in;
	real_A=real_A_in;
	real_Mi=real_Mi_in;
	real_logAge=real_logAge_in;
	real_feh=real_feh_in;

	cluster_weight=1;
}

iphas_obj::iphas_obj(float d_r_input, float d_i_input, float d_ha_input)
{
	d_r=d_r_input;
	d_i=d_i_input;
	d_ha=d_ha_input;
	cluster_weight=1;
}

iphas_obj::iphas_obj(float P1_input, float P2_input, float P3_input, float d_P1_input, float d_P2_input, float d_P3_input, float l_input, float b_input, string source)
{
	if (source=="2MASS")
	{
		r=-99;	
		i=-99;
		ha=-99;
		J=P1_input;
		H=P2_input;
		K=P3_input;
		d_J=d_P1_input;
		d_H=d_P2_input;
		d_K=d_P3_input;
		l=l_input * PI/180.;
		b=b_input* PI/180.;
	
		feh_sd=0.008;
		logT_sd=0.003;
		logg_sd=0.008;
		A_sd=sqrt(d_J*d_J+d_H*d_H)/4;
		dist_mod_sd=0.5*d_J;

		ri_sd=min(sqrt(d_r*d_r+d_i*d_i), 0.04f);
		rmag_sd=d_r;
	
		no_accept=0;
		cluster_weight=1;
	}
}



void iphas_obj::set_mag_weight(float mag_weight_input)
{
   mag_weight=mag_weight_input;
}

float iphas_obj::likelihood_eval(iso_obj test_iso, float test_A, float test_dist_mod)
{	
// Find probability of this parameter set - distance
	//cout << test_iso.Mi << " " << test_iso.logAge << " " << test_A << " " << test_dist_mod << " " ;
	// Find P(S|y,x,sigma_y)
		float current_prob1=0;

		if (r>-98){current_prob1+=-log(exp((r-r_max)*4)+1) ;}
		if (i>-98){current_prob1+=-log(exp((i-i_max)*4)+1) ;}
		if (ha>-98){current_prob1+=-log(exp((ha-ha_max)*4)+1) ;}

		if (J>-98){current_prob1+=-log(exp((J-J_max)*4)+1) ;}
		if (H>-98){current_prob1+=-log(exp((H-H_max)*4)+1) ;}
		if (K>-98){current_prob1+=-log(exp((K-K_max)*4)+1) ;}

	// Find p(y|x,sigma_y) 

		if (r>-98){current_prob1+=	-pow(pow(10,r-(test_iso.u*pow(test_A,2)+test_iso.v*test_A+test_iso.w)-test_dist_mod-test_iso.r0)-1,2)/(5*d_r*d_r) ;}
		if (i>-98){current_prob1+=	-pow(pow(10,i-(test_iso.u_i*pow(test_A,2)+test_iso.v_i*test_A+test_iso.w_i)-test_dist_mod-test_iso.i0)-1,2)/(5*d_i*d_i) ;}
		if (ha>-98){current_prob1+=	-pow(pow(10,ha-(test_iso.u_ha*pow(test_A,2)+test_iso.v_ha*test_A+test_iso.w_ha)-test_dist_mod-test_iso.ha0)-1,2)/(5*d_ha*d_ha);}

		if (J>-98){current_prob1+=	-pow(J-(test_iso.u_J*pow(test_A,2)+test_iso.v_J*test_A+test_iso.w_J)-test_dist_mod-test_iso.J0,2)/(2*d_J*d_K) ;}
		if (H>-98){current_prob1+=	-pow(H-(test_iso.u_H*pow(test_A,2)+test_iso.v_H*test_A+test_iso.w_H)-test_dist_mod-test_iso.H0,2)/(2*d_H*d_H) ;}
		if (K>-98){current_prob1+=	-pow(K-(test_iso.u_K*pow(test_A,2)+test_iso.v_K*test_A+test_iso.w_K)-test_dist_mod-test_iso.K0,2)/(2*d_K*d_K) ;}


		// IMF - Scalo type?
		current_prob1+=log(test_iso.IMF());
		current_prob1-=test_iso.logAge;

		current_prob1+=log_prior(test_dist_mod, test_iso.feh, l, b);

		current_prob1+=test_dist_mod*LN_TEN/5 - log(2*test_A*(test_iso.u-test_iso.u_i)+test_iso.v-test_iso.v_i);
		current_prob1+=log(test_iso.Jac);
//*/


	return current_prob1*cluster_weight;
} 


void iphas_obj::initial_guess(vector<iso_obj> &isochrones, vector<iso_obj> &guess_set, vector<bin_obj> &bin_mean)
{
	//IPHAS version
/*	for (int it=0; it<guess_set.size();it++)
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
*/
// 2MASS version
	int best_it=0, best_A=0, best_dist_mod=0;
	float best_diff=-9999;
	int prob_sum=0;
	float prob1;
	
	int index;
	for (int it=0; it<guess_set.size();it++)
	{
		index=pow(10, (K - ((J-K)-(guess_set[it].J0-guess_set[it].K0))*0.68 - guess_set[it].K0)/5.+1)/100.;
		if (int(index) < 150)
		{
//			if (abs(((J-K)-(guess_set[it].J0-guess_set[it].K0))/0.164 - 
//			A_mean[int(index)][0] )<best_diff)
//			{
//				best_diff=abs(((J-K)-(guess_set[it].J0-guess_set[it].K0))/0.164 - A_mean[int(index)][0]);
//				best_it=it;
//			cout << it  << " " << K << " " << J-K << " " << ((J-K)-(guess_set[it].J0-guess_set[it].K0))/0.164 << " " << pow(10, (K - ((J-K)-(guess_set[it].J0-guess_set[it].K0))*0.68 - guess_set[it].K0)/5.+1) 
//			<< " " << (((J-K)-(guess_set[it].J0-guess_set[it].K0))/0.164 - A_mean[int(index)][0]) << " new" << endl;
//			}
//			else
//			{
//				cout << it  << " " << K << " " << J-K << " " << ((J-K)-(guess_set[it].J0-guess_set[it].K0))/0.164 << " " << pow(10, (K - ((J-K)-(guess_set[it].J0-guess_set[it].K0))*0.68 - guess_set[it].K0)/5.+1) 
//				<< " " << (((J-K)-(guess_set[it].J0-guess_set[it].K0))/0.164 - A_mean[int(index)][0]) << " " << best_diff << endl;
//			}
			prob1=-pow((H-K)-(guess_set[it].H0-guess_set[it].K0) - ((J-K)-(guess_set[it].J0-guess_set[it].K0))*.390,2)/(d_H*d_H+d_K*d_K)
				-pow(( ((J-K)-(guess_set[it].J0-guess_set[it].K0))/0.164 - bin_mean[int(index)].last_mean_A )/(bin_mean[int(index)].last_sd_A),2 );
			if (prob1>best_diff)
			{
				best_diff=prob1;
				best_it=it;
			}
			//	cout << it << " " << pow((H-K)-(guess_set[it].H0-guess_set[it].K0) - ((J-K)-(guess_set[it].J0-guess_set[it].K0))*.390,2)/(d_H*d_H+d_K*d_K) << " " << prob1 << " " << best_it  << " " << best_diff<< endl;
		}

//		for (float A_test=((J-K)-(guess_set[it].J0-guess_set[it].K0))/0.164 - 3*(dJ+dK)/0.164;
//			A_test<((J-K)-(guess_set[it].J0-guess_set[it].K0))/0.164 + 3*(dJ+dK)/0.164; A_test+=(dJ+dK)/(3*0.164) )
//		{
//			for (float mod_test=K-A_test*0.112-guess_set[it].K0-dK*3; mod_test<mod_test=K-A_test*0.112-guess_set[it].K0+dK*3
//				mod_test+=dK/3.)
//			{

	}

	last_logT=guess_set[best_it].logT;	
	last_logg=guess_set[best_it].logg;
	last_iso=guess_set[best_it];

	last_A=((J-K)-(last_iso.J0-last_iso.K0))/0.164;
	last_dist_mod = K - last_A*0.112 -last_iso.K0;
	last_dist=pow(10,last_dist_mod/5+1);

	if (last_A<0){last_A=0.02;}

	last_bin = &bin_mean[floor(last_dist/100.)];
	last_bin->initial_add(this);

	//cout << last_dist << " " << last_A << " " << last_dist_mod << " " << best_it << endl;



	last_rmag=r;
	last_ri=(last_iso.r0-last_iso.i0)+(last_iso.u-last_iso.u_i)*pow(last_A,2) + (last_iso.v-last_iso.v_i)*last_A + (last_iso.w-last_iso.w_i);



// Find prob

	last_prob=likelihood_eval(last_iso, last_A, last_dist_mod);

	best_prob=last_prob+last_A_prob;
	best_iso=last_iso;
	best_A=last_A;
	best_dist_mod=last_dist_mod;
	best_it=0;
	
}




void iphas_obj::star_try1(vector<iso_obj> &isochrones, float &l, float &b, vector<bin_obj> &bin_mean )
{
	iso_obj test_iso;

	float test_rmag, test_ri;
	bin_obj* test_bin;

	float current_prob, current_A_prob, transition_prob=0;
	float sigma2_LN, mu_LN;

//-----------------------------------------------------------------------------------------------------------
// First isochrone posn.

	test_feh=last_iso.feh;//+gsl_ran_gaussian_ziggurat(rng_handle,feh_sd);//Z.Next()*feh_sd;
	if (gsl_ran_flat(rng_handle,0.,1.)<0.95)
	{
		test_logT=last_iso.logT+gsl_ran_gaussian_ziggurat(rng_handle,logT_sd);//Z.Next()*logT_sd;
		test_logg=last_iso.logg+gsl_ran_gaussian_ziggurat(rng_handle,logg_sd);//Z.Next()*logg_sd;
		try {test_iso=iso_get_Tg(test_feh, test_logT, test_logg, isochrones);}
		catch (int e){no_accept++; return;}
		test_A=last_A+gsl_ran_gaussian_ziggurat(rng_handle,A_sd);//VLN.Next(last_A, A_sd);//A_chain.back()+Z.Next()*A_sd;//
		test_dist_mod=last_dist_mod+gsl_ran_gaussian_ziggurat(rng_handle,dist_mod_sd);//+ Z.Next()*dist_mod_sd;
	}
	else
	{
		test_logT=last_iso.logT+gsl_ran_gaussian_ziggurat(rng_handle,5*logT_sd);//Z.Next()*logT_sd;
		test_logg=last_iso.logg+gsl_ran_gaussian_ziggurat(rng_handle,5*logg_sd);//Z.Next()*logg_sd;
		try {test_iso=iso_get_Tg(test_feh, test_logT, test_logg, isochrones);}
		catch (int e){no_accept++; return;}
		test_A=last_A+gsl_ran_gaussian_ziggurat(rng_handle,5*A_sd);//VLN.Next(last_A, A_sd);//A_chain.back()+Z.Next()*A_sd;//
		test_dist_mod=last_dist_mod+gsl_ran_gaussian_ziggurat(rng_handle,5*dist_mod_sd);//+ Z.Next()*dist_mod_sd;
	}
	while (test_dist_mod>15.88){test_dist_mod=last_dist_mod+gsl_ran_gaussian_ziggurat(rng_handle,dist_mod_sd);}
	test_dist=pow(10, test_dist_mod/5.+1);
	test_bin=&bin_mean[floor(test_dist/100.)];

	if (test_bin!=last_bin)
	{
		test_bin->try_add(this);
		last_bin->try_remove(this);
	}
	else { last_bin->try_update(this);}


	current_prob=likelihood_eval(test_iso, test_A, test_dist_mod);

	if (test_bin!=last_bin)
	{ 
		current_A_prob=test_bin->test_prob + last_bin->test_prob; 
		last_A_prob=test_bin->last_prob + last_bin->last_prob; 
	}
	else 
	{
		current_A_prob=test_bin->test_prob; 
		last_A_prob=test_bin->last_prob; 

	}

// Metropolis-Hastings algorithm step

	transition_prob=0;



	if (current_prob+current_A_prob-(last_prob+last_A_prob)+transition_prob>0 )		// New parameter set better => Accept
	{
		last_iso=test_iso;
		last_A=test_A;
		last_dist_mod=test_dist_mod;
		if (test_bin!=last_bin)
		{
			test_bin->accept();
			last_bin->accept();
		}
		else {last_bin->accept();}
		last_bin=test_bin;

//		last_rmag=test_rmag;
//		last_ri=test_ri;

		last_prob=current_prob;
		last_A_prob=current_A_prob;
		no_accept++;
		no_accept2++;



		if (current_prob+current_A_prob>best_prob)
		{
			best_prob=current_prob+current_A_prob;
			best_iso=test_iso;
			best_A=test_A;
			best_dist_mod=test_dist_mod;
			best_it=A_chain.size();
		}
	}

	else if (exp(current_prob+current_A_prob-(last_prob+last_A_prob)+transition_prob)>gsl_ran_flat(rng_handle, 0, 1) )	// New set worse => accept with P=P(new)/P(old)
	{
		last_iso=test_iso;
		last_A=test_A;
		last_dist_mod=test_dist_mod;
		if (test_bin!=last_bin)
		{
			test_bin->accept();
			last_bin->accept();
		}
		else {last_bin->accept();}
		last_bin=test_bin;

//		last_rmag=test_rmag;
//		last_ri=test_ri;

		last_prob=current_prob;
		last_A_prob=current_A_prob;

		no_accept++;
		no_accept2++;


	}
	else 
	{
		no_accept++;
		if (test_bin!=last_bin)
		{
			test_bin->reject();
			last_bin->reject();
		}
		else {last_bin->reject();}	
	}

	if (no_accept/100.==floor(no_accept/100.))
	{
		iso_obj_chain.push_back(last_iso);
		dist_mod_chain.push_back(last_dist_mod);
		A_chain.push_back(last_A);
		prob_chain.push_back(last_prob);

		rx_chain.push_back((last_iso.u*pow(last_A,2)+last_iso.v*last_A+last_iso.w)+last_dist_mod+last_iso.r0);
		ix_chain.push_back((last_iso.u_i*pow(last_A,2)+last_iso.v_i*last_A+last_iso.w_i)+last_dist_mod+last_iso.i0);
		hax_chain.push_back((last_iso.u_ha*pow(last_A,2)+last_iso.v_ha*last_A+last_iso.w_ha)+last_dist_mod+last_iso.ha0);
	}


}

void iphas_obj::mean_intervals(void)
{
// measure estimated means and ~credible intervals for each param. - then output

	float A_sum_in=0, A_sum2_in=0, d_sum_in=0, d_sum2_in=0, r_i0_sum_in=0, r_i0_sum2_in=0;
	float Mi_sum=0, Mi_sum2=0, logAge_sum=0, logAge_sum2=0, feh_sum=0, feh_sum2=0;
	float logT_sum=0, logT_sum2=0, logg_sum=0, logg_sum2=0;
	float prob_sum=0, A_prob_sum=0;
	float rx_sum=0, ix_sum=0, hax_sum=0;
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
	//	A_prob_sum+=A_prob_chain[A_prob_chain.size()-1];

		rx_sum+=rx_chain[rx_chain.size()-1];
		ix_sum+=ix_chain[rx_chain.size()-1];
		hax_sum+=hax_chain[rx_chain.size()-1];
	}

	A= A_sum_in/ceil(0.5*A_chain.size());
	dist = d_sum_in/ceil(0.5*dist_mod_chain.size());
	r_i0 =  r_i0_sum_in/ceil(0.5*iso_obj_chain.size());

	Mi=Mi_sum/ceil(0.5*iso_obj_chain.size());
	logAge=logAge_sum/ceil(0.5*iso_obj_chain.size());
	feh=feh_sum/ceil(0.5*iso_obj_chain.size());

	logT=logT_sum/ceil(0.5*iso_obj_chain.size());
	logg=logg_sum/ceil(0.5*iso_obj_chain.size());

	d_A=sqrt(A_sum2_in/ceil(0.5*A_chain.size()) - pow(A_sum_in/ceil(0.5*A_chain.size()),2));
	d_dist=1.*no_accept2/(1.*no_accept);//(last_iso.u*pow(last_A,2)+last_iso.v*last_A+last_iso.w)+last_dist_mod+last_iso.r0;//sqrt(d_sum2_in/ceil(0.5*A_chain.size()) - pow(d_sum_in/ceil(0.5*dist_mod_chain.size()),2));
	d_r_i0 =sqrt(r_i0_sum2_in/ceil(0.5*A_chain.size()) - pow(r_i0_sum_in/ceil(0.5*A_chain.size()),2));

	d_Mi=sqrt(Mi_sum2/ceil(0.5*iso_obj_chain.size())-pow(Mi,2));
	d_logAge=sqrt(logAge_sum2/ceil(0.5*iso_obj_chain.size())-pow(logAge,2));
	d_feh=sqrt(feh_sum2/ceil(0.5*iso_obj_chain.size())-pow(feh,2));
	
	distbin=floor(dist/100);

	mean_prob=prob_sum/ceil(0.5*prob_chain.size());
//	mean_A_prob=A_prob_sum/ceil(0.5*A_prob_chain.size());

	rx=rx_sum/ceil(0.5*rx_chain.size());
	ix=ix_sum/ceil(0.5*ix_chain.size());
	hax=hax_sum/ceil(0.5*hax_chain.size());
}


vector <float> iphas_obj::acl_calc(void)
{

	vector <float> acf (A_chain.size()-1, 0.);

	for (int it1=1; it1<A_chain.size(); it1++)
	{
		for (int lag=1; lag<it1; lag++)
		{
			acf[lag-1]+=(A_chain[it1]-A)*(A_chain[it1-lag]-A);
		}

	}

	for (int it1=0; it1<acf.size(); it1++)
	{
		acf[it1]/=(acf.size()-it1)*d_A;
	}

	return acf;
}

float log_prior(float test_dist_mod, float test_feh, float l, float b)
{
	float current_prob=0;
	float test_dist=pow(10,test_dist_mod/5+1);

	float /*cosb=1, sinb=0, cosl=-1,*/ R_gal;
	R_gal=sqrt(test_dist*test_dist*cos(b)*cos(b) + 64000000-16000*test_dist*cos(l)*cos(b));

	// density profile
	//	current_prob1+=-R_gal/2500 -test_dist*sinb/200;
	if (R_gal<13000){current_prob+=-R_gal/3000 -test_dist*sin(b)/200;}	// change back to 13000 for real data
	else {current_prob+=-R_gal/1200 +6.5 -test_dist*sin(b)/200;} 		// -6.5=13000/3000 - 13000/1200

	// metallicity profile
	//current_prob+=-pow(test_feh-(R_gal-8000.)*0.00007,2)/(2*0.5);
	//	current_prob1+=-pow(test_iso.feh,2)/(2*pow(0.05,2));

	//	// dist^2 term 
	current_prob+=2*log(test_dist);

	return current_prob;
}

