#include "sl_obj.h"

sl_obj::sl_obj(void)
{
	Cov_Mat.reserve(150);


	// set default MIN vals
	r_min=13.5; 
	i_min=12.0;
	ha_min=12.0;
	J_min=8.;
	H_min=7.;
	K_min=6.;
	// set default MAX vals
	r_max=21.;
	i_max=20.;
	ha_max=20.;
	J_max=17.;
	H_max=16.;
	K_max=15.5;

	sigma_fac=0.05;
	accepted=0;
	without_change=0;
	thin=200;
	rel_length=150;

	it_num=0.;

//	neighbour_slsl.clear();
//	define_cov_mat();

	vector<bin_obj> running_A_mean(150);

	for (int it=0; it<running_A_mean.size(); it++){running_A_mean[it].dist_bin=it;}

	sl_identifier=0;
}

sl_obj::sl_obj(string filename, float l_in, float b_in, string datatype, float s_R, float s_z, int identifier)
{
	Cov_Mat.reserve(150);

	cout << "Working on " << filename << endl;

	// Set up variables ---------------------------------------------------------------------
	l=l_in;
	b=b_in;
	rootname=filename;
	rootname.erase(rootname.size()-4);

	// set default MIN vals
	r_min=13.5; 
	i_min=12.0;
	ha_min=12.0;	
	J_min=8.;
	H_min=7.;
	K_min=6.;			
	// set default MAX vals
	r_max=21.;
	i_max=20.;
	ha_max=20.;
	J_max=17.;
	H_max=16.;
	K_max=15.5;



	//	global_previous_prob=0;
	//	previous_hyperprior_prob=0;
	//	current_hyperprior_prob=0;

	// Read in data	-------------------------------------------------------------------------

	if (datatype=="iphas" || datatype=="IPHAS"){star_cat=iphas_read(filename,r_min,i_min,ha_min,r_max,i_max,ha_max);}
	else if (datatype=="2MASS" || datatype=="2mass"){star_cat=TWOMASS_read(filename,J_min, H_min, K_min,J_max,H_max,K_max);}
	else {cout << "Unrecognised datatype: " << datatype << endl;}

	// Find expected A(d) -------------------------------------------------------------------
	cout << "J:" << J_min<<" - "<<J_max << endl;
	cout << "H:" << H_min<<" - "<<H_max << endl;
	cout << "K:" << K_min<<" - "<<K_max << endl;

	//vector<bin_obj2> backup_A_mean (150);
	backup_A_mean.resize(150);
	backup_A_mean=backup_A_mean_find(l, b, s_R, s_z, true);

	sigma_fac=0.025;
	accepted=0;
	without_change=0;
	thin=200;
	rel_length=150;

	it_num=0.;
//	neighbour_slsl.clear();

//	define_cov_mat();

	running_A_mean.resize(150);

	for (int it=0; it<running_A_mean.size(); it++){running_A_mean[it].dist_bin=it;}

	sl_identifier=identifier;
}



void sl_obj::output_write(float s_R, float s_z)
{

	ofstream A_out;
	string dummy_string;
	dummy_string=rootname+".td4";
	A_out.open(dummy_string.c_str(), ios::trunc);
	//A_out << "#\tdist\tA\tsigma_A\n";

	for (int x=0; x<running_A_mean.size(); x++)
	{
		A_out 	<< x*100 + 50 << "\t" << running_A_mean[x].final_A << "\t" << running_A_mean[x].final_sd <<"\t"
			<< running_A_mean[x].final_dA<<"\t"<<running_A_mean[x].final_dsd<<"\t"<<running_A_mean[x].last_prob<<"\t"<<running_A_mean[x].last_mean_A<<"\t"<<running_A_mean[x].last_n<<" "<<running_A_mean[x].last_lnA_sum<<" " <<running_A_mean[x].last_mu<<endl;
	}
	A_out.close();


	ofstream output;
	dummy_string=rootname+"-090.dat";
	output.open(dummy_string.c_str(), ios::trunc);
	output << "#\tr\ti\tha\tJ\tH\tK\tr_i0\tdist\tA\tdistbin\td_A\td_r_i0\td_dist\td_r\td_i\td_ha\tmag_weight\tprob\tA_prob\tMi\tlogAge\tfeh\td_Mi\td_lagAge\td_feh\tlogT\tlogg\trx\tix\thax\n" ;
	for (int y=0; y<star_cat.size(); y++)
	{
		
		output << star_cat[y].r << "\t" << star_cat[y].i << "\t" << star_cat[y].ha << "\t"
			<< star_cat[y].J << "\t" << star_cat[y].H << "\t"<< star_cat[y].K << "\t"
			<< star_cat[y].r_i0 << "\t" << star_cat[y].dist << "\t" << star_cat[y].A << "\t" 
			<< star_cat[y].distbin << "\t" << star_cat[y].d_A << "\t" << star_cat[y].d_r_i0 << "\t" << star_cat[y].d_dist << "\t" << star_cat[y].d_r << "\t"
			<< star_cat[y].d_i << "\t" << star_cat[y].d_ha << "\t" << star_cat[y].last_iso.Mi  << "\t" << star_cat[y].mean_prob  << "\t" 
			<< star_cat[y].mean_A_prob  << "\t" << star_cat[y].Mi  << "\t" << star_cat[y].logAge  << "\t" << star_cat[y].feh  << "\t" << star_cat[y].d_Mi  << "\t" 
			<< star_cat[y].d_logAge  << "\t" << star_cat[y].d_feh << "\t" << star_cat[y].logT << "\t" << star_cat[y].logg << "\t" << star_cat[y].rx << "\t" 
			<< star_cat[y].ix << "\t" << star_cat[y].hax << "\n";
	}
	output.close();

	ofstream rho_out;
	dummy_string=rootname+".rho";
	rho_out.open(dummy_string.c_str(), ios::trunc);
	//A_out << "#\tdist\tA\tsigma_A\n";

	vector <float> mean_rho, mean_rel;
	mean_rho=density_find(l, b, s_R, s_z, 100);
	//mean_rho=density_find(l, b, previous_s_R, previous_s_z, 100);
	mean_rel=backup_A_mean_find(l, b, s_R, s_z, false);

	for (int x=0; x<running_A_mean.size(); x++)
	{
		rho_out << x*100 << "\t" << running_A_mean[x].final_rho << "\t" << running_A_mean[x].final_drho <<"\t" << mean_rel[x]-mean_rel[x-1] << "\t" << mean_rel[x-1] << "\n";
	}
	rho_out.close();
	
}


void sl_obj::initial_guess(vector<iso_obj> &isochrones, vector<iso_obj> &guess_set, vector <LF> &LFs, float s_R, float s_z, float A_0)
{
	cout << "sightline " << sl_identifier << " has neighbours: ";
	for (int i=0; i<neighbour_slsl.size(); i++){cout << neighbour_slsl[i]->sl_identifier << " " ;}
	cout << "and higher neighbours: ";
	for (int i=0; i<higher_neighbour_slsl.size(); i++){cout << higher_neighbour_slsl[i]->sl_identifier << " ";}
	cout << endl;	

	define_cov_mat();
//	Trace file
		
	global_previous_prob=0;
	previous_hyperprior_prob=0;
	current_hyperprior_prob=0;
	previous_xsl_prob=0;
	current_xsl_prob=0;

	proposal_sd.resize(150, vector <float> (2));
	previous_rel.resize(150, vector <float> (4));
	internal_rel.resize(150, vector <float> (2));
	previous_internal_rel.resize(150, vector <float> (2));
	hyperprior_internal_rel.resize(150, vector <float> (2));


	vector<float> backup_rho_mean;
	backup_rho_mean=backup_rho_mean_find(l, b, s_R, s_z, A_0);


	for (int i=0; i<150; i++){last_m_vec[i]=log(backup_rho_mean[i] ) - Cov_Mat.coeffRef(i,i)*pow(last_s_vec[i],2)/2. ;}
	trial_rel=mvn_gen_internal_rel_no_neighbour();
	last_ln_vec=test_ln_vec;

	for (int it=0; it<running_A_mean.size(); it++)
	{
		running_A_mean[it].last_mean_rho=trial_rel[it][0];
		running_A_mean[it].last_sd_A=trial_rel[it][1];
	}
	initial_rho_to_A();



// Dump sources not in required region of c-c

	int it_stars=0;

	// IPHAS based filter
/*	while (it_stars<star_cat.size())
	{
		if (star_cat[it_stars].r-star_cat[it_stars].ha>guess_set[0].redline(star_cat[it_stars].r-star_cat[it_stars].i) || star_cat[it_stars].r-star_cat[it_stars].ha<guess_set[guess_set.size()-1].redline(star_cat[it_stars].r-star_cat[it_stars].i))
		{
			star_cat.erase(star_cat.begin()+it_stars);
		}
		else {it_stars++;}
	}*/

	// 2MASS based filter
	while (it_stars<star_cat.size())
	{
		if (star_cat[it_stars].J<0 || star_cat[it_stars].H<0 || star_cat[it_stars].K<0)
		{
			star_cat.erase(star_cat.begin()+it_stars);
		}
		else if (star_cat[it_stars].d_J<0 || star_cat[it_stars].d_H<0 || star_cat[it_stars].d_K<0)
		{
			star_cat.erase(star_cat.begin()+it_stars);
		}
		else if (star_cat[it_stars].J-star_cat[it_stars].H < (star_cat[it_stars].H - star_cat[it_stars].K)*1.56 -.05 )
		{
			star_cat.erase(star_cat.begin()+it_stars);
		}
		else if (star_cat[it_stars].J-star_cat[it_stars].H > (star_cat[it_stars].H - star_cat[it_stars].K)*1.56 +0.6 )
		{
			star_cat.erase(star_cat.begin()+it_stars);
		}
		else {it_stars++;}
	}			

	cout << "unfiltered:" << star_cat.size() << endl;

// Make initial guess

	it_stars=0;
	while (it_stars<star_cat.size())
	{
		star_cat[it_stars].initial_guess(isochrones, guess_set, running_A_mean);
		if (star_cat[it_stars].last_A<0){star_cat[it_stars].last_A = 0.02;} 
		it_stars++;
	}

	for (int it=0; it<running_A_mean.size(); it++)
	{
		running_A_mean[it].set_last_prob(); 
		running_A_mean[it].reject(); 		
		global_previous_prob+=running_A_mean[it].last_prob;
	}

	proposed_probs.resize(star_cat.size());

	previous_norm_prob=0;
	for (int it_LF=0; it_LF<LFs.size(); it_LF++)
	{
		previous_norm_prob+=-LFs[it_LF].LF_prob_last(running_A_mean, J_min, J_max)*(star_cat.size()+1);
	}

}



void sl_obj::dist_redMCMC(vector<iso_obj> &isochrones, vector <LF> &LFs)
{

	while (it_num<150000 )
	{
		update(isochrones, LFs);
	//	if (it_num/100==floor(it_num/100.)){hyperprior_update(LFs);}
	}
}

void sl_obj::update(vector<iso_obj> &isochrones, vector <LF> &LFs)
{


		global_current_prob=0;
		global_transition_prob=0;
		global_previous_prob=0;
		current_hyperprior_prob=0;

// First vary parameters for each star

		float dummy=0;
		float dummy3=0, dummy4=0;
	//	#pragma omp parallel for  num_threads(3) //reduction(+:dummy)
		star_prob=0;
		for (int it=0; it<star_cat.size(); it++)
		{
			/*if (gsl_ran_flat(rng_handle, 0, 1)>0.){*/star_cat[it].star_try1(isochrones, l, b, running_A_mean);//};
			star_prob+=star_cat[it].last_prob;
		}

		global_previous_prob=get_A_mean_last_prob();


// Now vary hyper-parameters

		move_on=false;	
		threshold=log( gsl_ran_flat(rng_handle, 0, 1) );
		trial_rel=mvn_gen_internal_rel();

		float sss=0;//0.01;
		//theta=gsl_ran_flat(rng_handle, 0, 2*PI);
		theta=gsl_ran_gaussian_ziggurat(rng_handle, sss);
		theta_min=theta-2*PI; 
		theta_max=theta;

		while (!move_on)
		{
			sum1=0;
			for (int i=0; i<rel_length; i++)
			{
//				running_A_mean[i].test_mean_rho=((running_A_mean[i].last_mean_rho)*cos(theta) + (trial_rel[i][0])*sin(theta));
//				if (running_A_mean[i].test_mean_rho<0){continue;}

				running_A_mean[i].test_mean_rho=gsl_ran_lognormal(rng_handle,log(running_A_mean[i].last_mean_rho)-pow(0.005,2)/2,0.005);

				sum1+=running_A_mean[i].test_mean_rho;
				//running_A_mean[i].test_sd_A=exp(log(running_A_mean[i].last_sd_A)*cos(theta) + log(trial_rel[i][1])*sin(theta));
				running_A_mean[i].test_sd_A=1.75/sqrt(8.*(i+.5)) * last_s_vec[i] * sum1;//running_A_mean[i].last_sd_A;//gsl_ran_lognormal(rng_handle, log( 1.75/sqrt(8.*(i+.5)) * last_s_vec[i] * sum1) - 0.05268 , 0.3246);
				test_ln_vec[i]=log(running_A_mean[i].test_mean_rho);
			}
			rho_to_A();

	// Find probability of this parameter set

			global_current_prob=get_A_mean_test_prob();

	// Normalisation term

			current_norm_prob=0;
			for (int it_LF=0; it_LF<LFs.size(); it_LF++)
			{
				current_norm_prob+=-LFs[it_LF].LF_prob_test(running_A_mean, J_min, J_max)*(star_cat.size()+1);
			}

	// Metropolis-Hastings algorithm step

			global_transition_prob=0;

			dummy=0;
	//		#pragma omp parallel for  num_threads(3) reduction(+:dummy)
			for (int i=1; i<150; i++)
			{
			// From new to old
			// mean_A
				dummy+=log(gsl_ran_lognormal_pdf(running_A_mean[i].last_mean_rho, log(running_A_mean[i].test_mean_rho)-pow(0.005,2)/2 ,0.005));
			//	dummy+=log(gsl_ran_lognormal_pdf(previous_internal_rel[it][1], log(internal_rel[it][1])-pow(proposal_sd[it][1],2)/2 ,proposal_sd[it][1]));
			// From old to new
			// mean_A
				dummy-=log(gsl_ran_lognormal_pdf(running_A_mean[i].test_mean_rho, log(running_A_mean[i].last_mean_rho)-pow(0.005,2)/2 ,0.005));
			//	dummy-=log(gsl_ran_lognormal_pdf(internal_rel[it][1], log(previous_internal_rel[it][1])-pow(proposal_sd[it][1],2)/2 ,proposal_sd[it][1]));
			}
			global_transition_prob=dummy;

			previous_hyperprior_prob=get_rho_last_prob();
			current_hyperprior_prob=get_rho_last_prob_new_rel();
	

// Accept or reject

			accepted++;
			if (global_current_prob+current_hyperprior_prob-global_previous_prob-previous_hyperprior_prob+global_transition_prob
				 + current_norm_prob-previous_norm_prob>log(gsl_ran_flat(rng_handle, 0, 1)))//threshold)// || sss<1E-7)		// New parameter set better => Accept
			{
			//	cout << "pass " << running_A_mean[0].A_chain.size() << " " << global_current_prob << " " << running_A_mean[50].last_mean_rho << " " << running_A_mean[50].test_mean_rho << " " << previous_norm_prob << " "  << endl;
				global_previous_prob=global_current_prob;
				previous_hyperprior_prob=current_hyperprior_prob;
				previous_xsl_prob=current_xsl_prob;
				previous_norm_prob=current_norm_prob;

				for (int it=0; it<running_A_mean.size(); it++){running_A_mean[it].accept();}
				last_ln_vec=test_ln_vec;

				move_on=true;
		

			}

			else 
			{
		//	cout << "fail " << running_A_mean[0].A_chain.size() << " " << global_current_prob << " " << global_previous_prob << " " << previous_hyperprior_prob << " " << current_hyperprior_prob << " " << previous_norm_prob << " " << current_norm_prob << " " << global_transition_prob  << endl;
			if (theta>0){theta_max=theta;}
			else {theta_min=theta;}
		//	theta=gsl_ran_flat(rng_handle, theta_min, theta_max);
			sss/=5;	
			theta=gsl_ran_gaussian_ziggurat(rng_handle, sss);
			for (int it=0; it<running_A_mean.size(); it++){running_A_mean[it].reject();}
			move_on=true;


			}
		}

		if (it_num/1000.==floor(it_num/1000.)){cout << sl_identifier << " " << it_num << " " << global_previous_prob << " " << running_A_mean[50].last_mean_rho << " " << running_A_mean[rel_length-1].last_mean_A << " " << previous_hyperprior_prob << " " << previous_norm_prob << " " << accepted << " " << accepted/it_num << endl;}


		if (floor(it_num/200.)==it_num/200)
		{

			for (int it=0; it<running_A_mean.size(); it++){running_A_mean[it].chain_push_back();}
		}
		it_num++;
}

float sl_obj::get_A_mean_test_prob(void)
{
	float dummy3=0;
//	#pragma omp parallel for  num_threads(3) reduction(+:dummy3)
	for (int it=0; it<running_A_mean.size(); it++)
	{
		running_A_mean[it].set_test_prob(); 
		dummy3+=running_A_mean[it].test_prob;
	}
//	if (isnan(dummy3)){cout << running_A_mean[149].test_mu << " " << running_A_mean[149].test_sigma << " "  << running_A_mean[149].test_mean_A << " " << running_A_mean[149].test_sd_A<< " "  << running_A_mean[148].test_mean_A << " " << running_A_mean[148].test_sd_A << endl;}
	return dummy3;
}

float sl_obj::get_A_mean_last_prob(void)
{
	float dummy3=0;
//	#pragma omp parallel for  num_threads(3) reduction(+:dummy3)
	for (int it=0; it<running_A_mean.size(); it++)
	{
		running_A_mean[it].set_last_prob(); 
		dummy3+=running_A_mean[it].last_prob;
	}
	return dummy3;
}

vector < vector <float> > sl_obj::mvn_gen_internal_rel(void)
{
	vector < vector <float> > new_rel(rel_length,vector <float> (2));
	Eigen::Matrix<float, 150, 1> int_vec, temp_vec, temp_vec2;
	float sum=0.;
	for (int i=0; i<150; i++){temp_vec[i]=0;}

	for (int i=0; i<rel_length; i++){test_z_dash[i]=gsl_ran_gaussian_ziggurat(rng_handle, 1.);}

	if (higher_neighbour_slsl.size()>0.5)
	{
		for (int i=0; i<neighbour_slsl.size(); i++){temp_vec+=cond_mu_Mat*neighbour_slsl[i]->last_s_inv.asDiagonal()*
										(neighbour_slsl[i]->last_ln_vec-neighbour_slsl[i]->last_m_vec);}

		temp_vec2=last_s_vec.asDiagonal()*temp_vec;

	//	for (int i=0; i<150; i++){cout << i << " " << temp_vec2[i] << " " << neighbour_slsl[0]->last_ln_vec[i] << " " << neighbour_slsl[0]->last_m_vec[i] << " " << cond_mu_Mat.coeffRef(i,i) << endl;}
		int_vec=last_s_vec.asDiagonal()*(chol_L_cond*test_z_dash) + last_m_vec+temp_vec2;
	}
	else {int_vec=last_s_vec.asDiagonal()*(chol_L_cond*test_z_dash) + last_m_vec;}

	for (int i=0; i<rel_length; i++)
	{
		new_rel[i][0]=exp(int_vec[i]);
	//	sum+=new_rel[i][0];
	//	new_rel[i][1]=gsl_ran_lognormal(rng_handle, -0.94, 0.246);//*new_rel[i][0];
	//	new_rel[i][1]=2;//*1.75/sqrt(8.*(i+.5)) * last_s_vec[i] * sum;//gsl_ran_lognormal(rng_handle, log( 1.75/sqrt(8.*(i+.5)) * last_s_vec[i] * sum) - 0.05268 , 0.3246);
	}
	return new_rel;
}

vector < vector <float> > sl_obj::mvn_gen_internal_rel_no_neighbour(void)
{
	vector < vector <float> > new_rel(rel_length,vector <float> (2));
	Eigen::Matrix<float, 150, 1> int_vec;
	float sum=0.;

	for (int i=0; i<rel_length; i++){test_z_dash[i]=gsl_ran_gaussian_ziggurat(rng_handle, 1.);}
	int_vec=last_s_vec.asDiagonal()*(chol_L_cond*test_z_dash) + last_m_vec;

	for (int i=0; i<rel_length; i++)
	{
		new_rel[i][0]=exp(int_vec[i]);
		sum+=new_rel[i][0];
	//	new_rel[i][1]=gsl_ran_lognormal(rng_handle, -0.94, 0.246);//*new_rel[i][0];
		new_rel[i][1]=1.75/sqrt(8.*(i+.5)) * last_s_vec[i] * sum;//gsl_ran_lognormal(rng_handle, log( 1.75/sqrt(8.*(i+.5)) * last_s_vec[i] * sum) - 0.05268 , 0.3246);
		test_ln_vec[i]=int_vec[i];
	//	if (isinf(new_rel[i][0]) || isinf(new_rel[i][1])){cout << "nn" << i << " " << int_vec[i] << " " << new_rel[i][0] << " " << new_rel[i][1] << " " << sum << endl;}
	}
	return new_rel;
}

void sl_obj::mvn_gen_internal_rel_from_z_dash(void)
{
	// set up
	Eigen::Matrix<float, 150, 1> int_vec, temp_vec, temp_vec2;
	float sum=0;
	for (int i=0; i<150; i++){temp_vec[i]=0;}

	// get new rhos

	if (higher_neighbour_slsl.size()>0.5)
	{
		for (int i=0; i<higher_neighbour_slsl.size(); i++){temp_vec+=higher_neighbour_slsl[i]->last_s_inv.asDiagonal()*
									(higher_neighbour_slsl[i]->last_ln_vec-higher_neighbour_slsl[i]->test_m_vec);}

		temp_vec2=last_s_vec.asDiagonal()*(cond_mu_Mat*temp_vec);
		int_vec=/*last_s_vec.asDiagonal()*(higher_chol_L*/(last_z_dash) + test_m_vec + temp_vec2;
	}
	else {int_vec=/*last_s_vec.asDiagonal()*(higher_chol_L*/(last_z_dash) + test_m_vec ;}

	for (int i=0; i<rel_length; i++)
	{
		running_A_mean[i].test_mean_rho=exp(int_vec[i]);
		sum+=running_A_mean[i].test_mean_rho;
		running_A_mean[i].test_sd_A=1.75/sqrt(8.*(i+.5)) * last_s_vec[i] * sum;//running_A_mean[i].last_sd_A;//gsl_ran_lognormal(rng_handle, log( 1.75/sqrt(8.*(i+.5)) * last_s_vec[i] * sum) - 0.05268 , 0.3246);
		test_ln_vec[i]=int_vec[i];
	}
	rho_to_A();
//	cout << sl_identifier << " " <<  running_A_mean[149].test_mean_A << " " << int_vec[149]<< " "  << test_m_vec[149] << " " << temp_vec2[149] << " " << last_z_dash[149] << endl;
}

Eigen::Matrix<float, 150, 1> sl_obj::get_last_z_dash(void)
{
	Eigen::Matrix<float, 150, 1> temp_vec, temp_vec2, temp_vec3, x;
	for (int i=0; i<150; i++){temp_vec[i]=0;}

	if (higher_neighbour_slsl.size()>0.5)
	{
		for (int i=0; i<higher_neighbour_slsl.size(); i++){temp_vec+=higher_neighbour_slsl[i]->last_s_inv.asDiagonal()*
									(higher_neighbour_slsl[i]->last_ln_vec-higher_neighbour_slsl[i]->last_m_vec);}

		temp_vec2=last_s_vec.asDiagonal()*(cond_mu_Mat*temp_vec);
		temp_vec3=last_ln_vec-(last_m_vec+temp_vec2);
	}
	else {temp_vec3=last_ln_vec-(last_m_vec);}


	//x=  chol_L_Inv * temp_vec3;
	//return last_s_Inv.asDiagonal() * chol_L_Inv * temp_vec2;
	return  /*higher_chol_L_Inv * last_s_inv.asDiagonal() */ temp_vec3;
}


//-----------------=================================================================================

float sl_obj::get_rho_last_prob(void)
{
	Eigen::Matrix<float, 150, 1> last_s_Inv, temp_vec, temp_vec2;
	float x;
	float sum=0;
	
	for (int i=0; i<rel_length; i++)
	{
		last_s_Inv[i]=1./last_s_vec[i];
		temp_vec[i]=log(running_A_mean[i].last_mean_rho);
		sum-=temp_vec[i];
	}

	temp_vec2=temp_vec-last_m_vec;

	x= -temp_vec2.transpose() * last_s_Inv.asDiagonal() * (Cov_Mat_Inv ) * last_s_Inv.asDiagonal() *  temp_vec2;
	
	return sum+x;
}

float sl_obj::get_rho_last_prob_new_rel(void)
{
	Eigen::Matrix<float, 150, 1> last_s_Inv, temp_vec, temp_vec2;
	float x;
	float sum=0;
	
	for (int i=0; i<rel_length; i++)
	{
		last_s_Inv[i]=1./last_s_vec[i];
		temp_vec[i]=log(running_A_mean[i].test_mean_rho);
		sum-=temp_vec[i];
	}

	temp_vec2=temp_vec-last_m_vec;

	x= -temp_vec2.transpose() * last_s_Inv.asDiagonal() * (Cov_Mat_Inv ) * last_s_Inv.asDiagonal() *  temp_vec2;
	
	return sum+x;
}

float sl_obj::get_rho_test_prob(void)
{
	Eigen::Matrix<float, 150, 1> last_s_Inv, temp_vec, temp_vec2;
	float x;
	float sum=0;
	
	for (int i=0; i<rel_length; i++)
	{
		last_s_Inv[i]=1./last_s_vec[i];
		temp_vec[i]=log(running_A_mean[i].last_mean_rho);
		sum-=temp_vec[i];
	}

	temp_vec2=temp_vec-test_m_vec;

	x= -temp_vec2.transpose() * last_s_Inv.asDiagonal() * (Cov_Mat_Inv ) * last_s_Inv.asDiagonal() *  temp_vec2;
	
	return sum+x;
}

float sl_obj::get_rho_last_prob_higher(void)
{
	Eigen::Matrix<float, 150, 1> temp_vec, temp_vec2, temp_vec3;
	float x;
	for (int i=0; i<150; i++){temp_vec[i]=0;}
	float sum=0;

	for (int i=0; i<higher_neighbour_slsl.size(); i++){temp_vec+=higher_neighbour_slsl[i]->last_s_inv.asDiagonal()*
								(higher_neighbour_slsl[i]->last_ln_vec-higher_neighbour_slsl[i]->last_m_vec);}

	temp_vec2=last_s_vec.asDiagonal()*(cond_mu_Mat*temp_vec);
	temp_vec3=last_ln_vec-(last_m_vec+temp_vec2);
	for (int i=0; i<rel_length; i++){sum-=last_ln_vec[i];}

	x= -temp_vec3.transpose() * last_s_inv.asDiagonal() * (higher_cond_Mat_Inv ) * last_s_inv.asDiagonal() *  temp_vec3;
	
	return sum+x;
}

float sl_obj::get_rho_test_prob_higher(void)
{
	Eigen::Matrix<float, 150, 1> temp_vec, temp_vec2, temp_vec3;
	float x;
	for (int i=0; i<150; i++){temp_vec[i]=0;}
	float sum=0;

	for (int i=0; i<higher_neighbour_slsl.size(); i++){temp_vec+=higher_neighbour_slsl[i]->last_s_inv.asDiagonal()*
								(higher_neighbour_slsl[i]->last_ln_vec-higher_neighbour_slsl[i]->test_m_vec);}

	temp_vec2=last_s_vec.asDiagonal()*(cond_mu_Mat*temp_vec);
	temp_vec3=last_ln_vec-(test_m_vec+temp_vec2);
	for (int i=0; i<rel_length; i++){sum-=last_ln_vec[i];}

	x= -temp_vec3.transpose() * last_s_inv.asDiagonal() * (higher_cond_Mat_Inv ) * last_s_inv.asDiagonal() *  temp_vec3;
	
	return sum+x;
}


//====================================================================================================



void sl_obj::rho_to_A(void)
{
	running_A_mean[0].rho_to_A(0.);
	for (int it=1; it<running_A_mean.size(); it++){running_A_mean[it].rho_to_A(running_A_mean[it-1].test_mean_A);}
}

void sl_obj::initial_rho_to_A(void)
{
	running_A_mean[0].initial_rho_to_A(0.);
	for (int it=1; it<running_A_mean.size(); it++){running_A_mean[it].initial_rho_to_A(running_A_mean[it-1].last_mean_A);}
}

void sl_obj::make_new_test_m_vec(float s_R, float s_z, float A_0)
{
	vector <float> test_rho;
	test_rho=backup_rho_mean_find(l,b,s_R, s_z, A_0);

	for (int i=0; i<150; i++)
	{
	//	cout << sl_identifier << " " << i << " " << test_m_vec[i] << " " << test_rho[i] << " " << l << " " << b << " " << s_R << " " << s_z << " " << A_0 <<  endl;
		test_m_vec[i]=log(test_rho[i]) - Cov_Mat.coeffRef(i,i)*pow(last_s_vec[i],2)/2. ;
	}
}


void sl_obj::mean_intervals(void)
{
//	#pragma omp parallel for  num_threads(3)
	for (int star_it=0; star_it<star_cat.size(); star_it++){star_cat[star_it].mean_intervals();}

//	#pragma omp parallel for  num_threads(3)
	for (int it=0; it<running_A_mean.size(); it++){running_A_mean[it].mean_intervals();}
}

void sl_obj::neighbour_set(sl_obj * neighbour)
{
	neighbour_slsl.push_back(neighbour);
}

void sl_obj::higher_neighbour_set(sl_obj * neighbour)
{
	higher_neighbour_slsl.push_back(neighbour);
}


void sl_obj::acl_calc(void)
{
//	ofstream acl_out;
//	string dummy_string;
//	dummy_string=rootname+".acl";
//	acl_out.open(dummy_string.c_str(), ios::trunc);
//	acl_out << "# lag acf_s_R acf_s_z acf_A_0" <<endl; 


//	vector <float> acl;
//	//vector <float> new_acl;
//	vector <float> new_acl(int(ceil(0.5*s_R_chain.size())), 0.);
//	vector <float> new_acl2(int(ceil(0.5*s_z_chain.size())), 0.);
//	vector <float> new_acl3(int(ceil(0.5*A_0_chain.size())), 0.);

//	for (int it1=floor(0.5*s_R_chain.size()); it1<s_R_chain.size(); it1++)
//	{
//		for (int lag=0; lag<it1-ceil(0.5*s_R_chain.size()); lag++)
//		{
//			new_acl[lag]+=(s_R_chain[it1]-s_R_mean)*(s_R_chain[it1-lag]-s_R_mean);
//			new_acl2[lag]+=(s_z_chain[it1]-s_z_mean)*(s_z_chain[it1-lag]-s_z_mean);
//			new_acl3[lag]+=(A_0_chain[it1]-A_0_mean)*(A_0_chain[it1-lag]-A_0_mean);
//		}

//	}

//	for (int it1=0; it1<new_acl.size(); it1++)
//	{
//		new_acl[it1]/=(new_acl.size());
//		new_acl2[it1]/=(new_acl2.size());
//		new_acl3[it1]/=(new_acl3.size());
//	}

////	new_acl=star_cat[0].acl_calc();
//	acl=new_acl;

////	for (int star_it=1; star_it<star_cat.size(); star_it++)
////	{
////		new_acl=star_cat[star_it].acl_calc();
////		for (int it=0; it<acl.size(); it++)
////		{
////			acl[it]+=new_acl[it];
////		}
////	}

//	for (int it=0; it<new_acl.size(); it++)
//	{
//		acl_out << it << " " << new_acl[it]<< " " << new_acl2[it]<< " " << new_acl3[it] << endl;
//	}
//	acl_out.close();

//	size_t const half_size = s_R_chain.size() / 2;
//	vector <float> unburnt_s_R(s_R_chain.begin()+half_size, s_R_chain.end());
//	vector <float> unburnt_s_z(s_z_chain.begin()+half_size, s_z_chain.end());
//	vector <float> unburnt_A_0(A_0_chain.begin()+half_size, A_0_chain.end());

//	cout << acl_block(unburnt_s_R) << " " << acl_block(unburnt_s_z) << " " << acl_block(unburnt_A_0) << " " << endl; 
}

float sl_obj::hyperprior_prob_get(vector < vector <float> > internal_rel)
{
	Eigen::Matrix<float, 150, 1> int_vec, temp_vec;
	for (int i=1; i<150; i++){int_vec[i]=internal_rel[i][0];}

	temp_vec=int_vec-last_m_vec;
	return temp_vec.transpose()*temp_vec;
}



void sl_obj::define_cov_mat(void)
{
	typedef Eigen::SparseMatrix<float> spMat;
//	vector <float> mean_rel;

//	vector <float> mean_rho(150,0);

//	mean_rel=backup_A_mean_find(l,b,2500, 125, true);
//	mean_rho[0]=mean_rel[0];

//	for (int i=1; i<150; i++){mean_rho[i]=mean_rel[i]-mean_rel[i-1];}

	typedef Eigen::Triplet<float> T;
	Eigen::SparseMatrix<float> CM(150,150), temp1(150,150);
	std::vector<T> tripletList, tripletList2, tripletList3, trips4;
	tripletList.reserve(150);
	tripletList2.reserve(150);
	tripletList3.reserve(150);

	for (int i=0; i<150; i++){tripletList.push_back(T(i, i, log(1+pow(10.,2*(-1-i/100.))) ) ) ;}
	for (int i=0; i<149; i++)
	{
		tripletList.push_back(T(i, i+1, log(1+pow(10.,(-1-i/100.)+(-1-(i+1)/100.)-2.)) ) ) ;
		tripletList.push_back(T(i+1, i, log(1+pow(10.,(-1-i/100.)+(-1-(i+1)/100.)-2.)) ) ) ;
	}
	CM.setFromTriplets(tripletList.begin(), tripletList.end());

	for (int i=0; i<150; i++)
	{
		last_s_vec[i]=5;
		last_s_inv[i]=1./last_s_vec[i];
//		last_m_vec[i]=log(mean_rho[i]) - CM.coeffRef(i,i)*pow(last_s_vec[i],2)/2. ;
	}

	Cov_Mat=CM;

	Eigen::SimplicialLLT<spMat> chol(CM);
	chol_L=chol.matrixL();

	Eigen::SparseMatrix<float> I(150,150);
	for (int i=0; i<150; i++){tripletList2.push_back(T(i, i, 1. ) ) ;}
	I.setFromTriplets(tripletList2.begin(), tripletList2.end());

	Cov_Mat_Inv = chol.compute(CM).solve(I);

	Eigen::Matrix<float, 150, 150> I_dense;
	I_dense= Eigen::MatrixXf::Identity(150, 150);

	Eigen::SimplicialLLT<spMat> chol2(chol_L);
	chol_L_Inv=chol2.compute(chol_L).solve(I);
	
	//chol_L_Inv=chol_L;
	//chol_L_Inv=chol_L_Inv.triangularView<Eigen::Lower>().solve(I_dense);

	//cout << "stuff: " << chol_L.coeffRef(10,10) << " " << chol_L_Inv.coeffRef(10,10) << endl;
	//cout << chol_L * chol_L_Inv << endl;
	Eigen::SparseMatrix<float> temp_rho(150,150);

	for (int i=0; i<150; i++){tripletList3.push_back(T(i, i, log(1+pow(10.,2*(-1-i/100.)-1.)) ) ) ;}
	temp_rho.setFromTriplets(tripletList3.begin(), tripletList3.end());
	rho_Mat=temp_rho;
//	for (int i=0; i<150; i++){cout << i << " " << rho_Mat.coeffRef(i,i) << " " << Cov_Mat.coeffRef(i,i) << endl;}

	cond_mu_Mat=rho_Mat*Cov_Mat_Inv;
	cond_Mat=(Cov_Mat - neighbour_slsl.size()*rho_Mat*Cov_Mat_Inv*rho_Mat).selfadjointView<Eigen::Lower>();
	higher_cond_Mat=(Cov_Mat - higher_neighbour_slsl.size()*rho_Mat*Cov_Mat_Inv*rho_Mat).selfadjointView<Eigen::Lower>();

	Eigen::SimplicialLLT<spMat> chol_cond(cond_Mat);
	chol_L_cond=chol_cond.matrixL();

	Eigen::SimplicialLLT<spMat> chol3(higher_cond_Mat);
	higher_cond_Mat_Inv=chol3.compute(higher_cond_Mat).solve(I);
	higher_chol_L=chol3.matrixL();

	higher_chol_L_Inv=higher_chol_L.transpose()*higher_cond_Mat_Inv;

	//cout << higher_chol_L*higher_chol_L_Inv << endl;

}


