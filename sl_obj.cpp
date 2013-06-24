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

	neighbour_sl=NULL;
	define_cov_mat();

	vector<bin_obj> running_A_mean(150);

}

sl_obj::sl_obj(string filename, float l_in, float b_in, string datatype)
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

	previous_s_R=2500.;
	previous_s_z=125.;

	//	global_previous_prob=0;
	//	previous_hyperprior_prob=0;
	//	current_hyperprior_prob=0;

	// Read in data	-------------------------------------------------------------------------

	if (datatype=="iphas" || datatype=="IPHAS"){star_cat=iphas_read(filename,r_min,i_min,ha_min,r_max,i_max,ha_max);}
	else if (datatype=="2MASS" || datatype=="2mass"){star_cat=TWOMASS_read(filename,J_min, H_min, K_min,J_max,H_max,K_max);}
	else {cout << "Unrecognised datatype: " << datatype << endl;}

	// Find expected A(d) -------------------------------------------------------------------

	//vector<bin_obj2> backup_A_mean (150);
	backup_A_mean.resize(150);
	backup_A_mean=backup_A_mean_find(l, b, previous_s_R, previous_s_z, true);

	sigma_fac=0.025;
	accepted=0;
	without_change=0;
	thin=200;
	rel_length=150;

	it_num=0.;
	neighbour_sl=NULL;

	define_cov_mat();

	running_A_mean.resize(150);
}



void sl_obj::output_write(void)
{
	ofstream A_out;
	string dummy_string;
	dummy_string=rootname+".td4";
	A_out.open(dummy_string.c_str(), ios::trunc);
	//A_out << "#\tdist\tA\tsigma_A\n";

	for (int x=0; x<running_A_mean.size(); x++)
	{
		A_out 	<< x*100 + 50 << "\t" << running_A_mean[x].final_A << "\t" << running_A_mean[x].final_sd <<"\t"
			<< running_A_mean[x].final_dA<<"\t"<<running_A_mean[x].final_dsd<<endl;
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
	mean_rho=density_find(l, b, s_R_mean, s_z_mean, 100);
	//mean_rho=density_find(l, b, previous_s_R, previous_s_z, 100);
	mean_rel=backup_A_mean_find(l, b, previous_s_R, previous_s_z, false);

	cout << s_R_mean << " " << s_z_mean << endl;

	for (int x=0; x<rho_mean.size(); x++)
	{
		rho_out << x*100 << "\t" << running_A_mean[x].final_rho << "\t" << running_A_mean[x].final_drho <<"\t" << mean_rel[x]-mean_rel[x-1] << "\t" << mean_rel[x-1] << "\n";
	}
	rho_out.close();
	
}


void sl_obj::initial_guess(vector<iso_obj> &isochrones, vector<iso_obj> &guess_set, vector <LF> &LFs)
{
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

// Start from backup_A_mean

	s_R_chain.push_back(previous_s_R);
	s_z_chain.push_back(previous_s_z);


	vector<float> backup_rho_mean;
	float rho_sum=0., Sch_max;
	Sch_max=SFD_read(l, b)*3.1;
	backup_rho_mean=backup_rho_mean_find(l, b, previous_s_R, previous_s_z, 1.);
	for (int it=0; it<backup_rho_mean.size(); it++){rho_sum+=backup_rho_mean[it];}

	for (int it=0; it<running_A_mean.size(); it++){running_A_mean[it].last_mean_rho=backup_rho_mean[it] * Sch_max/rho_sum ;}
	initial_rho_to_A();


	for (int i=0; i<150; i++)
	{

		previous_rel[i][0]=running_A_mean[i].last_mean_A;
		if (i==0) {previous_rel[i][1]=0.4;}//5*previous_rel[i][0];}
		else {previous_rel[i][1]=0.4;}//5*previous_rel[i][0];}//sqrt(pow(previous_rel[i-1][1],2)+pow(previous_rel[i][0]-previous_rel[i-1][0],2));}
		previous_rel[i][3]=sqrt(log(1+pow(previous_rel[i][1]/previous_rel[i][0],2)));
		previous_rel[i][2]=log(previous_rel[i][0])-pow(previous_rel[i][3],2)/2;

	}

	hyperprior_rel=previous_rel;

	previous_internal_rel[0][0]=previous_rel[0][0];//log(previous_rel[0][0]);//
	previous_internal_rel[0][1]=previous_rel[0][1];///previous_internal_rel[0][0];

	running_A_mean[0].test_mean_rho=previous_internal_rel[0][0];	

		
	for (int i=1; i<150; i++)
	{
		float mu_t, sig_t;
		previous_internal_rel[i][0]=previous_rel[i][0]-previous_rel[i-1][0];//log(previous_rel[i][0]-previous_rel[i-1][0]);//
		previous_internal_rel[i][1]=previous_rel[i][1];//sqrt(pow(previous_rel[i][1],2)-pow(previous_rel[i-1][1],2))/previous_internal_rel[i][0];

	}

	


	hyperprior_internal_rel=previous_internal_rel;

	global_A_chain.push_back(previous_rel);

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
		else if (star_cat[it_stars].J-star_cat[it_stars].H < (star_cat[it_stars].H - star_cat[it_stars].K)*1.55 -.05 )
		{
			star_cat.erase(star_cat.begin()+it_stars);
		}
		else if (star_cat[it_stars].J-star_cat[it_stars].H > (star_cat[it_stars].H - star_cat[it_stars].K)*1.55 +1.0 )
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
		star_cat[it_stars].initial_guess(isochrones, guess_set, previous_rel, running_A_mean);
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
		previous_norm_prob+=-LFs[it_LF].LF_prob(previous_rel)*(star_cat.size()+1);
	}



//	if (neighbour_sl)
//	{
//		for (int it=0; it<rel_length; it++)
//		{
//			previous_xsl_prob+=-pow( (log(previous_internal_rel[it][0])-log(1+pow(0.75/previous_internal_rel[it][0],2))/2)/(log(1+pow(0.75/previous_internal_rel[it][0],2))) 
//					- (log(neighbour_sl->previous_internal_rel[it][0])-log(1+pow(0.75/neighbour_sl->previous_internal_rel[it][0],2))/2)/(log(1+pow(0.75/neighbour_sl->previous_internal_rel[it][0],2))) ,2)/(2.*fBm_s);
//		//	previous_xsl_prob+=pow(previous_internal_rel[it][0] - neighbour_sl->previous_internal_rel[it][0],2)/(0.125);

//		}
//	}
}



void sl_obj::dist_redMCMC(vector<iso_obj> &isochrones, vector <LF> &LFs)
{

	while (it_num<150000 )
	{
		update(isochrones, LFs);
		if (it_num/100==floor(it_num/100.)){hyperprior_update();}
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
//		#pragma omp parallel for  num_threads(3) reduction(+:dummy)
		for (int it=0; it<star_cat.size(); it++)
		{
			/*if (gsl_ran_flat(rng_handle, 0, 1)>0.){*/star_cat[it].star_try1(isochrones, l, b, previous_rel, running_A_mean);//};
			dummy+=star_cat[it].last_A_prob;
		}
		global_previous_prob=dummy;

// Now vary hyper-parameters

		move_on=false;	
		threshold=log( gsl_ran_flat(rng_handle, 0, 1) );
		trial_rel=mvn_gen_internal_rel(previous_internal_rel, rel_length);

		float sss=0.01;
		//theta=gsl_ran_flat(rng_handle, 0, 2*PI);
		theta=gsl_ran_gaussian_ziggurat(rng_handle, sss);
		theta_min=theta-2*PI; 
		theta_max=theta;

		while (!move_on)
		{
			for (int i=0; i<rel_length; i++)
			{
				running_A_mean[i].test_mean_rho=exp(log(running_A_mean[i].last_mean_rho)*cos(theta) + log(trial_rel[i][0])*sin(theta));
				internal_rel[i][0] = running_A_mean[i].test_mean_rho;
				internal_rel[i][1] = exp(log(previous_internal_rel[i][1])*cos(theta) + log(trial_rel[i][1])*sin(theta));

//				cout << i << " " << previous_internal_rel[i][0] << " " << running_A_mean[i].last_mean_rho << " " << running_A_mean[i].test_n << endl;
			}
			new_rel=internal_to_external(internal_rel, rel_length);
			rho_to_A();

	// Find probability of this parameter set


			dummy=0;
	//		#pragma omp parallel for  num_threads(3) reduction(+:dummy)
			for (int it=0; it<star_cat.size(); it++)
			{
				proposed_probs[it]=star_cat[it].get_A_prob(star_cat[it].last_iso, star_cat[it].last_A, star_cat[it].last_dist_mod, new_rel);
				dummy+= proposed_probs[it];
			}
			global_current_prob=dummy;

	// Normalisation term

			current_norm_prob=0;
			for (int it_LF=0; it_LF<LFs.size(); it_LF++)
			{
				current_norm_prob+=-LFs[it_LF].LF_prob(new_rel)*(star_cat.size()+1);
			}

	// Neighbour term

	//		if (neighbour_sl)
	//		{
	//			for (int it=0; it<rel_length; it++)
	//			{
	//				current_xsl_prob+=-pow( (log(internal_rel[it][0])-log(1+pow(0.75/internal_rel[it][0],2))/2)/(log(1+pow(0.75/internal_rel[it][0],2))) - (log(neighbour_sl->internal_rel[it][0])-log(1+pow(0.75/neighbour_sl->internal_rel[it][0],2))/2)/(log(1+pow(0.75/neighbour_sl->internal_rel[it][0],2))) ,2)/(2./fBm_s);
	//			//	current_xsl_prob+=pow(internal_rel[it][0] - neighbour_sl->internal_rel[it][0],2)/(0.125);
	//			}
	//		}
	//		else if (!recv_neighbour_rel.empty())
	//		{
	//			for (int it=0; it<rel_length; it++)
	//			{
	//				current_xsl_prob+=-pow( (log(internal_rel[it][0])-log(1+pow(0.75/internal_rel[it][0],2))/2)/(log(1+pow(0.75/internal_rel[it][0],2)))- (log(neighbour_sl->internal_rel[it][0])-log(1+pow(0.75/neighbour_sl->internal_rel[it][0],2))/2)/(log(1+pow(0.75/neighbour_sl->internal_rel[it][0],2))) ,2)/(2./fBm_s);
	//			//	current_xsl_prob+=pow(internal_rel[it][0] -recv_neighbour_rel[it][0],2)/(0.125);
	//			}
	//		}	

	// Metropolis-Hastings algorithm step

			global_transition_prob=0;

			dummy=0;
	//		#pragma omp parallel for  num_threads(3) reduction(+:dummy)
//			for (int it=1; it<150; it++)
//			{
//			// From new to old
//			// mean_A
//				dummy+=log(gsl_ran_lognormal_pdf(previous_internal_rel[it][0], log(internal_rel[it][0])-pow(proposal_sd[it][0],2)/2 ,proposal_sd[it][0]));
//				dummy+=log(gsl_ran_lognormal_pdf(previous_internal_rel[it][1], log(internal_rel[it][1])-pow(proposal_sd[it][1],2)/2 ,proposal_sd[it][1]));
//			// From old to new
//			// mean_A
//				dummy-=log(gsl_ran_lognormal_pdf(internal_rel[it][0], log(previous_internal_rel[it][0])-pow(proposal_sd[it][0],2)/2 ,proposal_sd[it][0]));
//				dummy-=log(gsl_ran_lognormal_pdf(internal_rel[it][1], log(previous_internal_rel[it][1])-pow(proposal_sd[it][1],2)/2 ,proposal_sd[it][1]));
//			}
			global_transition_prob=dummy;
	

// Accept or reject

			accepted++;
			if (global_current_prob+current_hyperprior_prob-global_previous_prob-previous_hyperprior_prob+global_transition_prob+current_xsl_prob-previous_xsl_prob
				 + current_norm_prob-previous_norm_prob>threshold)		// New parameter set better => Accept
			{
				previous_rel=new_rel;
				previous_internal_rel=internal_rel;
				global_previous_prob=global_current_prob;
				previous_hyperprior_prob=current_hyperprior_prob;
				previous_xsl_prob=current_xsl_prob;
				previous_norm_prob=current_norm_prob;

				for (int it=0; it<running_A_mean.size(); it++){running_A_mean[it].accept();}

				for (int stars_it=0; stars_it<star_cat.size(); stars_it++){star_cat[stars_it].last_A_prob=proposed_probs[stars_it];}

				move_on=true;
		
		//		cout << global_A_chain.size() << " " << global_current_prob << " " << internal_rel[50][0] << " " << new_rel[rel_length-1][0] << " " << current_hyperprior_prob << " " << accepted/global_A_chain.size() << " " << theta << endl;
			}

			else 
			{
				//cout << "fail " << global_current_prob << " " << global_previous_prob << " " << global_transition_prob << " " << current_hyperprior_prob << " " << previous_hyperprior_prob << " " << new_rel[10][0] << " " << previous_rel[10][0] << " " << running_A_mean[10].test_mean_A << " " << running_A_mean[10].last_mean_A << endl;//*/
			if (theta>0){theta_max=theta;}
			else {theta_min=theta;}
		//	theta=gsl_ran_flat(rng_handle, theta_min, theta_max);
			sss/=10;	
			theta=gsl_ran_gaussian_ziggurat(rng_handle, sss);
			for (int it=0; it<running_A_mean.size(); it++){running_A_mean[it].reject();}


			}
		}

		if (neighbour_sl){if (it_num/1000.==floor(it_num/1000.)){cout << it_num << " " << global_previous_prob << " " << internal_rel[50][0] << " " << previous_rel[rel_length-1][0] << " " << previous_hyperprior_prob << " " << accepted << " " << accepted/it_num << " " << neighbour_sl->previous_rel[rel_length-1][0] << endl;}}
		else {if (it_num/1000.==floor(it_num/1000.)){cout << it_num << " " << global_previous_prob << " " << internal_rel[50][0] << " " << previous_rel[rel_length-1][0] << " " << previous_hyperprior_prob << " " << previous_norm_prob << " " << previous_s_R << " " << accepted << " " << accepted/it_num << endl;}}

//		if (it_num/10.==floor(it_num/10.))
//		{
//			ofstream trace1;
//			trace1.open("trace1.txt", ios::app);
//			trace1 <<it_num << " " << global_previous_prob << " " << internal_rel[50][0] << " " << previous_rel[80][0] << " " << previous_hyperprior_prob << " " << previous_norm_prob << " " << previous_s_R << " " << accepted << " " << accepted/it_num << " " << global_transition_prob << " " << dummy2 << " " << previous_xsl_prob << " " <<dummy3 << endl;
//			trace1.close();
//		}

		if (floor(it_num/100.)==it_num/100)
		{
			global_A_chain.push_back(previous_rel);
			for (int it=0; it<running_A_mean.size(); it++){running_A_mean[it].chain_push_back();}
		}
		it_num++;
}



vector < vector <float> > sl_obj::gen_internal_rel(vector < vector <float> > old_rel, int rel_length)
{
	vector < vector <float> > new_rel(rel_length,vector <float> (2));

	for (int it=0; it<rel_length; it++)
	{
		proposal_sd[it][0]=sigma_fac;
		proposal_sd[it][1]=sigma_fac/10;
	}

	for (int it=0; it<rel_length; it++)
	{	
		new_rel[it][0]=gsl_ran_lognormal(rng_handle,log(old_rel[it][0])-pow(proposal_sd[it][0],2)/2,proposal_sd[it][0]);
		new_rel[it][1]=gsl_ran_lognormal(rng_handle,log(old_rel[it][1])-pow(proposal_sd[it][1],2)/2,proposal_sd[it][1]);
	}
	return new_rel;
}

vector < vector <float> > sl_obj::mvn_gen_internal_rel(vector < vector <float> > old_rel, int rel_length)
{
	vector < vector <float> > new_rel(rel_length,vector <float> (2));
	Eigen::Matrix<float, 150, 1> int_vec, temp_vec;

	for (int i=0; i<rel_length; i++){temp_vec[i]=gsl_ran_gaussian_ziggurat(rng_handle, 1.);}
	int_vec=temp_vec.transpose()*Cov_Mat + Mean_vec.transpose();

	for (int i=0; i<rel_length; i++)
	{
		new_rel[i][0]=exp(int_vec[i]);
		new_rel[i][1]=0.4;//*new_rel[i][0];
	}
	return new_rel;
}

vector < vector <float> > sl_obj::internal_to_external(vector < vector <float> > int_rel, int rel_length)
{
	vector < vector <float> > ext_rel(rel_length,vector <float> (4));

	ext_rel[0][0]=int_rel[0][0];
	ext_rel[0][1]=0.4;//int_rel[0][1];
	ext_rel[0][3]=sqrt(log(1+pow(ext_rel[0][1]/ext_rel[0][0],2)));
	ext_rel[0][2]=log(ext_rel[0][0])-pow(ext_rel[0][3],2)/2;
	
	for (int it=1; it<rel_length; it++)
	{			

		ext_rel[it][0]=int_rel[it][0]+ext_rel[it-1][0];
		ext_rel[it][1]=int_rel[it][1];
		ext_rel[it][3]=sqrt(log(1+pow(ext_rel[it][1]/ext_rel[it][0],2)));
		ext_rel[it][2]=log(ext_rel[it][0])-pow(ext_rel[it][3],2)/2;
	}
	return ext_rel;
}

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


void sl_obj::hyperprior_update(void)
{
	float current_s_R, current_s_z;
	vector <vector <float> > current_hyperprior_rel (150, vector <float> (4)), current_hyperprior_internal_rel(150, vector <float> (2));
	vector <float> backup_rel;
	current_hyperprior_prob=0;


	current_s_R=previous_s_R+gsl_ran_gaussian_ziggurat(rng_handle,0.0010);
	current_s_z=previous_s_z+gsl_ran_gaussian_ziggurat(rng_handle,0.001);

	backup_rel=backup_A_mean_find(l,b,current_s_R, current_s_z, false);

	for (int it=0; it<150; it++)
	{
		current_hyperprior_rel[it][0]=backup_rel[it];
		current_hyperprior_rel[it][1]=hyperprior_rel[it][1];
		current_hyperprior_rel[it][3]=sqrt(log(1+pow(current_hyperprior_rel[it][1]/current_hyperprior_rel[it][0],2)));
	}

	current_hyperprior_internal_rel[0][0]=current_hyperprior_rel[0][0];
	current_hyperprior_internal_rel[0][1]=current_hyperprior_rel[0][1];
	for (int it=1; it<150; it++)
	{
		current_hyperprior_internal_rel[it][0]=current_hyperprior_rel[it][0]-current_hyperprior_rel[it-1][0];
		current_hyperprior_internal_rel[it][1]=current_hyperprior_rel[it][1];
	}

	if (current_hyperprior_prob>previous_hyperprior_prob)
	{
		hyperprior_rel=current_hyperprior_rel;
		hyperprior_internal_rel=current_hyperprior_internal_rel;
		previous_hyperprior_prob=current_hyperprior_prob;
		previous_s_R=current_s_R;
		previous_s_z=current_s_z;
	}
	else if (exp(current_hyperprior_prob-previous_hyperprior_prob)>gsl_ran_flat(rng_handle, 0, 1))
	{
		hyperprior_rel=current_hyperprior_rel;
		hyperprior_internal_rel=current_hyperprior_internal_rel;
		previous_hyperprior_prob=current_hyperprior_prob;
		previous_s_R=current_s_R;
		previous_s_z=current_s_z;
	}
	s_R_chain.push_back(previous_s_R);
	s_z_chain.push_back(previous_s_z);
}

void sl_obj::mean_intervals(void)
{
//	#pragma omp parallel for  num_threads(3)
	for (int star_it=0; star_it<star_cat.size(); star_it++){star_cat[star_it].mean_intervals();}

//	#pragma omp parallel for  num_threads(3)
	for (int it=0; it<running_A_mean.size(); it++){running_A_mean[it].mean_intervals();}

	float s_R_sum=0., s_z_sum=0;
	for (int it=floor(0.5*s_R_chain.size()); it<s_R_chain.size(); it++)
	{
		s_R_sum+=s_R_chain[it];
		s_z_sum+=s_z_chain[it];
	}
	s_R_mean=s_R_sum/ceil(0.5*s_R_chain.size());
	s_z_mean=s_z_sum/ceil(0.5*s_z_chain.size());
}

void sl_obj::neighbour_set(sl_obj * neighbour)
{
	neighbour_sl=neighbour;
}


void sl_obj::acl_calc(void)
{
	ofstream acl_out;
	string dummy_string;
	dummy_string=rootname+".acl";
	acl_out.open(dummy_string.c_str(), ios::trunc);
	acl_out << "# lag acf" <<endl; 


	vector <float> acl;
	vector <float> new_acl;

	acl=star_cat[0].acl_calc();

	for (int star_it=1; star_it<star_cat.size(); star_it++)
	{
		new_acl=star_cat[star_it].acl_calc();
		for (int it=0; it<acl.size(); it++)
		{
			acl[it]+=new_acl[it];
		}
	}

	for (int it=0; it<acl.size(); it++)
	{
		acl_out << it+1 << " " << acl[it] << endl;
	}
	acl_out.close();
}

float sl_obj::hyperprior_prob_get(vector < vector <float> > internal_rel)
{
	Eigen::Matrix<float, 150, 1> int_vec, temp_vec;
	for (int i=1; i<150; i++){int_vec[i]=internal_rel[i][0];}

	temp_vec=int_vec-Mean_vec;
	return temp_vec.transpose()*temp_vec;

}


void sl_obj::define_cov_mat(void)
{
	typedef Eigen::SparseMatrix<float> spMat;
	vector <float> mean_rel;

	vector <float> mean_rho(150,0);

	mean_rel=backup_A_mean_find(l,b,2500, 125, true);
	mean_rho[0]=mean_rel[0];

	for (int i=1; i<150; i++){mean_rho[i]=mean_rel[i]-mean_rel[i-1];}

	typedef Eigen::Triplet<float> T;
	Eigen::SparseMatrix<float> CM(150,150);
	std::vector<T> tripletList, tripletList2;
	tripletList.reserve(150);

	for (int i=0; i<150; i++){tripletList.push_back(T(i, i, log(1+pow(10.,2*(-1-i/100.))) ) ) ;}
	CM.setFromTriplets(tripletList.begin(), tripletList.end());
	for (int i=0; i<150; i++){Mean_vec[i]=log(mean_rho[i]) - CM.coeffRef(i,i) ;}

	float dummy=0;
	for (int i=0; i<150; i++){dummy+=exp(Mean_vec[i]+CM.coeffRef(i,i)/2);}
					//cout << i << " " << dummy << endl;}

	Cov_Mat=CM;
	//Cov_Mat.coeffRef(0,0)=1.;

	//Eigen::SimplicialLDLT<spMat> chol(CM);

	Eigen::SparseMatrix<float> chol1(150,150);
	for (int i=0; i<150; i++){tripletList2.push_back(T(i, i, sqrt(Cov_Mat.coeffRef(i,i)) ) ) ; }
	chol1.setFromTriplets(tripletList2.begin(), tripletList2.end());
	chol=chol1;

	Eigen::SimplicialLDLT<spMat> chol2(CM);

	//cout << Cov_Mat.nonZeros() << " CM " << Cov_Mat.coeffRef(0,0) << endl;
}


