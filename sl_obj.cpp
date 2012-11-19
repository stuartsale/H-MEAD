#include "sl_obj.h"

sl_obj::sl_obj(void)
{
	// set default MIN vals
	r_min=13.5; 
	i_min=12.0;
	ha_min=12.0;				
	// set default MAX vals
	r_max=21.;
	i_max=20.;
	ha_max=20.;

	sigma_fac=0.05;
	accepted=0;
	without_change=0;
	thin=200;
	rel_length=150;

	it_num=0.;

	neighbour_sl=NULL;

}

sl_obj::sl_obj(string filename, float l_in, float b_in, string datatype)
{
	// Set up variables ---------------------------------------------------------------------
	l=l_in;
	b=b_in;
	rootname=filename;
	rootname.erase(rootname.size()-4);

	// set default MIN vals
	r_min=13.5; 
	i_min=12.0;
	ha_min=12.0;				
	// set default MAX vals
	r_max=21.;
	i_max=20.;
	ha_max=20.;

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

	A_mean.resize(150);

	sigma_fac=0.025;
	accepted=0;
	without_change=0;
	thin=200;
	rel_length=150;

	it_num=0.;
	neighbour_sl=NULL;


}



void sl_obj::output_write(void)
{
	ofstream A_out;
	string dummy_string;
	dummy_string=rootname+".td4";
	A_out.open(dummy_string.c_str(), ios::trunc);
	//A_out << "#\tdist\tA\tsigma_A\n";

	for (int x=0; x<A_mean.size(); x++)
	{
		A_out << x*100 + 50 << "\t" << A_mean[x].mean_A << "\t" << A_mean[x].sigma <<"\t"<<A_mean[x].d_mean<<"\t"<<A_mean[x].d_sigma<<"\t"
			<<A_mean[x].size<<"\t"<<A_mean[x].error_measure<<"\t"<<A_mean[x].sum<<"\t"<<A_mean[x].diff<<"\t"<<A_mean[x].d_diff<<"\n";
	}
	A_out.close();


	ofstream output;
	dummy_string=rootname+"-090.dat";
	output.open(dummy_string.c_str(), ios::trunc);
	output << "#\tr\ti\tha\tr_i0\tdist\tA\tdistbin\td_A\td_r_i0\td_dist\td_r\td_i\td_ha\tmag_weight\tprob\tA_prob\tMi\tlogAge\tfeh\td_Mi\td_lagAge\td_feh\tlogT\tlogg\trx\tix\thax\n" ;
	for (int y=0; y<star_cat.size(); y++)
	{
		
		output << star_cat[y].r << "\t" << star_cat[y].i << "\t" << star_cat[y].ha << "\t"<< star_cat[y].r_i0 << "\t" << star_cat[y].dist << "\t" << star_cat[y].A << "\t" 
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
	//mean_rho=density_find(l, b, s_R_mean, s_z_mean, 100);
	mean_rho=density_find(l, b, previous_s_R, previous_s_z, 100);
	mean_rel=backup_A_mean_find(l, b, previous_s_R, previous_s_z, false);

	cout << s_R_mean << " " << s_z_mean << endl;

	//rho_out << 25 << "\t" << (A_mean[0].mean_A)/0.05 << "\t" << mean_rho[0] << "\n";
	rho_out << 25 << "\t" << (previous_rel[0][0])/0.05 << "\t0\t" << mean_rho[0] << "\t0\n";
	for (int x=1; x<A_mean.size(); x++)
	{
	//	rho_out << x*100 << "\t" << (A_mean[x].mean_A - A_mean[x-1].mean_A)/0.1 << "\t" << mean_rho[x] << "\n";
		rho_out << x*100 << "\t" << (previous_internal_rel[x][0]) << "\t" << (previous_rel[x-1][0]) <<"\t" << mean_rel[x]-mean_rel[x-1] << "\t" << mean_rel[x-1] << "\n";
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

	for (int i=0; i<150; i++)
	{

		previous_rel[i][0]=backup_A_mean[i];
		if (i==0) {previous_rel[i][1]=0.4;}//5*previous_rel[i][0];}
		else {previous_rel[i][1]=0.4;}//5*previous_rel[i][0];}//sqrt(pow(previous_rel[i-1][1],2)+pow(previous_rel[i][0]-previous_rel[i-1][0],2));}
		previous_rel[i][3]=sqrt(log(1+pow(previous_rel[i][1]/previous_rel[i][0],2)));
		previous_rel[i][2]=log(previous_rel[i][0])-pow(previous_rel[i][3],2)/2;

	//	previous_hyperprior_prob+=log(previous_rel[i][1]/(exp(pow(previous_rel[i][3],2))*pow(previous_rel[i][0],3)*previous_rel[i][3])) - 2*log(previous_rel[i][3]);
	}

	hyperprior_rel=previous_rel;

	previous_internal_rel[0][0]=previous_rel[0][0];//log(previous_rel[0][0]);//
	previous_internal_rel[0][1]=previous_rel[0][1];///previous_internal_rel[0][0];

	previous_hyperprior_prob+=-log(previous_internal_rel[0][0])- pow(log(previous_internal_rel[0][0]/previous_internal_rel[0][0]),2)/(2*pow(fBm_s,2))
		
	for (int i=1; i<150; i++)
	{
		float mu_t, sig_t;
		previous_internal_rel[i][0]=previous_rel[i][0]-previous_rel[i-1][0];//log(previous_rel[i][0]-previous_rel[i-1][0]);//
		previous_internal_rel[i][1]=previous_rel[i][1];//sqrt(pow(previous_rel[i][1],2)-pow(previous_rel[i-1][1],2))/previous_internal_rel[i][0];
		

	//	previous_hyperprior_prob+=-pow( ((log(previous_internal_rel[i][0])-(log(1+pow(1/*previous_internal_rel[i][1]/previous_internal_rel[i][0]*/,2))/2)) - (log(previous_internal_rel[i][0])-(log(1+pow(1/*previous_internal_rel[i][1]/previous_internal_rel[i][0]*/,2))/2))  )/(1) - (previous_rel[i-1][2]-hyperprior_rel[i-1][2])/(1/*previous_rel[i-1][3]*/),2)/(2.*fBm_s);
		mu_t=
		sig_t=

		previous_hyperprior_prob+=-log(previous_internal_rel[0][0]/previous_internal_rel[0][0])- pow(log(previous_internal_rel[0][0]/previous_internal_rel[0][0]),2)/(2*pow(fBm_s,2))
	}

	hyperprior_internal_rel=previous_internal_rel;

	global_A_chain.push_back(previous_rel);

// Dump sources not in required region of c-c

	int it_stars=0;
	while (it_stars<star_cat.size())
	{
		if (star_cat[it_stars].r-star_cat[it_stars].ha>guess_set[0].redline(star_cat[it_stars].r-star_cat[it_stars].i) || star_cat[it_stars].r-star_cat[it_stars].ha<guess_set[guess_set.size()-1].redline(star_cat[it_stars].r-star_cat[it_stars].i))
		{
			star_cat.erase(star_cat.begin()+it_stars);
		}
		else {it_stars++;}
	}

// Make initial guess

	it_stars=0;
	while (it_stars<star_cat.size())
	{
		star_cat[it_stars].initial_guess(isochrones, guess_set, previous_rel);
		if (star_cat[it_stars].last_A<0){star_cat[it_stars].last_A = 0.02;} 
		it_stars++;
	}


	for (int star_it=0; star_it<star_cat.size(); star_it++){global_previous_prob+=star_cat[star_it].last_A_prob;}

	proposed_probs.resize(star_cat.size());

	previous_norm_prob=0;
	for (int it_LF=0; it_LF<LFs.size(); it_LF++)
	{
		previous_norm_prob+=-LFs[it_LF].LF_prob(previous_rel)*(star_cat.size()+1);
	}

	if (neighbour_sl)
	{
		for (int it=0; it<rel_length; it++)
		{
			previous_xsl_prob+=-pow( (log(previous_internal_rel[it][0])-log(1+pow(0.75/previous_internal_rel[it][0],2))/2)/(log(1+pow(0.75/previous_internal_rel[it][0],2))) 
					- (log(neighbour_sl->previous_internal_rel[it][0])-log(1+pow(0.75/neighbour_sl->previous_internal_rel[it][0],2))/2)/(log(1+pow(0.75/neighbour_sl->previous_internal_rel[it][0],2))) ,2)/(2.*fBm_s);
		//	previous_xsl_prob+=pow(previous_internal_rel[it][0] - neighbour_sl->previous_internal_rel[it][0],2)/(0.125);

		}
	}
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
			/*if (gsl_ran_flat(rng_handle, 0, 1)>0.){*/star_cat[it].star_try1(isochrones, l, b, previous_rel);//};
			dummy+=star_cat[it].last_A_prob;
		}
		global_previous_prob=dummy;

	//	cout << star_cat[161].last_A << " " << star_cat[161].last_dist_mod << " " << star_cat[161].last_prob << " " << star_cat[161].last_iso.logT << " " << star_cat[161].last_iso.logg << " " << star_cat[161].last_iso.r0-star_cat[161].last_iso.i0 << " " << star_cat[161].last_iso.r0-star_cat[161].last_iso.ha0 << " " << it_num  << " " << star_cat[161].last_iso.Mi << " " << log(star_cat[161].last_iso.Jac) << " " << log(star_cat[161].last_iso.IMF())  << " " 
//<< star_cat[161].last_iso.r0+star_cat[161].last_iso.u*pow(star_cat[161].last_A,2)+star_cat[161].last_iso.v*star_cat[161].last_A+star_cat[161].last_iso.w+star_cat[161].last_dist_mod  << " " 
//<< star_cat[161].last_iso.i0+star_cat[161].last_iso.u_i*pow(star_cat[161].last_A,2)+star_cat[161].last_iso.v_i*star_cat[161].last_A+star_cat[161].last_iso.w_i+star_cat[161].last_dist_mod  << " " << star_cat[161].last_iso.ha0+star_cat[161].last_iso.u_ha*pow(star_cat[161].last_A,2)+star_cat[161].last_iso.v_ha*star_cat[161].last_A+star_cat[161].last_iso.w_ha+star_cat[161].last_dist_mod << endl;

// Now vary hyper-parameters

		vector < vector <float> > new_rel(rel_length,vector <float> (4));

		for (int it=0; it<rel_length; it++)
		{
			proposal_sd[it][0]=sigma_fac;
			proposal_sd[it][1]=sigma_fac/10;
		}
		

		for (int it=0; it<rel_length; it++)
		{	
			internal_rel[it][0]=gsl_ran_lognormal(rng_handle,log(previous_internal_rel[it][0])-pow(proposal_sd[it][0],2)/2,proposal_sd[it][0]);
			//while (internal_rel[it][0]>0.5){internal_rel[it][0]=gsl_ran_lognormal(rng_handle,log(previous_internal_rel[it][0])-pow(proposal_sd[it][0],2)/2,proposal_sd[it][0]);}
			internal_rel[it][1]=gsl_ran_lognormal(rng_handle,log(previous_internal_rel[it][1])-pow(proposal_sd[it][1],2)/2,proposal_sd[it][1]);
			
                        //current_hyperprior_prob+=log(gsl_ran_lognormal_pdf(internal_rel[it][0],log(hyperprior_internal_rel[it][0]),1.));
		}

		new_rel[0][0]=internal_rel[0][0];
		new_rel[0][1]=internal_rel[0][1];//*internal_rel[0][0];
		new_rel[0][3]=sqrt(log(1+pow(new_rel[0][1]/new_rel[0][0],2)));
		new_rel[0][2]=log(new_rel[0][0])-pow(new_rel[0][3],2)/2;
		
		for (int it=1; it<rel_length; it++)
		{			

			new_rel[it][0]=internal_rel[it][0]+new_rel[it-1][0];//exp(internal_rel[it][0])+new_rel[it-1][0];//
			new_rel[it][1]=internal_rel[it][1];//sqrt(pow(internal_rel[it][1]*internal_rel[it][0],2)+pow(new_rel[it-1][1],2));
			new_rel[it][3]=sqrt(log(1+pow(new_rel[it][1]/new_rel[it][0],2)));
			new_rel[it][2]=log(new_rel[it][0])-pow(new_rel[it][3],2)/2;

			current_hyperprior_prob+=-pow( ((log(internal_rel[it][0])-(log(1+pow(1/*internal_rel[it][1]/internal_rel[it][0]*/,2))/2)) - (log(hyperprior_internal_rel[it][0])-(log(1+pow(1/*hyperprior_internal_rel[it][1]/hyperprior_internal_rel[it][0]*/,2))/2))  )/(1) - (new_rel[it-1][2]-hyperprior_rel[it-1][2])/(/*hyperprior_rel[it-1][3]*/1),2)/(2.*fBm_s);
			//current_hyperprior_prob+=-pow( ((log(internal_rel[it][0])-(log(1+pow(1/*internal_rel[it][1]/internal_rel[it][0]*/,2))/2)) - (log(hyperprior_internal_rel[it][0])-(log(1+pow(1/*hyperprior_internal_rel[it][1]/hyperprior_internal_rel[it][0]*/,2))/2))  )/(1) - ((log(internal_rel[it-1][0])-(log(1+pow(1/*internal_rel[it-1][1]/internal_rel[it-1][0]*/,2))/2)) - (log(hyperprior_internal_rel[it-1][0])-(log(1+pow(1/*hyperprior_internal_rel[it-1][1]/hyperprior_internal_rel[it-1][0]*/,2))/2))  )/(1) ,2)/(2.*fBm_s);
		}

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

		if (neighbour_sl)
		{
			for (int it=0; it<rel_length; it++)
			{
				current_xsl_prob+=-pow( (log(internal_rel[it][0])-log(1+pow(0.75/internal_rel[it][0],2))/2)/(log(1+pow(0.75/internal_rel[it][0],2))) - (log(neighbour_sl->internal_rel[it][0])-log(1+pow(0.75/neighbour_sl->internal_rel[it][0],2))/2)/(log(1+pow(0.75/neighbour_sl->internal_rel[it][0],2))) ,2)/(2./fBm_s);
			//	current_xsl_prob+=pow(internal_rel[it][0] - neighbour_sl->internal_rel[it][0],2)/(0.125);
			}
		}
		else if (!recv_neighbour_rel.empty())
		{
			for (int it=0; it<rel_length; it++)
			{
				current_xsl_prob+=-pow( (log(internal_rel[it][0])-log(1+pow(0.75/internal_rel[it][0],2))/2)/(log(1+pow(0.75/internal_rel[it][0],2)))- (log(neighbour_sl->internal_rel[it][0])-log(1+pow(0.75/neighbour_sl->internal_rel[it][0],2))/2)/(log(1+pow(0.75/neighbour_sl->internal_rel[it][0],2))) ,2)/(2./fBm_s);
			//	current_xsl_prob+=pow(internal_rel[it][0] -recv_neighbour_rel[it][0],2)/(0.125);
			}
		}	

// Metropolis-Hastings algorithm step

		global_transition_prob=0;

		dummy=0;
//		#pragma omp parallel for  num_threads(3) reduction(+:dummy)
		for (int it=1; it<150; it++)
		{
		// From new to old
		// mean_A
			dummy+=log(gsl_ran_lognormal_pdf(previous_internal_rel[it][0], log(internal_rel[it][0])-pow(proposal_sd[it][0],2)/2 ,proposal_sd[it][0]));
			dummy+=log(gsl_ran_lognormal_pdf(previous_internal_rel[it][1], log(internal_rel[it][1])-pow(proposal_sd[it][1],2)/2 ,proposal_sd[it][1]));
		// From old to new
		// mean_A
			dummy-=log(gsl_ran_lognormal_pdf(internal_rel[it][0], log(previous_internal_rel[it][0])-pow(proposal_sd[it][0],2)/2 ,proposal_sd[it][0]));
			dummy-=log(gsl_ran_lognormal_pdf(internal_rel[it][1], log(previous_internal_rel[it][1])-pow(proposal_sd[it][1],2)/2 ,proposal_sd[it][1]));
		}
		global_transition_prob=dummy;

// Accept or reject

		float dummy2;
		dummy2=new_rel[149][0]-previous_rel[149][0];	
		int dummy3;
	
		if (global_current_prob+current_hyperprior_prob-global_previous_prob-previous_hyperprior_prob+global_transition_prob+current_xsl_prob-previous_xsl_prob
			 + current_norm_prob-previous_norm_prob>0)		// New parameter set better => Accept
		{
			previous_rel=new_rel;
			previous_internal_rel=internal_rel;
			global_previous_prob=global_current_prob;
			previous_hyperprior_prob=current_hyperprior_prob;
			previous_xsl_prob=current_xsl_prob;
			previous_norm_prob=current_norm_prob;
			without_change=0;
			accepted++;
			dummy3=1;

			for (int stars_it=0; stars_it<star_cat.size(); stars_it++){star_cat[stars_it].last_A_prob=proposed_probs[stars_it];}
		
		//	cout << global_A_chain.size() << " " << global_current_prob << " " << internal_rel[50][0] << " " << new_rel[rel_length-1][0] << " " << current_hyperprior_prob << " " << accepted/global_A_chain.size() << endl;
		}

		else if (exp(global_current_prob+current_hyperprior_prob-global_previous_prob-previous_hyperprior_prob+global_transition_prob+current_xsl_prob-previous_xsl_prob
			+ current_norm_prob-previous_norm_prob)>gsl_ran_flat(rng_handle, 0, 1))	// New set worse => accept with P=P(new)/P(old)
		{
			previous_rel=new_rel;
			previous_internal_rel=internal_rel;
			global_previous_prob=global_current_prob;
			previous_hyperprior_prob=current_hyperprior_prob;
			previous_xsl_prob=current_xsl_prob;
			previous_norm_prob=current_norm_prob;
			without_change=0;
			accepted++;
			dummy3=1;

			for (int stars_it=0; stars_it<star_cat.size(); stars_it++){star_cat[stars_it].last_A_prob=proposed_probs[stars_it];}

		//	cout << global_A_chain.size() << " " << global_current_prob << " " << internal_rel[50][0] << " " << new_rel[rel_length-1][0] << " " << current_hyperprior_prob << " " << accepted/global_A_chain.size() <<  endl;
		}
		else 
		{
			without_change++;
			dummy3=0;
		//	cout << "fail " << global_current_prob << " " << global_previous_prob << " " << global_transition_prob << " " << current_hyperprior_prob << " " << previous_hyperprior_prob << " " << star_cat.size() << endl;//*/

		}
		if (neighbour_sl){if (it_num/1000.==floor(it_num/1000.)){cout << it_num << " " << global_previous_prob << " " << internal_rel[50][0] << " " << previous_rel[rel_length-1][0] << " " << previous_hyperprior_prob << " " << accepted << " " << accepted/it_num << " " << neighbour_sl->previous_rel[rel_length-1][0] << endl;}}
		else {if (it_num/1000.==floor(it_num/1000.)){cout << it_num << " " << global_previous_prob << " " << internal_rel[50][0] << " " << previous_rel[rel_length-1][0] << " " << previous_hyperprior_prob << " " << previous_norm_prob << " " << previous_s_R << " " << accepted << " " << accepted/it_num << endl;}}

		if (it_num/10.==floor(it_num/10.))
		{
			ofstream trace1;
			trace1.open("trace1.txt", ios::app);
			trace1 <<it_num << " " << global_previous_prob << " " << internal_rel[50][0] << " " << previous_rel[80][0] << " " << previous_hyperprior_prob << " " << previous_norm_prob << " " << previous_s_R << " " << accepted << " " << accepted/it_num << " " << global_transition_prob << " " << dummy2 << " " << previous_xsl_prob << " " <<dummy3 << endl;
			trace1.close();
		}

		if (floor(it_num/100.)==it_num/100){global_A_chain.push_back(previous_rel);}
		it_num++;
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
		current_hyperprior_rel[it][2]=log(current_hyperprior_rel[it][0])-pow(current_hyperprior_rel[it][3],2)/2;
	}

	current_hyperprior_internal_rel[0][0]=current_hyperprior_rel[0][0];
	current_hyperprior_internal_rel[0][1]=current_hyperprior_rel[0][1];
	for (int it=1; it<150; it++)
	{
		current_hyperprior_internal_rel[it][0]=current_hyperprior_rel[it][0]-current_hyperprior_rel[it-1][0];
		current_hyperprior_internal_rel[it][1]=current_hyperprior_rel[it][1];
		current_hyperprior_prob+=-pow( ((log(previous_internal_rel[it][0])-(log(1+pow(1/*previous_internal_rel[it][1]/previous_internal_rel[it][0]*/,2))/2)) - (log(current_hyperprior_internal_rel[it][0])-(log(1+pow(1/*current_hyperprior_internal_rel[it][1]/current_hyperprior_internal_rel[it][0]*/,2))/2))  )/(1) - (previous_rel[it-1][2]-current_hyperprior_rel[it-1][2])/(/*current_hyperprior_rel[it-1][3]*/1),2)/(2.*fBm_s);
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

	//#pragma omp parallel for  num_threads(3)
	for (int it=0; it<150; it++)
	{
		float A_sum=0., sigma_sum=0.;
		for (int m=floor(0.50*global_A_chain.size()); m<global_A_chain.size(); m++)
		{
			A_sum+=global_A_chain[m][it][0];
			sigma_sum+=log(global_A_chain[m][it][1]);
		}
		A_mean[it].mean_A=A_sum/ceil(0.50*global_A_chain.size());
		A_mean[it].sigma=exp(sigma_sum/ceil(0.50*global_A_chain.size()));

		vector <float> A_diffs, sigma_diffs;
		for (int m=floor(0.50*global_A_chain.size()); m<global_A_chain.size(); m++)
		{
			A_diffs.push_back(abs(global_A_chain[m][it][0]-A_mean[it].mean_A));
			sigma_diffs.push_back(abs(global_A_chain[m][it][1]-A_mean[it].sigma));
		}
		sort(A_diffs.begin(),A_diffs.end());
		sort(sigma_diffs.begin(), sigma_diffs.end());

		A_mean[it].d_mean=A_diffs[int(0.682*A_diffs.size())];
		A_mean[it].d_sigma=sigma_diffs[int(0.682*A_diffs.size())];

		A_mean[it].error_measure=sqrt(pow(A_mean[it].d_mean/A_mean[it].mean_A,2)+pow(A_mean[it].d_sigma/A_mean[it].sigma,2));
	}

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


