#include "sl_obj.h"

sl_obj::sl_obj(string filename, double l_in, double b_in)
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
	r_max=0.;
	i_max=0.;
	ha_max=0.;

	//	global_previous_prob=0;
	//	previous_hyperprior_prob=0;
	//	current_hyperprior_prob=0;

	// Read in data	-------------------------------------------------------------------------

	star_cat=iphas_read(filename,r_min,i_min,ha_min,r_max,i_max,ha_max);

	// Find expected A(d) -------------------------------------------------------------------

	//vector<bin_obj2> backup_A_mean (150);
	backup_A_mean.resize(150);
	backup_A_mean=backup_A_mean_find(l, b);
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
}



void sl_obj::dist_redMCMC(vector<iso_obj> &isochrones, vector<iso_obj> &guess_set, double ri_min, double ri_max)
{

	double sigma_fac=0.05, accepted=0;
// Set up

	int without_change=0;
	int thin=200;

	vector <vector <vector <double> > > global_A_chain;

	double global_previous_prob=0;
	double previous_hyperprior_prob=0, current_hyperprior_prob=0;
	double global_current_prob, global_transition_prob;
	double sigma2_LN, mu_LN;

	vector < vector<double> > proposal_sd (150, vector <double> (2));
	vector < vector <double> > previous_rel (150, vector <double> (4));
	vector < vector <double> > internal_rel (150, vector <double> (2));
	vector < vector <double> > previous_internal_rel (150, vector <double> (2));
	vector < vector <double> > first_internal_rel (150, vector <double> (2));

	vector <double> initial_dists;
	int rel_length=150;

// Start from backup_A_mean


	for (int i=0; i<150; i++)
	{

		previous_rel[i][0]=backup_A_mean[i].mean_A;
		if (i==0) {previous_rel[i][1]=0.4;}//5*previous_rel[i][0];}
		else {previous_rel[i][1]=0.4;}//5*previous_rel[i][0];}//sqrt(pow(previous_rel[i-1][1],2)+pow(previous_rel[i][0]-previous_rel[i-1][0],2));}
		previous_rel[i][3]=sqrt(log(1+pow(previous_rel[i][1]/previous_rel[i][0],2)));
		previous_rel[i][2]=log(previous_rel[i][0])-pow(previous_rel[i][3],2)/2;

		previous_hyperprior_prob+=log(previous_rel[i][1]/(exp(pow(previous_rel[i][3],2))*pow(previous_rel[i][0],3)*previous_rel[i][3])) - 2*log(previous_rel[i][3]);
	}

	previous_internal_rel[0][0]=previous_rel[0][0];//log(previous_rel[0][0]);//
	previous_internal_rel[0][1]=previous_rel[0][1];///previous_internal_rel[0][0];

	previous_hyperprior_prob+=log(gsl_ran_lognormal_pdf(previous_internal_rel[0][0],log(previous_internal_rel[0][0])-0.75,1.5));
		
//	previous_hyperprior_prob+=-1*log(previous_internal_rel[0][0]);//-previous_internal_rel[0][0];//
	for (int i=1; i<150; i++)
	{
		previous_internal_rel[i][0]=previous_rel[i][0]-previous_rel[i-1][0];//log(previous_rel[i][0]-previous_rel[i-1][0]);//
		previous_internal_rel[i][1]=previous_rel[i][1];//sqrt(pow(previous_rel[i][1],2)-pow(previous_rel[i-1][1],2))/previous_internal_rel[i][0];
		

		previous_hyperprior_prob+=log(gsl_ran_lognormal_pdf(previous_internal_rel[i][0],log(previous_internal_rel[i][0])-0.75,1.5));
	//	previous_hyperprior_prob+=-1*log(previous_internal_rel[i][0]);//-log(previous_internal_rel[i][0]);//
	//	previous_hyperprior_prob+=log(gsl_ran_lognormal_pdf(sqrt(pow(previous_internal_rel[i][1]/previous_internal_rel[i-1][1],2)*(i+1)/i-1.)*previous_rel[i-1][0]/previous_internal_rel[i][0],1.38629, 1.6651));
	//	previous_hyperprior_prob+=log(gsl_ran_lognormal_pdf(previous_internal_rel[i][1],0.34657359,  0.832554611)); 
	//	previous_hyperprior_prob+=log(gsl_ran_lognormal_pdf(previous_internal_rel[i][0]/previous_internal_rel[i-1][0], log(previous_internal_rel[i][0]/previous_internal_rel[i-1][0])-0.02, 0.02));
	}

	first_internal_rel=previous_internal_rel;

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
		if (star_cat[it_stars].last_A>0 && star_cat[it_stars].last_iso.r0-star_cat[it_stars].last_iso.i0> ri_min  && star_cat[it_stars].last_iso.r0-star_cat[it_stars].last_iso.i0<ri_max){initial_dists.push_back(star_cat[it_stars].last_dist_mod);}
		else if (star_cat[it_stars].last_iso.r0-star_cat[it_stars].last_iso.i0> ri_min  && star_cat[it_stars].last_iso.r0-star_cat[it_stars].last_iso.i0<ri_max){initial_dists.push_back(star_cat[it_stars].last_dist_mod);star_cat[it_stars].last_A = 0.02;} 
		it_stars++;
	}

	sort(initial_dists.begin(), initial_dists.end());

//	cout << "95th percentile on dist is: " << pow(10,initial_dists[int(initial_dists.size()*0.95)]/5+1) << " kpc away" << endl;
//	cout << "5th percentile on dist is: " << pow(10,initial_dists[int(initial_dists.size()*0.05)]/5+1) << " kpc away" << endl;

	for (int star_it=0; star_it<star_cat.size(); star_it++){global_previous_prob+=star_cat[star_it].last_A_prob;}
	
// Loop through this section

	vector <double> proposed_probs(star_cat.size());

	float it_num=0.;

	while (it_num<150000 )
	{
		global_current_prob=0;
		global_transition_prob=0;
		global_previous_prob=0;
		current_hyperprior_prob=0;

// First vary parameters for each star

		#pragma omp parallel for  num_threads(3) reduction(+:global_previous_prob)
		for (int it=0; it<star_cat.size(); it++)
		{
			/*if (gsl_ran_flat(rng_handle, 0, 1)>0.){*/star_cat[it].star_try1(isochrones, l, b, previous_rel);//};
			global_previous_prob+=star_cat[it].last_A_prob;
		}

	//	cout << star_cat[161].last_A << " " << star_cat[161].last_dist_mod << " " << star_cat[161].last_prob << " " << star_cat[161].last_iso.logT << " " << star_cat[161].last_iso.logg << " " << star_cat[161].last_iso.r0-star_cat[161].last_iso.i0 << " " << star_cat[161].last_iso.r0-star_cat[161].last_iso.ha0 << " " << it_num  << " " << star_cat[161].last_iso.Mi << " " << log(star_cat[161].last_iso.Jac) << " " << log(star_cat[161].last_iso.IMF())  << " " 
//<< star_cat[161].last_iso.r0+star_cat[161].last_iso.u*pow(star_cat[161].last_A,2)+star_cat[161].last_iso.v*star_cat[161].last_A+star_cat[161].last_iso.w+star_cat[161].last_dist_mod  << " " 
//<< star_cat[161].last_iso.i0+star_cat[161].last_iso.u_i*pow(star_cat[161].last_A,2)+star_cat[161].last_iso.v_i*star_cat[161].last_A+star_cat[161].last_iso.w_i+star_cat[161].last_dist_mod  << " " << star_cat[161].last_iso.ha0+star_cat[161].last_iso.u_ha*pow(star_cat[161].last_A,2)+star_cat[161].last_iso.v_ha*star_cat[161].last_A+star_cat[161].last_iso.w_ha+star_cat[161].last_dist_mod << endl;

// Now vary hyper-parameters

		vector < vector <double> > new_rel(rel_length,vector <double> (4));

		for (int it=0; it<rel_length; it++)
		{
			proposal_sd[it][0]=sigma_fac;
			proposal_sd[it][1]=sigma_fac/10;
		}
		

		for (int it=0; it<rel_length; it++)
		{	
			internal_rel[it][0]=gsl_ran_lognormal(rng_handle,log(previous_internal_rel[it][0])-pow(proposal_sd[it][0],2)/2,proposal_sd[it][0]);
			while (internal_rel[it][0]>0.5){internal_rel[it][0]=gsl_ran_lognormal(rng_handle,log(previous_internal_rel[it][0])-pow(proposal_sd[it][0],2)/2,proposal_sd[it][0]);}
			internal_rel[it][1]=gsl_ran_lognormal(rng_handle,log(previous_internal_rel[it][1])-pow(proposal_sd[it][1],2)/2,proposal_sd[it][1]);
			
			current_hyperprior_prob+=log(gsl_ran_lognormal_pdf(internal_rel[it][0],log(first_internal_rel[it][0])-0.75,1.5));

		}

		new_rel[0][0]=internal_rel[0][0];
		new_rel[0][1]=0.4;//pow(new_rel[0][0], 0.4)/3;//internal_rel[0][1];//*internal_rel[0][0];
		new_rel[0][3]=sqrt(log(1+pow(new_rel[0][1]/new_rel[0][0],2)));
		new_rel[0][2]=log(new_rel[0][0])-pow(new_rel[0][3],2)/2;
		current_hyperprior_prob+=+log(new_rel[0][1]/(exp(pow(new_rel[0][3],2))*pow(new_rel[0][0],3)*new_rel[0][3]))- 2*log(new_rel[0][3]);
		
		for (int it=1; it<rel_length; it++)
		{			

			new_rel[it][0]=internal_rel[it][0]+new_rel[it-1][0];//exp(internal_rel[it][0])+new_rel[it-1][0];//
			new_rel[it][1]=0.4;//pow(new_rel[it][0], 0.4)/3;//internal_rel[it][1];//sqrt(pow(internal_rel[it][1]*internal_rel[it][0],2)+pow(new_rel[it-1][1],2));
			new_rel[it][3]=sqrt(log(1+pow(new_rel[it][1]/new_rel[it][0],2)));
			new_rel[it][2]=log(new_rel[it][0])-pow(new_rel[it][3],2)/2;

			//if (new_rel[it][1]>2*new_rel[it][0]){current_hyperprior_prob-=1E6;}

			//current_hyperprior_prob+=log(gsl_ran_lognormal_pdf(internal_rel[it][1],0.34657359,  0.832554611));
	//		current_hyperprior_prob+=log(gsl_ran_lognormal_pdf(internal_rel[it][0]/internal_rel[it-1][0], log(first_internal_rel[it][0]/first_internal_rel[it-1][0])-0.0002, 0.02));
			current_hyperprior_prob+=+log(new_rel[it][1]/(exp(pow(new_rel[it][3],2))*pow(new_rel[it][0],3)*new_rel[it][3]))- 2*log(new_rel[it][3]);

		}

// Find probability of this parameter set


		#pragma omp parallel for  num_threads(3) reduction(+:global_current_prob)
		for (int it=0; it<star_cat.size(); it++)
		{
			proposed_probs[it]=star_cat[it].get_A_prob(star_cat[it].last_iso, star_cat[it].last_A, star_cat[it].last_dist_mod, new_rel);
			global_current_prob+= proposed_probs[it];
		}

// Metropolis-Hastings algorithm step

		global_transition_prob=0;

		#pragma omp parallel for  num_threads(3) reduction(+:global_transition_prob)
		for (int it=1; it<rel_length; it++)
		{
		// From new to old
		// mean_A
			global_transition_prob+=log(gsl_ran_lognormal_pdf(previous_internal_rel[it][0], log(internal_rel[it][0])-pow(proposal_sd[it][0],2)/2 ,proposal_sd[it][0]));
		//	global_transition_prob+=log(gsl_ran_lognormal_pdf(previous_internal_rel[it][1], log(internal_rel[it][1])-pow(proposal_sd[it][1],2)/2 ,proposal_sd[it][1]));
		// From old to new
		// mean_A
			global_transition_prob-=log(gsl_ran_lognormal_pdf(internal_rel[it][0], log(previous_internal_rel[it][0])-pow(proposal_sd[it][0],2)/2 ,proposal_sd[it][0]));
		//	global_transition_prob-=log(gsl_ran_lognormal_pdf(internal_rel[it][1], log(previous_internal_rel[it][1])-pow(proposal_sd[it][1],2)/2 ,proposal_sd[it][1]));
		}	

// Accept or reject
	
	
		if (global_current_prob+current_hyperprior_prob-global_previous_prob-previous_hyperprior_prob+global_transition_prob>0)		// New parameter set better => Accept
		{
			previous_rel=new_rel;
			previous_internal_rel=internal_rel;
			global_previous_prob=global_current_prob;
			previous_hyperprior_prob=current_hyperprior_prob;
			without_change=0;
			accepted++;

			for (int stars_it=0; stars_it<star_cat.size(); stars_it++){star_cat[stars_it].last_A_prob=proposed_probs[stars_it];}
		
		//	cout << global_A_chain.size() << " " << global_current_prob << " " << internal_rel[50][0] << " " << new_rel[rel_length-1][0] << " " << current_hyperprior_prob << " " << accepted/global_A_chain.size() << endl;
		}

		else if (exp(global_current_prob+current_hyperprior_prob-global_previous_prob-previous_hyperprior_prob+global_transition_prob)>gsl_ran_flat(rng_handle, 0, 1))	// New set worse => accept with P=P(new)/P(old)
		{
			previous_rel=new_rel;
			previous_internal_rel=internal_rel;
			global_previous_prob=global_current_prob;
			previous_hyperprior_prob=current_hyperprior_prob;
			without_change=0;
			accepted++;

		//	cout << global_A_chain.size() << " " << global_current_prob << " " << internal_rel[50][0] << " " << new_rel[rel_length-1][0] << " " << current_hyperprior_prob << " " << accepted/global_A_chain.size() <<  endl;
		}
		else 
		{
			without_change++;
		//	cout << "fail " << global_current_prob << " " << global_previous_prob << " " << global_transition_prob << " " << current_hyperprior_prob << " " << previous_hyperprior_prob << " " << star_cat.size() << endl;//*/

		}
		if (it_num/1000.==floor(it_num/1000.)){cout << it_num << " " << global_previous_prob << " " << internal_rel[50][0] << " " << previous_rel[rel_length-1][0] << " " << previous_hyperprior_prob << " " << accepted << " " << accepted/it_num << endl;}

		if (floor(it_num/50.)==it_num/50){global_A_chain.push_back(previous_rel);}
		it_num++;
	//	cout << it_num << " " << previous_rel[90][0] << " " << global_current_prob+current_hyperprior_prob <<  endl;
	}
	
	#pragma omp parallel for  num_threads(3)
	for (int star_it=0; star_it<star_cat.size(); star_it++){star_cat[star_it].mean_intervals();}

	if (pow(10,initial_dists[int(initial_dists.size()*0.95)]/5+1)<15000)
	{
		rel_length=150;//floor(pow(10,initial_dists[int(initial_dists.size()*0.95)]/5+1)/100);
		A_mean.resize(rel_length);
	}

	#pragma omp parallel for  num_threads(3)
	for (int it=0; it<rel_length; it++)
	{
		double A_sum=0., sigma_sum=0.;
		for (int m=floor(0.70*global_A_chain.size()); m<global_A_chain.size(); m++)
		{
			A_sum+=global_A_chain[m][it][0];
			sigma_sum+=log(global_A_chain[m][it][1]);
		}
		A_mean[it].mean_A=A_sum/ceil(0.30*global_A_chain.size());
		A_mean[it].sigma=exp(sigma_sum/ceil(0.30*global_A_chain.size()));

		vector <double> A_diffs, sigma_diffs;
		for (int m=floor(0.70*global_A_chain.size()); m<global_A_chain.size(); m++)
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

}


