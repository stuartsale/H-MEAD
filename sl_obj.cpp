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
	batch_accepted=0;
	without_change=0;
	thin=200;
	rel_length=150;

	it_num=0.;

	neighbour_sl=NULL;
}

sl_obj::sl_obj(string filename, float l_in, float b_in, string datatype, float r_max_in, float i_max_in, float ha_max_in)
{
	// Set up variables ---------------------------------------------------------------------
	l=l_in;
	b=b_in;
	unsigned found=filename.find_last_of("/");
	if (found==string::npos){rootname=filename;}
	else{rootname=filename.substr(found+1);}
	rootname.erase(rootname.size()-4);

	// set default MIN vals
	r_min=13.5; 
	i_min=12.0;
	ha_min=12.0;				
	// set default MAX vals
	r_max=r_max_in;
	i_max=i_max_in;
	ha_max=ha_max_in;

	//	global_previous_prob=0;
	//	previous_hyperprior_prob=0;
	//	current_hyperprior_prob=0;

	// Read in data	-------------------------------------------------------------------------

	if (datatype=="iphas" || datatype=="IPHAS"){star_cat=iphas_read(filename,r_min,i_min,ha_min,r_max,i_max,ha_max);}
	else if (datatype=="2MASS" || datatype=="2mass"){star_cat=TWOMASS_read(filename,J_min, H_min, K_min,J_max,H_max,K_max);}
	else {cout << "Unrecognised datatype: " << datatype << endl;}

	float l_min=361;
	float l_max=-1;
	float b_max=-181;
	float b_min=181;

	for (int it=0; it<star_cat.size(); it++)
	{
		if (star_cat[it].l>l_max){l_max=star_cat[it].l;}
		if (star_cat[it].l<l_min){l_min=star_cat[it].l;}
		if (star_cat[it].b>b_max){b_max=star_cat[it].b;}
		if (star_cat[it].b<b_min){b_min=star_cat[it].b;}
	}

	dl=(l_max-l_min)*180./PI;
	db=(b_max-b_min)*180./PI;

	// Find expected A(d) -------------------------------------------------------------------

	//vector<bin_obj2> backup_A_mean (150);
	backup_A_mean.resize(150);
	backup_A_mean=backup_A_mean_find(l, b);

	A_mean.resize(150);

	sigma_fac=0.025;
	accepted=0;
	batch_accepted=0;
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
	dummy_string=hmd_dir+rootname+".td4";
	//cout << "out: " << dummy_string << endl;
	A_out.open(dummy_string.c_str(), ios::trunc);
	//A_out << "#\tdist\tA\tsigma_A\n";

	for (int x=0; x<A_mean.size(); x++)
	{
		A_out << x*100 + 50 << "\t" << A_mean[x].mean_A << "\t" << A_mean[x].sigma <<"\t"<<A_mean[x].d_mean<<"\t"<<A_mean[x].d_sigma<<"\t"
			<<A_mean[x].size<<"\t"<<A_mean[x].error_measure<<"\t"<<A_mean[x].sum<<"\t"<<A_mean[x].diff<<"\t"<<A_mean[x].d_diff<<"\t"<<A_mean[x].median_A<<"\t"<<A_mean[x].median_sigma<<"\n";
	}
	A_out.close();


	ofstream output;
	dummy_string=hmd_dir+rootname+"-090.dat";
	output.open(dummy_string.c_str(), ios::trunc);
	output << "#\tr\ti\tha\tr_i0\tdist\tA\tdistbin\td_A\td_r_i0\td_dist\td_r\td_i\td_ha\tmag_weight\tprob\tA_prob\tMi\tlogAge\tfeh\td_Mi\td_lagAge\td_feh\tlogT\tlogg\trx\tix\thax\tupdate_prop\n" ;
	for (int y=0; y<star_cat.size(); y++)
	{
		
		output << star_cat[y].r << "\t" << star_cat[y].i << "\t" << star_cat[y].ha << "\t"<< star_cat[y].r_i0 << "\t" << star_cat[y].dist << "\t" << star_cat[y].A << "\t" 
			<< star_cat[y].distbin << "\t" << star_cat[y].d_A << "\t" << star_cat[y].d_r_i0 << "\t" << star_cat[y].d_dist << "\t" << star_cat[y].d_r << "\t"
			<< star_cat[y].d_i << "\t" << star_cat[y].d_ha << "\t" << star_cat[y].last_iso.Mi  << "\t" << star_cat[y].mean_prob  << "\t" 
			<< star_cat[y].mean_A_prob  << "\t" << star_cat[y].Mi  << "\t" << star_cat[y].logAge  << "\t" << star_cat[y].feh  << "\t" << star_cat[y].d_Mi  << "\t" 
			<< star_cat[y].d_logAge  << "\t" << star_cat[y].d_feh << "\t" << star_cat[y].logT << "\t" << star_cat[y].logg << "\t" << star_cat[y].rx << "\t" 
			<< star_cat[y].ix << "\t" << star_cat[y].hax << "\t" << star_cat[y].update_prop << "\n";
	}
	output.close();
	
	float sum1=0;
	ofstream rho_out;
	dummy_string=hmd_dir+rootname+".rho";
	rho_out.open(dummy_string.c_str(), ios::trunc);
	rho_out << "#\tdist\trho\td_rho\tA\n" ;
	for (int x=0; x<rho_final.size(); x++)
	{
		sum1+= rho_final[x][0];
		rho_out << x*100 + 50 << "\t" << rho_final[x][0] << "\t" << rho_final[x][1] << "\t" << sum1 <<"\n";
	}
	rho_out.close();


	ofstream samp_out;
	dummy_string=hmd_dir+rootname+".samp";	
	samp_out.open(dummy_string.c_str(), ios::trunc);
	for (int x=0; x<global_A_chain[0][0].size(); x++)
	{
		samp_out << x*100 + 50 << "\t";
		for (int y=0; y<20; y++){samp_out << global_A_chain[int((0.70+y*0.3/20)*global_A_chain.size())][0][x] << "\t";}
		samp_out << endl;
	}
	samp_out.close();

	ofstream srho_out;
	dummy_string=hmd_dir+rootname+".srho";	
	srho_out.open(dummy_string.c_str(), ios::trunc);
	for (int x=0; x<global_A_chain[0][0].size(); x++)
	{
		srho_out << x*100 + 50 << "\t";
		if (x==0){for (int y=0; y<20; y++){srho_out << global_A_chain[int((0.70+y*0.3/20)*global_A_chain.size())][0][x] << "\t";}}
		else {for (int y=0; y<20; y++){srho_out << global_A_chain[int((0.70+y*0.3/20)*global_A_chain.size())][0][x]-
						global_A_chain[int((0.70+y*0.3/20)*global_A_chain.size())][0][x-1]  << "\t";}}
		srho_out << endl;
	}
	srho_out.close();
}


void sl_obj::initial_guess(vector<iso_obj> &isochrones, vector<iso_obj> &guess_set, vector <LF> &LFs)
{
	for (int lf_it=0; lf_it<LFs.size(); lf_it++){LFs[lf_it].set_prior_lf(l,b, dl, db); LFs[lf_it].precompute_Aminmax();}

	global_previous_prob=0;
	previous_hyperprior_prob=0;
	current_hyperprior_prob=0;
	previous_xsl_prob=0;
	current_xsl_prob=0;

	proposal_sd.resize(150, vector <float> (2));
	previous_rel.resize(4, vector <float> (150));
	internal_rel.resize(150, vector <float> (2));
	previous_internal_rel.resize(150, vector <float> (2));
	first_internal_rel.resize(150, vector <float> (2));

// Start from backup_A_mean


	for (int i=0; i<150; i++)
	{

		previous_rel[0][i]=backup_A_mean[i].mean_A;
		previous_rel[1][i]=sqrt(0.5*previous_rel[0][i]);
		previous_rel[3][i]=sqrt(log(1+pow(previous_rel[1][i]/previous_rel[0][i],2)));
		previous_rel[2][i]=log(previous_rel[0][i])-pow(previous_rel[3][i],2)/2;
	}

	previous_internal_rel[0][0]=previous_rel[0][0];//log(previous_rel[0][0]);//
	previous_internal_rel[0][1]=0.5;//previous_rel[1][0]/previous_internal_rel[0][0];

	//previous_hyperprior_prob+=log(gsl_ran_lognormal_pdf(previous_internal_rel[0][0],log(previous_internal_rel[0][0])-0.75,1.5));
//        previous_hyperprior_prob+=log(gsl_ran_lognormal_pdf(previous_internal_rel[0][0],log(previous_internal_rel[0][0]),3.5));
        previous_hyperprior_prob+=log(gsl_ran_lognormal_pdf(previous_internal_rel[0][1],log(0.4)-pow(1.5,2)/2.,1.5));
//	previous_hyperprior_prob+=-log(previous_internal_rel[0][1]);
//	previous_hyperprior_prob+=-log(previous_internal_rel[0][0]);

		

	for (int i=1; i<150; i++)
	{

		previous_internal_rel[i][0]=previous_rel[0][i]-previous_rel[0][i-1];//log(previous_rel[0][i]-previous_rel[0][i-1]);//
		previous_internal_rel[i][1]=0.5;//sqrt(pow(previous_rel[1][i],2)-pow(previous_rel[1][i-1],2))/previous_internal_rel[0][i];
	


//                previous_hyperprior_prob+=log(gsl_ran_lognormal_pdf(previous_internal_rel[i][0],log(previous_internal_rel[i][0]),3.5));
                previous_hyperprior_prob+=log(gsl_ran_lognormal_pdf(previous_internal_rel[i][1],log(0.4)-pow(1.5,2)/2.,1.5));
//		previous_hyperprior_prob+=-log(previous_internal_rel[i][1]);
//		previous_hyperprior_prob+=-log(previous_internal_rel[i][0]);
	}

	first_internal_rel=previous_internal_rel;
	first_rel=previous_rel;

	global_A_chain.push_back(previous_rel);

// Dump sources not in required region of c-c

	int it_stars=0;
	while (it_stars<star_cat.size())
	{
		if (star_cat[it_stars].r-star_cat[it_stars].ha>guess_set[0].redline(star_cat[it_stars].r-star_cat[it_stars].i) || star_cat[it_stars].r-star_cat[it_stars].ha<guess_set[guess_set.size()-1].redline(star_cat[it_stars].r-star_cat[it_stars].i)-0.05)
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

	for (int it=0; it<10000; it++)
	{
		for (int it=0; it<star_cat.size(); it++)
		{
			star_cat[it].star_try1(isochrones, l, b, previous_rel);
		}
	}

	previous_norm_prob=0;
	for (int it_LF=0; it_LF<LFs.size(); it_LF++)
	{
		previous_norm_prob+=-log(LFs[it_LF].LF_prob2(previous_rel)+LFs[it_LF].beta)*(star_cat.size()+LFs[it_LF].alpha);
	}

	for (int it=0; it<rel_length; it++)
	{
		proposal_sd[it][0]=2*sigma_fac/(star_cat.size()/200.);
		proposal_sd[it][1]=sigma_fac/2/(star_cat.size()/200.);
	}

}



void sl_obj::dist_redMCMC(vector<iso_obj> &isochrones, vector <LF> &LFs)
{
	while (it_num<150000 )
	{
		update(isochrones, LFs);
	}
}

void sl_obj::update(vector<iso_obj> &isochrones, vector <LF> &LFs)
{
		global_current_prob=0;
		global_transition_prob=0;
		global_previous_prob=0;
		current_hyperprior_prob=0;
		last_part_prior1=0;
		test_part_prior1=0;
		last_part_prior2=0;
		test_part_prior2=0;

// First vary parameters for each star

		float dummy=0;
		int max_dist_bin1=0;
		int max_dist_bin2=0;
		vector <int> dist_bin_vec;
//		#pragma omp parallel for  num_threads(3) reduction(+:dummy)
		for (int it=0; it<star_cat.size(); it++)
		{
			/*if (gsl_ran_flat(rng_handle, 0, 1)>0.){*/star_cat[it].star_try1(isochrones, l, b, previous_rel);//};
			dummy+=star_cat[it].last_A_prob;
//			if (int(pow(10,star_cat[it].last_dist_mod/5+1)/100.)>max_dist_bin){max_dist_bin=int(pow(10,star_cat[it].last_dist_mod/5+1)/100.);}
			dist_bin_vec.push_back(int(pow(10,star_cat[it].last_dist_mod/5+1)/100.));
		}
		global_previous_prob=dummy;
//		max_dist_bin=min(150, max_dist_bin);
		sort(dist_bin_vec.begin(), dist_bin_vec.end());
//		cout << dist_bin_vec.back() << " " <<dist_bin_vec.back()*0.75 << " " << dist_bin_vec[int(0.9*dist_bin_vec.size())] << " " << dist_bin_vec[int(0.9*dist_bin_vec.size())] << " " << dist_bin_vec[int(0.75*dist_bin_vec.size())]  << endl;

		max_dist_bin1=min(dist_bin_vec[int(0.9*dist_bin_vec.size())],150);
		max_dist_bin2=min(dist_bin_vec[int(0.95*dist_bin_vec.size())],150);

// Now vary hyper-parameters

		vector < vector <float> > new_rel(4,vector <float> (rel_length));

		for (int it=0; it<max_dist_bin1; it++)
		{	
			internal_rel[it][0]=gsl_ran_lognormal(rng_handle,log(previous_internal_rel[it][0])-pow(proposal_sd[it][0],2)/2,proposal_sd[it][0]);
			internal_rel[it][1]=gsl_ran_lognormal(rng_handle,log(previous_internal_rel[it][1])-pow(proposal_sd[it][1],2)/2,proposal_sd[it][1]);

			
//                        test_part_prior1+=log(gsl_ran_lognormal_pdf(internal_rel[it][0],log(first_internal_rel[it][0])-0.125,.5));
//                        last_part_prior1+=log(gsl_ran_lognormal_pdf(previous_internal_rel[it][0],log(first_internal_rel[it][0])-0.125,.5));
                        current_hyperprior_prob+=log(gsl_ran_lognormal_pdf(internal_rel[it][1],log(0.4)-pow(1.5,2)/2.,1.5));
		//	current_hyperprior_prob+=-log(internal_rel[it][1]);
//			current_hyperprior_prob+=-log(internal_rel[it][0]);
		}
//		last_part_prior1*=0*rel_length/(max_dist_bin1);
//		test_part_prior1*=0*rel_length/(max_dist_bin1);

		new_rel[0][0]=internal_rel[0][0];
		new_rel[1][0]=sqrt(internal_rel[0][1]*internal_rel[0][0]);//internal_rel[0][1];//*internal_rel[0][0];
		new_rel[3][0]=sqrt(log(1+pow(new_rel[1][0]/new_rel[0][0],2)));
		new_rel[2][0]=log(new_rel[0][0])-pow(new_rel[3][0],2)/2;
		
		for (int it=1; it<max_dist_bin1; it++)
		{			
			new_rel[0][it]=internal_rel[it][0]+new_rel[0][it-1];
			new_rel[1][it]=sqrt(pow(new_rel[1][it-1],2)+internal_rel[it][0]*internal_rel[it][1]);
			new_rel[3][it]=sqrt(log(1+pow(new_rel[1][it]/new_rel[0][it],2)));
			new_rel[2][it]=log(new_rel[0][it])-pow(new_rel[3][it],2)/2;

		}	
		
		float ratio_fac1;
		float ratio_fac2;
		for (int it=max_dist_bin1; it<rel_length; it++)
		{	
			internal_rel[it][0]=gsl_ran_lognormal(rng_handle,log(previous_internal_rel[it][0])-pow(proposal_sd[it][0],2)/2,proposal_sd[it][0]);
			internal_rel[it][1]=gsl_ran_lognormal(rng_handle,log(previous_internal_rel[it][1])-pow(proposal_sd[it][1],2)/2,proposal_sd[it][1]);

			ratio_fac1=(new_rel[0][max_dist_bin1-1]-new_rel[0][max(max_dist_bin1-20,0)])/(first_rel[0][max_dist_bin1-1]-first_rel[0][max(max_dist_bin1-20,0)]);
			ratio_fac2=(previous_rel[0][max_dist_bin1-1]-previous_rel[0][max(max_dist_bin1-20,0)])/(first_rel[0][max_dist_bin1-1]-first_rel[0][max(max_dist_bin1-20,0)]);
//                        test_part_prior2+=log(gsl_ran_lognormal_pdf(internal_rel[it][0],log(first_internal_rel[it][0]*ratio_fac1)-0.125/ratio_fac1,.5/ratio_fac1));
//                        last_part_prior2+=log(gsl_ran_lognormal_pdf(previous_internal_rel[it][0],log(first_internal_rel[it][0]*ratio_fac2)-0.125/ratio_fac2,.5/ratio_fac2));
			test_part_prior2+=-pow(internal_rel[it][0]-first_internal_rel[it][0]*ratio_fac1,2)/(2.5*first_internal_rel[it][0]+1E-9);
			last_part_prior2+=-pow(previous_internal_rel[it][0]-first_internal_rel[it][0]*ratio_fac2,2)/(2.5*first_internal_rel[it][0]+1E-9);
                        current_hyperprior_prob+=log(gsl_ran_lognormal_pdf(internal_rel[it][1],log(0.4)-pow(1.5,2)/2.,1.5));
	//		current_hyperprior_prob+=-log(internal_rel[it][1]);
//			current_hyperprior_prob+=-log(internal_rel[it][0]);
		}

		last_part_prior2*=(star_cat.size()*1.)*rel_length/max(rel_length-max_dist_bin1,1);
		test_part_prior2*=(star_cat.size()*1.)*rel_length/max(rel_length-max_dist_bin1,1);


	
		for (int it=max_dist_bin1; it<rel_length; it++)
		{			
			new_rel[0][it]=internal_rel[it][0]+new_rel[0][it-1];
			new_rel[1][it]=sqrt(pow(new_rel[1][it-1],2)+internal_rel[it][0]*internal_rel[it][1]);
			new_rel[3][it]=sqrt(log(1+pow(new_rel[1][it]/new_rel[0][it],2)));
			new_rel[2][it]=log(new_rel[0][it])-pow(new_rel[3][it],2)/2;
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
			current_norm_prob+=-log(LFs[it_LF].LF_prob2(new_rel)+LFs[it_LF].beta)*(star_cat.size()+LFs[it_LF].alpha);
		}

// Metropolis-Hastings algorithm step

		global_transition_prob=0;

		dummy=0;
//		#pragma omp parallel for  num_threads(3) reduction(+:dummy)
		for (int it=0; it<150; it++)
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
	
	
		if (global_current_prob+current_hyperprior_prob-global_previous_prob-previous_hyperprior_prob+global_transition_prob+current_xsl_prob-previous_xsl_prob
			 + current_norm_prob-previous_norm_prob+test_part_prior1-last_part_prior1+test_part_prior2-last_part_prior2>0)		// New parameter set better => Accept
		{
			previous_rel=new_rel;
			previous_internal_rel=internal_rel;
			global_previous_prob=global_current_prob;
			previous_hyperprior_prob=current_hyperprior_prob;
			previous_xsl_prob=current_xsl_prob;
			previous_norm_prob=current_norm_prob;
			without_change=0;
			accepted++;
			batch_accepted++;

			for (int stars_it=0; stars_it<star_cat.size(); stars_it++){star_cat[stars_it].last_A_prob=proposed_probs[stars_it];}
		
		//	cout << global_A_chain.size() << " " << global_current_prob << " " << internal_rel[50][0] << " " << new_rel[rel_length-1][0] << " " << current_hyperprior_prob << " " << accepted/global_A_chain.size() << endl;
		}

		else if (exp(global_current_prob+current_hyperprior_prob-global_previous_prob-previous_hyperprior_prob+global_transition_prob+current_xsl_prob-previous_xsl_prob
			+ current_norm_prob-previous_norm_prob+test_part_prior1-last_part_prior1+test_part_prior2-last_part_prior2)>gsl_ran_flat(rng_handle, 0, 1))	// New set worse => accept with P=P(new)/P(old)
		{
			previous_rel=new_rel;
			previous_internal_rel=internal_rel;
			global_previous_prob=global_current_prob;
			previous_hyperprior_prob=current_hyperprior_prob;
			previous_xsl_prob=current_xsl_prob;
			previous_norm_prob=current_norm_prob;
			without_change=0;
			accepted++;
			batch_accepted++;

		//	cout << global_A_chain.size() << " " << global_current_prob << " " << internal_rel[50][0] << " " << new_rel[rel_length-1][0] << " " << current_hyperprior_prob << " " << accepted/global_A_chain.size() <<  endl;
		}
		else 
		{
			without_change++;
//			cout << "fail " << (global_current_prob+current_hyperprior_prob-global_previous_prob-previous_hyperprior_prob+global_transition_prob+current_xsl_prob-previous_xsl_prob
//			+ current_norm_prob-previous_norm_prob+test_part_prior2-last_part_prior2) << " " <<
// global_current_prob << " " << global_previous_prob << " " << global_transition_prob << " " << current_hyperprior_prob << " " << previous_hyperprior_prob << " " << star_cat.size() << " " << new_rel[1][50] << " " << previous_rel[1][50] << endl;//*/

		}
		if (neighbour_sl){if (it_num/1000.==floor(it_num/1000.)){cout << it_num << " " << global_previous_prob << " " << internal_rel[50][0] << " " << previous_rel[0][rel_length-1] << " " << previous_hyperprior_prob << " " << accepted << " " << accepted/it_num << " " << neighbour_sl->previous_rel[0][rel_length-1] << endl;}}
		else {if (it_num/149999.==floor(it_num/149999.)){cout << it_num << " " << global_previous_prob << " " << internal_rel[50][0] << " " << previous_rel[0][rel_length-1] << " " << previous_hyperprior_prob << " " << previous_norm_prob << " " << max_dist_bin1 << " " << max_dist_bin2 << " " << last_part_prior2 << " " << accepted << " " << accepted/it_num << endl;}}

		if (floor(it_num/100.)==it_num/100){global_A_chain.push_back(previous_rel);}
		it_num++;


	// Adaptive step

	if (it_num/1000.==floor(it_num/1000.) && it_num<75000)
	{
		if (batch_accepted/1000.>0.234)
		{
			for (int it=0; it<rel_length; it++)
			{
				proposal_sd[it][0]*=1.1;
				proposal_sd[it][1]*=1.1;
			}
		}
		else
		{
			for (int it=0; it<rel_length; it++)
			{
				proposal_sd[it][0]/=1.1;
				proposal_sd[it][1]/=1.1;
			}
		}
	batch_accepted=0;

	for (int it=0; it<star_cat.size(); it++){star_cat[it].adaptive_proposal_update(1000);}
	}


}


void sl_obj::mean_intervals(void)
{
//	#pragma omp parallel for  num_threads(3)
	for (int star_it=0; star_it<star_cat.size(); star_it++){star_cat[star_it].mean_intervals();}

	//#pragma omp parallel for  num_threads(3)
	for (int it=0; it<150; it++)
	{
		vector <float> unburnt_A, unburnt_sigma;

		float A_sum=0., sigma_sum=0.;
		for (int m=floor(0.70*global_A_chain.size()); m<global_A_chain.size(); m++)
		{
			A_sum+=global_A_chain[m][0][it];
			sigma_sum+=log(global_A_chain[m][1][it]);
			unburnt_A.push_back(global_A_chain[m][0][it]);
			unburnt_sigma.push_back(global_A_chain[m][1][it]);
		}
		A_mean[it].mean_A=A_sum/ceil(0.30*global_A_chain.size());
		A_mean[it].sigma=exp(sigma_sum/ceil(0.30*global_A_chain.size()));

		sort(unburnt_A.begin(), unburnt_A.end());
		sort(unburnt_sigma.begin(), unburnt_sigma.end());
		A_mean[it].median_A=unburnt_A[int(0.5*unburnt_A.size())];
		A_mean[it].median_sigma=unburnt_sigma[int(0.5*unburnt_sigma.size())];

		vector <float> A_diffs, sigma_diffs;
		for (int m=floor(0.70*global_A_chain.size()); m<global_A_chain.size(); m++)
		{
			A_diffs.push_back(abs(global_A_chain[m][0][it]-A_mean[it].mean_A));
			sigma_diffs.push_back(abs(global_A_chain[m][1][it]-A_mean[it].sigma));
		}
		sort(A_diffs.begin(),A_diffs.end());
		sort(sigma_diffs.begin(), sigma_diffs.end());

		A_mean[it].d_mean=A_diffs[int(0.682*A_diffs.size())];
		A_mean[it].d_sigma=sigma_diffs[int(0.682*A_diffs.size())];

		A_mean[it].error_measure=sqrt(pow(A_mean[it].d_mean/A_mean[it].mean_A,2)+pow(A_mean[it].d_sigma/A_mean[it].sigma,2));
	}
	
	rho_final.resize(150);
	rho_final[0].resize(2);
	vector <float> rho_diffs;

	float rho_sum=0.;
	for (int m=floor(0.70*global_A_chain.size()); m<global_A_chain.size(); m++)
	{
		rho_sum+=global_A_chain[m][0][0];
	}
	rho_final[0][0]=rho_sum/ceil(0.30*global_A_chain.size());
	for (int m=floor(0.70*global_A_chain.size()); m<global_A_chain.size(); m++)
	{
		rho_diffs.push_back(abs(global_A_chain[m][0][0]-rho_final[0][0]) );
	}
	sort(rho_diffs.begin(),rho_diffs.end());
	rho_final[0][1]=rho_diffs[int(0.682*rho_diffs.size())];

	for (int it=1; it<150; it++)
	{
		rho_final[it].resize(2);
		vector <float> rho_diffs;

		float rho_sum=0.;
		for (int m=floor(0.70*global_A_chain.size()); m<global_A_chain.size(); m++)
		{
			rho_sum+=global_A_chain[m][0][it]-global_A_chain[m][0][it-1];
		}
		rho_final[it][0]=rho_sum/ceil(0.30*global_A_chain.size());
		for (int m=floor(0.70*global_A_chain.size()); m<global_A_chain.size(); m++)
		{
			rho_diffs.push_back(abs(global_A_chain[m][0][it]-global_A_chain[m][0][it-1]-rho_final[it][0]) );
		}
		sort(rho_diffs.begin(),rho_diffs.end());
		rho_final[it][1]=rho_diffs[int(0.682*rho_diffs.size())];
	}
	

}

void sl_obj::neighbour_set(sl_obj * neighbour)
{
	neighbour_sl=neighbour;
}


void sl_obj::acl_calc(void)
{
	ofstream acl_out;
	string dummy_string;
	dummy_string=hmd_dir+rootname+".acl";
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


