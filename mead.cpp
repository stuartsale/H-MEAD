#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdlib> 
#include <sstream>
#include <ctime>
#include "helper.h"
#include "bin_obj.h"
#include "iso_obj.h"
#include "iphas_obj.h"
#include "omp.h"
using namespace std;

#ifndef PI
#define PI 3.14159265358979323846
#endif
#ifndef log2PIover2
#define log2PIover2 0.918938533
#endif

//Reads in IPHAS datas
//vector<iphas_obj> iphas_read(string filename); 
//vector<iphas_obj> iphas_read(string filename,double &i_min,double &r_min,double &ha_min);		
vector<iphas_obj> iphas_read(string filename,double &r_min,double &i_min,double &ha_min,double &r_max, double &i_max, double &ha_max);		
vector<iso_obj> iso_read(const string &filename);
vector<iso_obj> iso_read_Tg(const string &filename);
vector<bin_obj2> backup_A_mean_find(double l_gal, double b_gal);		

// Gets stdout of a system cmd
// Used to read in Schlegel Galactic reddening
// Code 'borrowed' from somewhere on the net
//
// (BM notes: replace this with CFITSIO access to Schlegel ?)
string getStdoutFromCommand(string cmd);

// Dump me baby, yeah!
void output_write(string filename);

//converts double => C++ string
string stringify(double x);	


double r_min, i_min, ha_min, r_max, i_max, ha_max;
vector <vector <vector <double> > > lookup_table;

gsl_rng* rng_handle;

//--------------------------------
// MAIN 
//--------------------------------

//Mead::Mead(string ifname){

int main(int argc, char* argv[]) 
{	

	//initialize up random number generator
    	//gsl_rng_env_setup ();
	 rng_handle = gsl_rng_alloc (gsl_rng_taus2);

	//seed the random no generator
	gsl_rng_set(rng_handle, time(0));

	vector<iphas_obj> colours;		
	vector<bin_obj2> A_mean;
	vector<bin_obj2> backup_A_mean (150);				
	A_mean.reserve(251); 
//	backup_A_mean.reserve(251); 


	//initialize up random number generator
    	gsl_rng_env_setup ();
	gsl_rng* rng_handle = gsl_rng_alloc (gsl_rng_taus2);

	//seed the random no generator
	//gsl_rng_set(rng_handle, time(0));

	lookup_table=lookup_creator();

// set default MIN vals
	r_min=13.5; 
	i_min=12.0;
	ha_min=12.0;				
// set default MAX vals
	r_max=0.;
	i_max=0.;
	ha_max=0.;

      //Reading in isochrones data into a vector

	//vector<iso_obj> isochrones=iso_read("padova-iso_reg.dat");
	vector<iso_obj> isochrones=iso_read_Tg("padova-iso_tefflogg3.dat");

	vector<iso_obj> guess_set;
//	guess_set.push_back(iso_get_Tg(0.,3.574 ,5.00 , isochrones));	//M1
//	guess_set.push_back(iso_get_Tg(0.,3.591 ,4.95 , isochrones));	//M0
//	guess_set.push_back(iso_get_Tg(0.,3.602 ,4.90 , isochrones));	//K7
//	guess_set.push_back(iso_get_Tg(0.,3.643 ,4.65 , isochrones));	//K5
	guess_set.push_back(iso_get_Tg(0.,3.663 ,4.57 , isochrones));	//K4
	guess_set.push_back(iso_get_Tg(0.,3.672 ,4.56 , isochrones));	//K3
	guess_set.push_back(iso_get_Tg(0.,3.686 ,4.55 , isochrones));	//K2
	guess_set.push_back(iso_get_Tg(0.,3.695 ,4.55 , isochrones));	//K1
	guess_set.push_back(iso_get_Tg(0.,3.703 ,4.57 , isochrones));	//K0
	guess_set.push_back(iso_get_Tg(0.,3.720 ,4.55 , isochrones));	//G8
	guess_set.push_back(iso_get_Tg(0.,3.740 ,4.49 , isochrones));	//G5
	guess_set.push_back(iso_get_Tg(0.,3.763 ,4.40 , isochrones));	//G2
	guess_set.push_back(iso_get_Tg(0.,3.774 ,4.39 , isochrones));	//G0
	guess_set.push_back(iso_get_Tg(0.,3.789 ,4.35 , isochrones));	//F8
	guess_set.push_back(iso_get_Tg(0.,3.813 ,4.34 , isochrones));	//F5
	guess_set.push_back(iso_get_Tg(0.,3.845 ,4.26 , isochrones));	//F2
	guess_set.push_back(iso_get_Tg(0.,3.863 ,4.28 , isochrones));	//F0
	guess_set.push_back(iso_get_Tg(0.,3.903 ,4.26 , isochrones));	//A7
	guess_set.push_back(iso_get_Tg(0.,3.924 ,4.22 , isochrones));	//A5
	guess_set.push_back(iso_get_Tg(0.,3.949 ,4.20 , isochrones));	//A3
	guess_set.push_back(iso_get_Tg(0.,3.961 ,4.16 , isochrones));	//A2


//	while (1.0648*guess_set[guess_set.size()-1].Mi<2.060){guess_set.push_back(iso_get(0., 1.0648*guess_set[guess_set.size()-1].Mi, 8.5, isochrones));}
//	guess_set.push_back(iso_get(0., 2.060, 8.5, isochrones));

   // Read in IPHAS data
	string iphas_filename=argv[1];		
	colours=iphas_read(iphas_filename,r_min,i_min,ha_min,r_max,i_max,ha_max);
	//r_max-=0.2;
	//i_max-=0.2;
	//ha_max-=0.2;
	r_max=21.5;
	i_max=20.5;
	ha_max=20.5;
	//cout << "r_min=" << r_min << " i_min=" << i_min << " ha_min=" << ha_min << endl; 
	//cout << "r_max=" << r_max << " i_max=" << i_max << " ha_max=" << ha_max << endl; 


   //
   // do this after iphas_read so that limiting mags are accurate.

	backup_A_mean=backup_A_mean_find(atof(argv[2]), atof(argv[3]));

//	iphas_obj test_obj(14.821,14.462,14.581, 0.002,0.002,0.002,180,0);
//	test_obj.dist_redMCMC(isochrones, guess_set 180.,0.);

  
   // governs how aggresively objects near mag limits 
   // are downweighted - see section 4.6 in my thesis
	//double x_value=2.0;	
			
// first run through


//	cout << "backup size:" << backup_A_mean.size() << endl;

//	cout << "Real: " << real_prob(colours, isochrones, guess_set, atof(argv[2]), atof(argv[3]), backup_A_mean, -0.0272, 0.53) << endl;


	clock_t start;
	start=time(NULL);

	A_mean=dist_redMCMC(colours, isochrones, guess_set, atof(argv[2]), atof(argv[3]), backup_A_mean, -0.0272, 0.53);			// }

	//cout << "total time: " << (time(NULL)-start) <<"s\n";
   	
// Write results to file
	output_write(iphas_filename, A_mean, colours);

	return 0;
}



//-------------------------------------------------------------------
//-------------------------------------------------------------------
//       HELPER FUNCTIONS
//-------------------------------------------------------------------
//-------------------------------------------------------------------

/*double real_prob(vector<iphas_obj> &stars, vector<iso_obj> &isochrones, vector<iso_obj> &guess_set, double l, double b, vector <bin_obj2> backup_A_mean, double ri_min, double ri_max)
{
	vector < vector <double> > previous_rel (150, vector <double> (4));

	for (int i=0; i<150; i++)
	{
		previous_rel[i][0]=backup_A_mean[i].mean_A;
		previous_rel[i][1]=0.2;//backup_A_mean[i].sigma;
		previous_rel[i][3]=sqrt(log(1+pow(previous_rel[i][1]/previous_rel[i][0],2)));
		previous_rel[i][2]=log(previous_rel[i][0])-pow(previous_rel[i][3],2)/2;
	}

// Dump sources not in required region of c-c

	int it_stars=0;
	while (it_stars<stars.size())
	{
		if (stars[it_stars].r-stars[it_stars].ha>guess_set[0].redline(stars[it_stars].r-stars[it_stars].i) || stars[it_stars].r-stars[it_stars].ha<guess_set[guess_set.size()-1].redline(stars[it_stars].r-stars[it_stars].i))
		{
			stars.erase(stars.begin()+it_stars);
		}
		else {it_stars++;}
	}
	
// Make initial guess

	it_stars=0;
	while (it_stars<stars.size())
	{
		stars[it_stars].initial_guess(isochrones, guess_set, previous_rel);
		if (stars[it_stars].last_A>0 && stars[it_stars].last_iso.r0-stars[it_stars].last_iso.i0> ri_min  && stars[it_stars].last_iso.r0-stars[it_stars].last_iso.i0<ri_max){it_stars++;}
		else {stars.erase(stars.begin()+it_stars);}
	}

	cout << "start-------------------\n\n"<<endl;
	double global_previous_prob=0;
	for (int star_it=0; star_it<stars.size(); star_it++)
	{
		iso_obj test_iso;
		test_iso=iso_get(stars[star_it].real_feh, stars[star_it].real_Mi, log10(stars[star_it].real_logAge)+6, isochrones);

		stars[star_it].last_prob=stars[star_it].prob_eval(test_iso, stars[star_it].real_A, 5*log10(stars[star_it].real_dist/10), previous_rel);
		global_previous_prob+=stars[star_it].last_prob;
//		cout << stars[star_it].last_prob << " " <<stars[star_it].real_Mi << " " << stars[star_it].real_logAge << " " << stars[star_it].real_A << " " << stars[star_it].real_dist << " " << stars[star_it].r << " " << stars[star_it].i << " " << stars[star_it].ha << " " << test_iso.r0 << " " << test_iso.i0 << " " << test_iso.ha0 << " " << log(test_iso.IMF()) << endl;
	}

	return global_previous_prob;

}*/

vector <bin_obj2> dist_redMCMC(vector<iphas_obj> &stars, vector<iso_obj> &isochrones, vector<iso_obj> &guess_set, double l, double b, vector <bin_obj2> backup_A_mean, double ri_min, double ri_max)
{



	vector <bin_obj2> first_bins(150);

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
	while (it_stars<stars.size())
	{
		if (stars[it_stars].r-stars[it_stars].ha>guess_set[0].redline(stars[it_stars].r-stars[it_stars].i) || stars[it_stars].r-stars[it_stars].ha<guess_set[guess_set.size()-1].redline(stars[it_stars].r-stars[it_stars].i))
		{
			stars.erase(stars.begin()+it_stars);
		}
		else {it_stars++;}
	}

// Make initial guess

	it_stars=0;
	while (it_stars<stars.size())
	{
		stars[it_stars].initial_guess(isochrones, guess_set, previous_rel);
		if (stars[it_stars].last_A>0 && stars[it_stars].last_iso.r0-stars[it_stars].last_iso.i0> ri_min  && stars[it_stars].last_iso.r0-stars[it_stars].last_iso.i0<ri_max){initial_dists.push_back(stars[it_stars].last_dist_mod);}
		else if (stars[it_stars].last_iso.r0-stars[it_stars].last_iso.i0> ri_min  && stars[it_stars].last_iso.r0-stars[it_stars].last_iso.i0<ri_max){initial_dists.push_back(stars[it_stars].last_dist_mod);stars[it_stars].last_A = 0.02;} 
		it_stars++;
	}

	sort(initial_dists.begin(), initial_dists.end());

//	cout << "95th percentile on dist is: " << pow(10,initial_dists[int(initial_dists.size()*0.95)]/5+1) << " kpc away" << endl;
//	cout << "5th percentile on dist is: " << pow(10,initial_dists[int(initial_dists.size()*0.05)]/5+1) << " kpc away" << endl;

	for (int star_it=0; star_it<stars.size(); star_it++){global_previous_prob+=stars[star_it].last_prob;}
	
// Loop through this section

	vector <double> proposed_probs(stars.size());

	float it_num=0.;

	while (it_num<150000 )
	{
		global_current_prob=0;
		global_transition_prob=0;
		global_previous_prob=0;
		current_hyperprior_prob=0;

// First vary parameters for each star

		#pragma omp parallel for  num_threads(3) reduction(+:global_previous_prob)
		for (int it=0; it<stars.size(); it++)
		{
			/*if (gsl_ran_flat(rng_handle, 0, 1)>0.){*/stars[it].star_try1(isochrones, l, b, previous_rel);//};
			global_previous_prob+=stars[it].last_prob;
		}

	//	cout << stars[161].last_A << " " << stars[161].last_dist_mod << " " << stars[161].last_prob << " " << stars[161].last_iso.logT << " " << stars[161].last_iso.logg << " " << stars[161].last_iso.r0-stars[161].last_iso.i0 << " " << stars[161].last_iso.r0-stars[161].last_iso.ha0 << " " << it_num  << " " << stars[161].last_iso.Mi << " " << log(stars[161].last_iso.Jac) << " " << log(stars[161].last_iso.IMF())  << " " 
//<< stars[161].last_iso.r0+stars[161].last_iso.u*pow(stars[161].last_A,2)+stars[161].last_iso.v*stars[161].last_A+stars[161].last_iso.w+stars[161].last_dist_mod  << " " 
//<< stars[161].last_iso.i0+stars[161].last_iso.u_i*pow(stars[161].last_A,2)+stars[161].last_iso.v_i*stars[161].last_A+stars[161].last_iso.w_i+stars[161].last_dist_mod  << " " << stars[161].last_iso.ha0+stars[161].last_iso.u_ha*pow(stars[161].last_A,2)+stars[161].last_iso.v_ha*stars[161].last_A+stars[161].last_iso.w_ha+stars[161].last_dist_mod << endl;

// Now vary hyper-parameters

		vector < vector <double> > new_rel(rel_length,vector <double> (4));

		for (int it=0; it<rel_length; it++)
		{
			proposal_sd[it][0]=sigma_fac/10;
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
		for (int it=0; it<stars.size(); it++)
		{
			proposed_probs[it]=stars[it].prob_eval(stars[it].last_iso, stars[it].last_A, stars[it].last_dist_mod, new_rel);
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

			for (int stars_it=0; stars_it<stars.size(); stars_it++){stars[stars_it].last_prob=proposed_probs[stars_it];}
		
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
		//	cout << "fail " << global_current_prob << " " << global_previous_prob << " " << global_transition_prob << " " << current_hyperprior_prob << " " << previous_hyperprior_prob << " " << stars.size() << endl;//*/

		}
		if (it_num/1000.==floor(it_num/1000.)){cout << it_num << " " << global_previous_prob << " " << internal_rel[50][0] << " " << previous_rel[rel_length-1][0] << " " << previous_hyperprior_prob << " " << accepted << " " << accepted/it_num << endl;}

		if (floor(it_num/50.)==it_num/50){global_A_chain.push_back(previous_rel);}
		it_num++;
	//	cout << it_num << " " << previous_rel[90][0] << " " << global_current_prob+current_hyperprior_prob <<  endl;
	}
	
	#pragma omp parallel for  num_threads(3)
	for (int star_it=0; star_it<stars.size(); star_it++){stars[star_it].mean_intervals();}

	if (pow(10,initial_dists[int(initial_dists.size()*0.95)]/5+1)<15000)
	{
		rel_length=150;//floor(pow(10,initial_dists[int(initial_dists.size()*0.95)]/5+1)/100);
		first_bins.resize(rel_length);
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
		first_bins[it].mean_A=A_sum/ceil(0.30*global_A_chain.size());
		first_bins[it].sigma=exp(sigma_sum/ceil(0.30*global_A_chain.size()));

		vector <double> A_diffs, sigma_diffs;
		for (int m=floor(0.70*global_A_chain.size()); m<global_A_chain.size(); m++)
		{
			A_diffs.push_back(abs(global_A_chain[m][it][0]-first_bins[it].mean_A));
			sigma_diffs.push_back(abs(global_A_chain[m][it][1]-first_bins[it].sigma));
		}
		sort(A_diffs.begin(),A_diffs.end());
		sort(sigma_diffs.begin(), sigma_diffs.end());

		first_bins[it].d_mean=A_diffs[int(0.682*A_diffs.size())];
		first_bins[it].d_sigma=sigma_diffs[int(0.682*A_diffs.size())];

		first_bins[it].error_measure=sqrt(pow(first_bins[it].d_mean/first_bins[it].mean_A,2)+pow(first_bins[it].d_sigma/first_bins[it].sigma,2));
	}

	//vector<double> sizes (150,0);
	//for (int it=0; it<stars.size(); it++)
	//{	
	//	if (stars[it].distbin<rel_length)
	//	{
	//		first_bins[stars[it].distbin].size++;
//			first_bins[stars[it].distbin].sum+=stars[it].A/pow(stars[it].d_A,2);
	//	}
	//}
	return first_bins;
	
}




string stringify(double x)		
{  
	ostringstream o;
	if(!(o << x)){throw 66;}
	return o.str();
}


vector<iso_obj> iso_read(const string &filename)		// Function to read in calibration data
{
	ifstream input1;
	input1.open(filename.c_str());
	if(!input1) { //output file couldn't be opened
		cerr << "Error: file could not be opened" << endl;
		exit(1);
	}

	vector<iso_obj> totalfile;		// making a vector to store the calibration_obj
	
	while (!input1.eof())				// Running down file - reading it in
	{
		string str; 	
		getline(input1, str);			// grab a line
		string temp;
		stringstream sstest(str);
		sstest>>temp;
		if (temp!="#")				// check the line isn't commented out
		{		
			double buf;
			stringstream ss(str);		// turn that line into a stringstream
		
			vector<double> fromfile;	//vector to put contents of line into
		
			while (ss>>buf){			// Includes implicit conversion from string to double
				fromfile.push_back(buf);	
			}
			if (fromfile.size()==8)		// check there's something in the line
			{
         			iso_obj objnew(fromfile[0], fromfile[1], fromfile[2], fromfile[3], fromfile[4], fromfile[5], fromfile[6], fromfile[7],1);
				totalfile.push_back(objnew);
			}		
		}
	}
	return totalfile;
	input1.close();
}

vector<iso_obj> iso_read_Tg(const string &filename)		// Function to read in calibration data
{
	ifstream input1;
	input1.open(filename.c_str());
	if(!input1) { //output file couldn't be opened
		cerr << "Error: file could not be opened" << endl;
		exit(1);
	}

	vector<iso_obj> totalfile;		// making a vector to store the calibration_obj
	
	while (!input1.eof())				// Running down file - reading it in
	{
		double jac;
		string str; 	
		getline(input1, str);			// grab a line
		string temp;
		stringstream sstest(str);
		sstest>>temp;
		if (temp!="#")				// check the line isn't commented out
		{		
			double buf;
			stringstream ss(str);		// turn that line into a stringstream
		
			vector<double> fromfile;	//vector to put contents of line into
		
			while (ss>>buf){			// Includes implicit conversion from string to double
				fromfile.push_back(buf);	
			}
			if (fromfile.size()==14)		// check there's something in the line
			{
				if (fromfile[5]!=0 && fromfile[6]!=0 && fromfile[7]!=0 && fromfile[8]!=0)
				{
					jac=abs(fromfile[5]*fromfile[8]-fromfile[6]*fromfile[7]);
				}
				else {jac=0;}
         			iso_obj objnew(fromfile[0], fromfile[1], fromfile[2], fromfile[3], fromfile[4], fromfile[9], fromfile[10], fromfile[11], jac);
				totalfile.push_back(objnew);
			}		
		}
	}
	return totalfile;
	input1.close();
}

// function to read in IPHAS data, works in much the same manner as calibration_read
vector<iphas_obj> iphas_read(string filename,double &r_min1,double &i_min1,double &ha_min1,double &r_max1, double &i_max1, double &ha_max1)		
{						
   //-----------------------
   // Pre-cond: filename is a string containing the filename of the IPHAS 
   // catalogue data for the region examined.
   //--------------------
   // Post-cond: iphas_colours is returned.
   // As a vector of iphas_obj, it
   // contains all sources in the IPHAS catalogue
   // iff they match these conditions:
   // (a) Classified as stellar or _probably_ stellar
   //     in all three bands (Ha, r' and i')
   //----------------------------------------------

	vector<iphas_obj> iphas_colours;
	ifstream iphas_data;
	iphas_data.open(filename.c_str());
	if(!iphas_data) { //output file couldn't be opened
		cerr << "Error: file could not be opened \n";
		exit(1);
	}	
	while (!iphas_data.eof())				// Running down file
	{
		string str1; 	
		getline(iphas_data, str1);
		string temp;
		stringstream sstest(str1);
		sstest>>temp;
		if (temp!="#")
		{		
			double buffer;
			stringstream ss1(str1);
		
			vector<double> infromfile;
		
			while (ss1>>buffer){			// Includes implicit conversion from string to double
				infromfile.push_back(buffer);
			}// correct length		-r_class is stellar or prob stellar-	--i_class is stellar or prob stellar-----	-ha class is stellar or prob stellar---		--------r,i,Ha photometry non-zero----------------		-----r,i,Ha photometric errors non-zero-------------------		---------small RA & DEC offsets between r and i, and r and Ha-------------------
			if (infromfile.size()==30 && (infromfile[6]==-1 || infromfile[6]==-2) && (infromfile[11]==-1 || infromfile[11]==-2) && (infromfile[16]==-1 || infromfile[16]==-2) && (infromfile[4]!=0) && (infromfile[9]!=0) && (infromfile[14]!=0) && (infromfile[5]!=0) && (infromfile[10]!=0) && (infromfile[15]!=0) && infromfile[18]<=1.0 && infromfile[19]<=1.0 && infromfile[20]<=1.0 && infromfile[21]<=1.0) 	//selecting only stellar or probably stellar objects and those with small RA & DEC offsets
			{
            //iphas_obj::iphas_obj(double r_input, double i_input, double ha_input, double d_r_input,double d_i_input, double d_ha_input, double l_input, double b_input)
   		        	if(infromfile[4] > r_max1) { r_max1 = infromfile[4];}
            			if(infromfile[9] > i_max1) { i_max1 = infromfile[9];}
            			if(infromfile[14] > ha_max1) { ha_max1 = infromfile[14];}
			        iphas_obj next_obj(infromfile[4], infromfile[9], infromfile[14], sqrt(pow(infromfile[5],2)+pow(0.0016165105,2)), sqrt(pow(infromfile[10],2)+pow(0.0016165105,2)), sqrt(pow(infromfile[15],2)+pow(0.0016165105,2)), infromfile[0], infromfile[1], infromfile[22], infromfile[23]*1.0857, infromfile[28], infromfile[27], infromfile[29]);	// making the new iphas_obj

				iphas_colours.push_back(next_obj);																											// pushing it into the vetor
			}

			if (infromfile.size()==22 && (infromfile[6]==-1 || infromfile[6]==-2) && (infromfile[11]==-1 || infromfile[11]==-2) && (infromfile[16]==-1 || infromfile[16]==-2) && (infromfile[4]!=0) && (infromfile[9]!=0) && (infromfile[14]!=0) && (infromfile[5]!=0) && (infromfile[10]!=0) && (infromfile[15]!=0) && infromfile[18]<=1.0 && infromfile[19]<=1.0 && infromfile[20]<=1.0 && infromfile[21]<=1.0) 	//selecting only stellar or probably stellar objects and those with small RA & DEC offsets
			{
            //iphas_obj::iphas_obj(double r_input, double i_input, double ha_input, double d_r_input,double d_i_input, double d_ha_input, double l_input, double b_input)
   		        	if(infromfile[4] > r_max1) { r_max1 = infromfile[4];}
            			if(infromfile[9] > i_max1) { i_max1 = infromfile[9];}
            			if(infromfile[14] > ha_max1) { ha_max1 = infromfile[14];}
			        iphas_obj next_obj(infromfile[4], infromfile[9], infromfile[14], sqrt(pow(infromfile[5],2)+pow(0.0016165105,2)), sqrt(pow(infromfile[10],2)+pow(0.0016165105,2)), sqrt(pow(infromfile[15],2)+pow(0.0016165105,2)), infromfile[0], infromfile[1]);	// making the new iphas_obj

				iphas_colours.push_back(next_obj);																											// pushing it into the vetor
			}

			else if (infromfile.size()==22)		// if sources are saturated alter bright limits
			{
				if (infromfile[6]==-9 && infromfile[4]>r_min1){r_min1=infromfile[4];}		//r'
				if (infromfile[11]==-9 && infromfile[9]>i_min1){i_min1=infromfile[9];}		//i'
				if (infromfile[16]==-9 && infromfile[14]>ha_min1){ha_min1=infromfile[14];} 	//Ha
			}
		}	
	}
	iphas_data.close();
	return  iphas_colours;
}

void output_write(string filename, vector<bin_obj2> A_mean, vector<iphas_obj> colours)		// function for outputing results
{
	filename.erase(filename.size()-4);
	ofstream A_out;
	string dummy_string;
	dummy_string=filename+".td4";
	A_out.open(dummy_string.c_str(), ios::trunc);
	//A_out << "#\tdist\tA\tsigma_A\n";
	int last_good=-2;
	for (int x=0; x<A_mean.size(); x++)
	{
	/*	if (A_mean[x].size>=10. && A_mean[x].error_measure<=0.15)
		{
			if (last_good>=0 && last_good!=x-1)
			{
				for (int y=last_good; y<x; y++)
				{
					A_out << y*100 + 50 << "\t" << A_mean[y].mean_A << "\t" << A_mean[y].sigma <<"\t"<<A_mean[y].d_mean<<"\t"<<A_mean[y].d_sigma<<"\t"<<A_mean[y].size<<"\t"<<A_mean[y].error_measure<<"\t"<<A_mean[y].sum<<"\t"<<A_mean[y].diff<<"\t"<<A_mean[y].d_diff<<"\n";				
				}
			}
			else
			{
				A_out << x*100 + 50 << "\t" << A_mean[x].mean_A << "\t" << A_mean[x].sigma <<"\t"<<A_mean[x].d_mean<<"\t"<<A_mean[x].d_sigma<<"\t"<<A_mean[x].size<<"\t"<<A_mean[x].error_measure<<"\t"<<A_mean[x].sum<<"\t"<<A_mean[x].diff<<"\t"<<A_mean[x].d_diff<<"\n";
			}
			last_good=x;
		}*/
	A_out << x*100 + 50 << "\t" << A_mean[x].mean_A << "\t" << A_mean[x].sigma <<"\t"<<A_mean[x].d_mean<<"\t"<<A_mean[x].d_sigma<<"\t"<<A_mean[x].size<<"\t"<<A_mean[x].error_measure<<"\t"<<A_mean[x].sum<<"\t"<<A_mean[x].diff<<"\t"<<A_mean[x].d_diff<<"\n";
	}
	A_out.close();

	vector < vector <double> > previous_rel (150, vector <double> (4));

	ofstream output;
	dummy_string=filename+"-090.dat";
	output.open(dummy_string.c_str(), ios::trunc);
	output << "#\tr\ti\tha\tr_i0\tdist\tA\tdistbin\td_A\td_r_i0\td_dist\td_r\td_i\td_ha\tmag_weight\tprob\tA_prob\tMi\tlogAge\tfeh\td_Mi\td_lagAge\td_feh\tlogT\tlogg\trx\tix\thax\n" ;
	for (int y=0; y<colours.size(); y++)
	{
		
		output << colours[y].r << "\t" << colours[y].i << "\t" << colours[y].ha << "\t"<< colours[y].r_i0 << "\t" << colours[y].dist << "\t" << colours[y].A << "\t" << colours[y].distbin << "\t" << colours[y].d_A << "\t" << colours[y].d_r_i0 << "\t" << colours[y].d_dist << "\t"
		 << colours[y].d_r << "\t" << colours[y].d_i << "\t" << colours[y].d_ha << "\t" << colours[y].last_iso.Mi  << "\t" << colours[y].mean_prob  << "\t" << colours[y].mean_A_prob  << "\t" << colours[y].Mi  << "\t" << colours[y].logAge  << "\t" << colours[y].feh  << "\t" << colours[y].d_Mi  << "\t" << colours[y].d_logAge  << "\t" << colours[y].d_feh << "\t"
		 << colours[y].logT << "\t" << colours[y].logg << "\t" << colours[y].rx << "\t" << colours[y].ix << "\t" << colours[y].hax << "\n";
	}
	output.close();//*/

/*	ofstream short_out;
	dummy_string=filename+"-short.dat";
	short_out.open(dummy_string.c_str(), ios::trunc);
	short_out << "#\tl\tb\tdist\tA\tr_i0\tlumclass\n";
	for (int z=0; z<colours.size(); z++)
	{
		short_out << colours[z].l << "\t" << colours[z].b << "\t" << colours[z].dist << "\t" << colours[z].A << "\t" << colours[z].r_i0 << "\t" << colours[z].lumclass << "\n";
	}S
	short_out.close();*/
} 

// Faster CDF of normal dist
// Based on Abromowitz and Stegun Handbook of Mathematical Functions
// And http://inside.mines.edu/~ckarlsso/codes/codefiles/cdf.cpp
//

double cdf_normal_fast(double x, double mu, double sigma)
{
  const double b1 =  0.319381530;
  const double b2 = -0.356563782;
  const double b3 =  1.781477937;
  const double b4 = -1.821255978;
  const double b5 =  1.330274429;
  const double p  =  0.2316419;
  double Zx;

  x=(x-mu)/sigma;
  Zx=exp(-x*x/2)/sqrt(2);

  if(x >= 0.0) {
      double t = 1.0 / (1.0 + p*x);
      return (1.0 - Zx*t* 
      (t*(t*(t*(t*b5 + b4) + b3) + b2) + b1));
  } 
  else { 
      double t = 1.0 / ( 1.0 - p * x );
      return ( Zx*t* 
      (t*(t*(t*(t*b5 + b4) + b3) + b2) + b1));
  }
}

// Normal CDF for x << mu
// Based on http://www.rfcafe.com/references/mathematical/erf-erfc.htm

double cdf_normal_smallx(double x, double mu, double sigma)
{
	double q=(mu-x)/pow(sigma,1);
	return -log(q*sqrt(2*PI)) + log(1-pow(q,-2)) -pow(q,2)/2;
}


// Faster Inverse CDF of a normal distn
// Based on Abromowitz and Stegun Handbook of Mathematical Functions
// And 

double inv_cdf_normal_fast(double p, double mu, double sigma)
{
    const double c[] = {2.515517, 0.802853, 0.010328};
    const double d[] = {1.432788, 0.189269, 0.001308};
    double t;

    if (p>0.5){t=sqrt(-2.*log(1-p));}
    else {t=sqrt(-2.*log(p));}

    return mu + sigma*(t - ((c[2]*t + c[1])*t + c[0]) / 
           (((d[2]*t + d[1])*t + d[0])*t + 1.0));
}


// cumulative distribution function of a normal distribution
double cdf_normal (double x, double mu, double sigma)		
{
	return ((1+erf((x-mu)/(sqrt(2)*sigma)))/2);
}

// Inverse of the cdf of a normal distribution
//double inv_cdf_normal(double p, double mu, double, sigma)
//{
//	return mu + sigma*sqrt(2)* inv_erf(2*p-1);
//}



// gets stdout following the use of a command - used to read in Schlegel Galactic reddening - code 'borrowed' from somewhere on the net
string getStdoutFromCommand(string cmd)				
{
	// setup
	string str;
	FILE *stream;
	char buffer[100];

	// do it
	stream = popen(cmd.c_str(), "r");
	while ( fgets(buffer, 100, stream) != NULL )
		str.append(buffer);
	pclose(stream);

	// exit
	return str;
}


vector<bin_obj2> backup_A_mean_find(double l_gal, double b_gal)
{
	double Sch_max, density_dust, A_6250;
	vector<bin_obj2> backup_A_mean (150);

	// retrieve Schlegel et al limit
	string Sch_string="./CodeC/lambert_getval CodeC/SFD_dust_4096_ngp.fits CodeC/SFD_dust_4096_sgp.fits 1 "; 	
	Sch_string.append(stringify(l_gal));
	Sch_string.append(" ");
	Sch_string.append(stringify(b_gal));
	Sch_max=atof(getStdoutFromCommand(Sch_string).c_str())*2.944;		// 2.944 to convert E(B-V) given by Schlegel to A_6250
	Sch_max=4.0;
//   cout << "Sch_max = " << Sch_max << endl;

	// integrate dust density to ~infinity, used to normalise the dust distribution so that at infinity it gives the Schlegel value
	double dust_inf=0;
	for (double d=0; d<50000; d+=10)
	{
		dust_inf+=exp(-sqrt(pow(8080.,2)+pow(d*cos(b_gal*PI/180.),2)-2.*8080.*d*cos(b_gal*PI/180.)*cos(l_gal*PI/180.))/2500 - fabs(d*sin(b_gal*PI/180.)+17)/125)*10;	// Dust scale height and lengh from Marshall et al 2006
	}

	double const_term=Sch_max/dust_inf;

	A_6250=0;
	for (double d=0.0; d<=15001.0; d+=10.0)
	{
		density_dust=exp(-sqrt(pow(8080.,2)+pow(d*cos(b_gal*PI/180.),2)-2.*8080.*d*cos(b_gal*PI/180.)*cos(l_gal*PI/180.))/2500 - fabs(d*sin(b_gal*PI/180.)+17)/125);
		A_6250+=const_term*10*density_dust;		// max/total_int * delta_d * rho(d)

		if (d/100!=int(d/100) && d/50==int(d/50))				
		{
         		//cout << "d/100=" << d/100 << " A=" << A_6250 << " " << backup_A_mean.size() << endl; 
										// also make backup_A_mean at this point
			backup_A_mean[int((d-50)/100)].mean_A=A_6250;
			backup_A_mean[int((d-50)/100)].sigma=0.1*A_6250;
			backup_A_mean[int((d-50)/100)].d_mean=0.1;
		}
	}
	return backup_A_mean;   
}

double int_lookup(double A_max, double A_mean, double sd)
{
	if (A_mean>=9.95){A_mean=9.94;}
	if (A_mean<0.05){A_mean=0.05;}
	if (A_max>=9.95){A_max=9.94;}
	if (A_max<-1.95){A_max=-1.95;}
	if (sd>=1.95){sd=1.94;}
	if (sd<0.05){sd=0.05;}
	return lookup_table[int(floor((A_max+2)*10.-0.5))][int(floor(A_mean*10.-0.5))][int(floor(sd*10.-0.5))];
}

double integral_func (double *A_test, size_t dim, void *params)
{
	params_struct *p;
	p=(params_struct *)params;
	double xi=sqrt(log(1+pow(p->sigma/(p->A_mean),2)));
	double mu= log(p->A_mean)-xi/2;
	//double cdf= gsl_cdf_lognormal_P(p->A_max, (mu), xi);
	return gsl_ran_lognormal_pdf(*A_test, (mu),  xi) / (1 + exp(4*(*A_test-p->A_max)))  ;
}


vector <vector <vector <double> > > lookup_creator(void)
{
	vector <vector <vector <double> > > dummy_table(120, vector <vector <double> > (100, vector <double> (20, 0)));

	struct params_struct {double A_max; double A_mean; double sigma;};

	double res=1, err;

	double low[1]={0.};
	double hi[1]={10.};

	params_struct int_params;

	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (1);

	for (int A_max_it=0; A_max_it<120; A_max_it++)
	{
		int_params.A_max=A_max_it*0.1+0.1-2.;
	//	cout << A_max_it << endl;

		for (int A_mean_it=0; A_mean_it<100; A_mean_it++)
		{
			int_params.A_mean=A_mean_it*0.1+0.1;

			for (int sd_it=0; sd_it<20; sd_it++)
			{
				int_params.sigma=sd_it*0.1+0.1;
				gsl_monte_function F={&integral_func, 1, &int_params};
				gsl_monte_vegas_init(s);

				gsl_monte_vegas_integrate (&F, low, hi, 1, 100, rng_handle, s, &res, &err);
			
				dummy_table[A_max_it][A_mean_it][sd_it]=res;
			}
		}
	}

	gsl_monte_vegas_free(s);

	return dummy_table;
}

