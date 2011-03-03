#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdlib> 
#include <sstream>
#include <time.h>
#include <ctime>
#include "helper.h"
#include "bin_obj.h"
#include "iso_obj.h"
#include "iphas_obj.h"
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

Uniform U;
Normal Z;
VariLogNormal VLN;


//--------------------------------
// MAIN 
//--------------------------------

//Mead::Mead(string ifname){

int main(int argc, char* argv[]) 
{	

	vector<iphas_obj> colours;		
	vector<bin_obj2> A_mean;
	vector<bin_obj2> backup_A_mean (150);				
	A_mean.reserve(251); 
//	backup_A_mean.reserve(251); 

// Setting up random number generators

	srand((unsigned)time(0));
	double s=double(rand())/(double(RAND_MAX)+1.0);	//Generating a random seed for the MT generator in range 0-1
	//MT urng(s);
	MultWithCarry urng(s);
	Random::Set(urng);
	srand ( time(NULL) );

	//initialize up random number generator
    	gsl_rng_env_setup ();
	gsl_rng* rng_handle = gsl_rng_alloc (gsl_rng_taus2);

	//seed the random no generator
	gsl_rng_set(rng_handle, time(0));


// set default MIN vals
	r_min=13.5; 
	i_min=12.0;
	ha_min=12.0;				
// set default MAX vals
	r_max=0.;
	i_max=0.;
	ha_max=0.;

      //Reading in isochrones data into a vector

	vector<iso_obj> isochrones=iso_read("padova-iso_reg.dat");

	vector<iso_obj> guess_set;
	guess_set.push_back(iso_get(0., 0.738, 8.5, isochrones));
	while (1.0648*guess_set[guess_set.size()-1].Mi<2.015){guess_set.push_back(iso_get(0., 1.0648*guess_set[guess_set.size()-1].Mi, 8.5, isochrones));}
	guess_set.push_back(iso_get(0., 2.015, 8.5, isochrones));

   // Read in IPHAS data
	string iphas_filename=argv[1];		
	colours=iphas_read(iphas_filename,r_min,i_min,ha_min,r_max,i_max,ha_max);
	cout << "r_min=" << r_min << " i_min=" << i_min << " ha_min=" << ha_min << endl; 
	cout << "r_max=" << r_max << " i_max=" << i_max << " ha_max=" << ha_max << endl; 


   //
   // do this after iphas_read so that limiting mags are accurate.

	backup_A_mean=backup_A_mean_find(atof(argv[2]), atof(argv[3]));

//	iphas_obj test_obj(14.821,14.462,14.581, 0.002,0.002,0.002,180,0);
//	test_obj.dist_redMCMC(isochrones, guess_set 180.,0.);

  
   // governs how aggresively objects near mag limits 
   // are downweighted - see section 4.6 in my thesis
	double x_value=2.0;				

// first run through

	cout << "backup size:" << backup_A_mean.size() << endl;

	cout << "Real: " << real_prob(colours, isochrones, guess_set, atof(argv[2]), atof(argv[3]), backup_A_mean, -0.0272, 0.53) << endl;
	A_mean=dist_redMCMC(colours, isochrones, guess_set, atof(argv[2]), atof(argv[3]), backup_A_mean, -0.0272, 0.53);			// }
   	
// Write results to file
	output_write(iphas_filename, A_mean, colours);

	return 0;
}



//-------------------------------------------------------------------
//-------------------------------------------------------------------
//       HELPER FUNCTIONS
//-------------------------------------------------------------------
//-------------------------------------------------------------------

double real_prob(vector<iphas_obj> &stars, vector<iso_obj> &isochrones, vector<iso_obj> &guess_set, double l, double b, vector <bin_obj2> backup_A_mean, double ri_min, double ri_max)
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

}

vector <bin_obj2> dist_redMCMC(vector<iphas_obj> &stars, vector<iso_obj> &isochrones, vector<iso_obj> &guess_set, double l, double b, vector <bin_obj2> backup_A_mean, double ri_min, double ri_max)
{

	vector <bin_obj2> first_bins(150);

	double sigma_fac=0.02;
// Set up

	int without_change=0;
	int thin=100;

	vector <vector <vector <double> > > global_A_chain;

	double global_previous_prob=0;
	double previous_hyperprior_prob=0, current_hyperprior_prob=0;
	double global_current_prob, global_transition_prob;
	double sigma2_LN, mu_LN;

	vector < vector<double> > proposal_sd (150, vector <double> (2));
	vector < vector <double> > previous_rel (150, vector <double> (4));
	vector < vector <double> > internal_rel (150, vector <double> (2));
	vector < vector <double> > previous_internal_rel (150, vector <double> (2));

// Start from backup_A_mean



	for (int i=0; i<150; i++)
	{
		previous_rel[i][0]=backup_A_mean[i].mean_A;
		previous_rel[i][1]=0.2;//backup_A_mean[i].sigma;//0.2;//
		previous_rel[i][3]=sqrt(log(1+pow(previous_rel[i][1]/previous_rel[i][0],2)));
		previous_rel[i][2]=log(previous_rel[i][0])-pow(previous_rel[i][3],2)/2;
	}

	previous_internal_rel[0][0]=log(previous_rel[0][0]);//previous_rel[0][0];//
	previous_internal_rel[0][1]=previous_rel[0][1];
		
	//previous_hyperprior_prob+=-previous_internal_rel[0][0];//-log(previous_internal_rel[0][0]);//
	//previous_hyperprior_prob+=-1*log(previous_internal_rel[0][1]);//-previous_internal_rel[0][1];//
	for (int i=1; i<150; i++)
	{
		previous_internal_rel[i][0]=log(previous_rel[i][0]-previous_rel[i-1][0]);//previous_rel[i][0]-previous_rel[i-1][0];//
		previous_internal_rel[i][1]=previous_rel[i][1];
		
	//	previous_hyperprior_prob+=-previous_internal_rel[i][0];//-log(previous_internal_rel[i][0]);//
	//	previous_hyperprior_prob+=-1*log(previous_internal_rel[i][1]);//-previous_internal_rel[i][1];//
	}

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
		if (stars[it_stars].last_A>0 && stars[it_stars].last_iso.r0-stars[it_stars].last_iso.i0> ri_min  && stars[it_stars].last_iso.r0-stars[it_stars].last_iso.i0<ri_max){it_stars++;}
		else {stars.erase(stars.begin()+it_stars);}
	}

	for (int star_it=0; star_it<stars.size(); star_it++){global_previous_prob+=stars[star_it].last_prob;}
	
// Loop through this section

	vector <double> proposed_probs(stars.size());

	while ((global_A_chain.size()<150000 /*|| accepted1<2000 || accepted2<2000 || accepted3<2000*/) && without_change<4000)
	{
		global_current_prob=0;
		global_transition_prob=0;
		global_previous_prob=0;
		current_hyperprior_prob=0;

		//if (global_A_chain.size()/1000.==floor(global_A_chain.size()/1000.)){cout << global_A_chain.size() << endl;}
		//cout << global_A_chain.size() << endl;
	//	total++;
// First vary parameters for each star

		for (int it=0; it<stars.size(); it++){if (U.Next()>0.9975){stars[it].star_try1(isochrones, l, b, previous_rel);};}

// Now vary hyper-parameters

		for (int it=0; it<stars.size(); it++){global_previous_prob+=stars[it].last_prob;}

		vector < vector <double> > new_rel(150,vector <double> (4));

		for (int it=0; it<150; it++)
		{
			proposal_sd[it][0]=sigma_fac;
			proposal_sd[it][1]=sigma_fac;
		}
		

		for (int it=0; it<150; it++)
		{	
		//	internal_rel[it][0]=VLN.Next(previous_internal_rel[it][0], proposal_sd[it][0]*previous_internal_rel[it][0]);
			internal_rel[it][0]=previous_internal_rel[it][0]+Z.Next()*proposal_sd[it][0];//*previous_internal_rel[it][0]);
		//	while (internal_rel[it][0]<0){internal_rel[it][0]=previous_internal_rel[it][0]+Z.Next()*proposal_sd[it][0];}
			internal_rel[it][1]=0.2;//VLN.Next(previous_internal_rel[it][1], proposal_sd[it][1]*previous_internal_rel[it][1]); //
			
		//	current_hyperprior_prob+=-internal_rel[it][0];//-log(internal_rel[it][0]);//
		//	current_hyperprior_prob+=-1*log(internal_rel[it][1]);//-internal_rel[it][1];//
		}

		new_rel[0][0]=exp(internal_rel[0][0]);//internal_rel[0][0];//
		new_rel[0][1]=internal_rel[0][1];
		new_rel[0][3]=sqrt(log(1+pow(new_rel[0][1]/new_rel[0][0],2)));
		new_rel[0][2]=log(new_rel[0][0])-pow(new_rel[0][3],2)/2;
		for (int it=1; it<150; it++)
		{
			new_rel[it][0]=exp(internal_rel[it][0])+new_rel[it-1][0];//internal_rel[it][0]+new_rel[it-1][0];//
			new_rel[it][1]=internal_rel[it][1];
			new_rel[it][3]=sqrt(log(1+pow(new_rel[it][1]/new_rel[it][0],2)));
			new_rel[it][2]=log(new_rel[it][0])-pow(new_rel[it][3],2)/2;
		}

// Find probability of this parameter set

		for (int it=0; it<stars.size(); it++)
		{
			proposed_probs[it]=stars[it].prob_eval(stars[it].last_iso, stars[it].last_A, stars[it].last_dist_mod, new_rel);
			global_current_prob+= proposed_probs[it];
		}

// Metropolis-Hastings algorithm step

		global_transition_prob=0;

		for (int it=1; it<150; it++)
		{
		// From new to old
		// mean_A
		//	sigma2_LN=log(1+pow(proposal_sd[it][0],2));
		//	mu_LN=log(internal_rel[it][0])-sigma2_LN/2;
		//	global_transition_prob+=-log(previous_internal_rel[it][0]) - pow(log(previous_internal_rel[it][0])-mu_LN,2)/(2*sigma2_LN);
	//		cout << mu_LN << " " << sigma2_LN << " " << internal_rel[it][0] << " " << it << " " <<-log(previous_internal_rel[it][0]) - pow(log(previous_internal_rel[it][0])-mu_LN,2)/(2*sigma2_LN)<< endl;
		// mean_A - normal sampler
		//	global_transition_prob+=-log(1-cdf_normal_fast(0, internal_rel[it][0], proposal_sd[it][0]));
		// sigma
		//	sigma2_LN=log(1+pow(proposal_sd[it][1],2));
		//	mu_LN=log(new_rel[it][1])-sigma2_LN/2;
		//	global_transition_prob+=-log(previous_rel[it][1]) - pow(log(previous_rel[it][1])-mu_LN,2)/(2*sigma2_LN);

			global_transition_prob+=-pow(internal_rel[it][0]-previous_internal_rel[it][0],2)/(2*pow(proposal_sd[it][0],2));

		// From old to new
		// mean_A
		//	sigma2_LN=log(1+pow(proposal_sd[it][0],2));
		//	mu_LN=log(previous_internal_rel[it][0])-sigma2_LN/2;
		//	global_transition_prob-=-log(internal_rel[it][0]) - pow(log(internal_rel[it][0])-mu_LN,2)/(2*sigma2_LN);
		// mean_A - normal sampler
		//	global_transition_prob-=-log(1-cdf_normal_fast(0, previous_internal_rel[it][0], proposal_sd[it][0]));
		// sigma
		//	sigma2_LN=log(1+pow(proposal_sd[it][1],2));
		//	mu_LN=log(previous_rel[it][1])-sigma2_LN/2;
		//	global_transition_prob-=-log(new_rel[it][1] ) - pow(log(new_rel[it][1])-mu_LN,2)/(2*sigma2_LN);	

			global_transition_prob-=-pow(previous_rel[it][0]-previous_rel[it-1][0]-previous_internal_rel[it][0],2)/(2*pow(proposal_sd[it][0],2));
		}	

// Accept or reject
	

	
		if (global_current_prob+current_hyperprior_prob-global_previous_prob-previous_hyperprior_prob+global_transition_prob>0)		// New parameter set better => Accept
		{
			global_A_chain.push_back(new_rel);
			previous_rel=new_rel;
		//	previous_internal_rel=internal_rel;
			global_previous_prob=global_current_prob;
			previous_hyperprior_prob=current_hyperprior_prob;
			without_change=0;

			for (int stars_it=0; stars_it<stars.size(); stars_it++){stars[stars_it].last_prob=proposed_probs[stars_it];}
		
			cout << global_A_chain.size() << " " << global_current_prob << " " << internal_rel[50][1] << " " << new_rel[149][0] << " " << current_hyperprior_prob << endl;
		}

		else if (exp(global_current_prob+current_hyperprior_prob-global_previous_prob-previous_hyperprior_prob+global_transition_prob)>U.Next())	// New set worse => accept with P=P(new)/P(old)
		{
			global_A_chain.push_back(new_rel);
			previous_rel=new_rel;
		//	previous_internal_rel=internal_rel;
			global_previous_prob=global_current_prob;
			previous_hyperprior_prob=current_hyperprior_prob;
			without_change=0;

			cout << global_A_chain.size() << " " << global_current_prob << " " << internal_rel[50][1] << " " << new_rel[149][0] << " " << current_hyperprior_prob <<  endl;
		}
		else 
		{
			without_change++;
		/*	cout << "fail " << global_current_prob << " " << global_previous_prob << " " << global_transition_prob << " " << stars.size() << endl;*/
			global_A_chain.push_back(previous_rel);
		}
	}

	for (int star_it=0; star_it<stars.size(); star_it++){stars[star_it].mean_intervals();}

	for (int it=0; it<150; it++)
	{
		double A_sum=0., sigma_sum=0.;
		for (int m=floor(0.50*global_A_chain.size()); m<global_A_chain.size(); m++)
		{
			A_sum+=global_A_chain[m][it][0];
			sigma_sum+=log(global_A_chain[m][it][1]);
		}
		first_bins[it].mean_A=A_sum/ceil(0.50*global_A_chain.size());
		first_bins[it].sigma=exp(sigma_sum/ceil(0.50*global_A_chain.size()));

		vector <double> A_diffs, sigma_diffs;
		for (int m=floor(0.50*global_A_chain.size()); m<global_A_chain.size(); m++)
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
	for (int it=0; it<stars.size(); it++)
	{	
//		first_bins[stars[it].distbin].size++;
//		first_bins[stars[it].distbin].sum+=stars[it].A/pow(stars[it].d_A,2);
	}
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
         			iso_obj objnew(fromfile[0], fromfile[1], fromfile[2], fromfile[3], fromfile[4], fromfile[5], fromfile[6], fromfile[7]);
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
	output << "#\tr\ti\tha\tr_i0\tdist\tA\tdistbin\td_A\td_r_i0\td_dist\td_r\td_i\td_ha\tmag_weight\tprob\n" ;
	for (int y=0; y<colours.size(); y++)
	{
		
		output << colours[y].r << "\t" << colours[y].i << "\t" << colours[y].ha << "\t"<< colours[y].r_i0 << "\t" << colours[y].dist << "\t" << colours[y].A << "\t" << colours[y].distbin << "\t" << colours[y].d_A << "\t" << colours[y].d_r_i0 << "\t" << colours[y].d_dist << "\t" << colours[y].d_r << "\t" << colours[y].d_i << "\t" << colours[y].d_ha << "\t" << colours[y].mag_weight  << "\t" << colours[y].last_prob << "\n";
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
	int q=(mu-x)/pow(sigma,2);
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
   cout << "Sch_max = " << Sch_max << endl;

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
			backup_A_mean[int((d-50)/100)].mean_A=d/10000.;//A_6250;
			backup_A_mean[int((d-50)/100)].sigma=0.1*A_6250;
			backup_A_mean[int((d-50)/100)].d_mean=0.1;
		}
	}
	return backup_A_mean;   
}

	

