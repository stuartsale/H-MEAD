#include <cstdlib> 
#include <ctime>
//#include "helper.h"
//#include "bin_obj.h"
//#include "iso_obj.h"
//#include "iphas_obj.h"
#include "sl_obj.h"
//#include "omp.h"
//#include <mpi.h>
using namespace std;

#ifndef PI
#define PI 3.14159265358979323846
#endif
#ifndef log2PIover2
#define log2PIover2 0.918938533
#endif

//Reads in IPHAS datas




double r_min, i_min, ha_min, r_max, i_max, ha_max;

vector <vector <vector <double> > > lookup_table;

vector <sl_obj> slsl;

	// Wider disc params

		float previous_s_R, previous_s_z, previous_A_0;
		float test_s_R, test_s_z, test_A_0;
		vector <float> s_R_chain, s_z_chain, A_0_chain;
		float s_R_mean, s_z_mean, A_0_mean;
		float s_R_sd, s_z_sd, A_0_sd;
		float previous_rho_prob, test_rho_prob;

void hyperprior_update_all(vector <LF> &LFs);
void mean_intervals(void);
void acl_calc(void);
void neighbour_find(vector<sl_obj> &sl_list);
	float last_asis_prob, test_asis_prob, last_norm_prob_sum, test_norm_prob_sum;

	float last_log_sum, test_log_sum;
	float last_theta_prior_prob, test_theta_prior_prob;

int gal_update=0, gal_update2=0;

gsl_rng* rng_handle;

string config_dir;

//--------------------------------
// MAIN 
//--------------------------------

//Mead::Mead(string ifname){

int main(int argc, char* argv[]) 
{	
	// Set up

	char hostname[1024];
	gethostname(hostname, 1023);
	string hostname_s=hostname;
	if (hostname_s=="orion")
	{
		cout << "hostname is " << hostname_s << endl;
		config_dir="/home/stuart/work/work-oxford/distance_red/config/";
	}
	else
	{
		cout << "hostname is " << hostname_s << endl;	
		config_dir="/usersVol1/sale/distance_red/config/";
	}


	// Initialise MPI

//	int numprocs, rank, namelen, rc;
//	char processor_name[MPI_MAX_PROCESSOR_NAME];
//	MPI_Status Stat;

//	rc=MPI_Init(&argc, &argv);
//	if (rc!=MPI_SUCCESS)
//	{
//		printf ("Error starting MPI program. Terminating.\n");
//		MPI_Abort(MPI_COMM_WORLD, rc);
//	}

//	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	MPI_Get_processor_name(processor_name, &namelen);

//	printf("Process %d on %s out of %d\n", rank, processor_name, numprocs);

	//initialize up random number generator
    	//gsl_rng_env_setup ();
	 rng_handle = gsl_rng_alloc (gsl_rng_taus2);

	//seed the random no generator
	gsl_rng_set(rng_handle, time(0));


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
	r_max=21.;
	i_max=20.;
	ha_max=20.;

      //Reading in isochrones data into a vector

	//vector<iso_obj> isochrones=iso_read("padova-iso_reg.dat");
	//vector<iso_obj> isochrones=iso_read_Tg("padova-iso_tefflogg3.dat");
	vector<iso_obj> isochrones=iso_read_Tg_2MASS(config_dir+"padova-iso_tefflogg-JHK.dat");

	vector<iso_obj> guess_set;
	guess_set.push_back(iso_get_Tg(0.,3.574 ,5.00 , isochrones));	//M1
	guess_set.push_back(iso_get_Tg(0.,3.591 ,4.95 , isochrones));	//M0
	guess_set.push_back(iso_get_Tg(0.,3.602 ,4.90 , isochrones));	//K7
	guess_set.push_back(iso_get_Tg(0.,3.643 ,4.65 , isochrones));	//K5
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
	guess_set.push_back(iso_get_Tg(0.,3.982 ,4.07 , isochrones));	//A0 
	guess_set.push_back(iso_get_Tg(0.,4.061 ,4.07 , isochrones));	//B8
	guess_set.push_back(iso_get_Tg(0.,4.188 ,4.10 , isochrones));	//B5 
	guess_set.push_back(iso_get_Tg(0.,4.362 ,4.06 , isochrones));	//B2 	
//	guess_set.push_back(iso_get_Tg(0.,4.498 ,4.00 , isochrones));	//B0 

	guess_set.push_back(iso_get_Tg(0.,3.580 ,1.41 , isochrones));	//M1 III
	guess_set.push_back(iso_get_Tg(0.,3.591 ,1.63 , isochrones));	//M0 III	
	guess_set.push_back(iso_get_Tg(0.,3.602 ,1.93 , isochrones));	//K5 III	
	guess_set.push_back(iso_get_Tg(0.,3.613 ,2.16 , isochrones));	//K4 III	
	guess_set.push_back(iso_get_Tg(0.,3.628 ,2.36 , isochrones));	//K3 III	
	guess_set.push_back(iso_get_Tg(0.,3.648 ,2.63 , isochrones));	//K2 III	
	guess_set.push_back(iso_get_Tg(0.,3.663 ,2.78 , isochrones));	//K1 III	
	guess_set.push_back(iso_get_Tg(0.,3.681 ,2.89 , isochrones));	//K0 III		
	guess_set.push_back(iso_get_Tg(0.,3.695 ,2.95 , isochrones));	//G8 III		
	guess_set.push_back(iso_get_Tg(0.,3.712 ,3.07 , isochrones));	//G5 III		
	guess_set.push_back(iso_get_Tg(0.,3.740 ,3.20 , isochrones));	//G2 III		


	vector <LF> lfs;

//	LF lfzero("config/iphas_LFs/lfp0000_r.dat");
//	lfs.push_back(lfzero);
	LF lfzero(config_dir+"2MASS_LFs/lfp0000_J.dat");
	lfs.push_back(lfzero);

//	while (1.0648*guess_set[guess_set.size()-1].Mi<2.060){guess_set.push_back(iso_get(0., 1.0648*guess_set[guess_set.size()-1].Mi, 8.5, isochrones));}
//	guess_set.push_back(iso_get(0., 2.060, 8.5, isochrones));

// Set up initial disc model 

	previous_s_R=2500.;
	previous_s_z=125.;
	previous_A_0=0.19;

	cout << previous_s_R << " " <<previous_s_z << " " <<previous_A_0 << endl;

// Read in config file

	vector <vector <string> > config_file;
	config_file=config_read(argv[1]);

	slsl.resize(config_file.size());

	for(int it_conf=0; it_conf<config_file.size(); it_conf++)
	{

   // Read in data

		sl_obj sl1( config_file[it_conf][0],atof(config_file[it_conf][1].c_str()), atof(config_file[it_conf][2].c_str()), config_file[it_conf][3], previous_s_R, previous_s_z, it_conf);
		//slsl.push_back(sl1);
		slsl[it_conf]=sl1;
	}
	neighbour_find(slsl);
	for(int it_conf=0; it_conf<config_file.size(); it_conf++)
	{
	//	if (it_conf!=0){slsl[it_conf].neighbour_set(&slsl[it_conf-1]);}
		slsl[it_conf].initial_guess(isochrones, guess_set, lfs, previous_s_R, previous_s_z, previous_A_0);
	}



	clock_t start;
	start=time(NULL);
		ofstream trace1;
		trace1.open((slsl[0].rootname+".trc").c_str(), ios::trunc);

	while (slsl[0].it_num<240000)
	{
		#pragma omp parallel for //num_threads(2)
		for (int it_conf=0; it_conf<config_file.size(); it_conf++)
		{
			slsl[it_conf].update(isochrones, lfs);
		}
		hyperprior_update_all(lfs);

		if (slsl[0].it_num/10.==floor(slsl[0].it_num/10.))
		{
			trace1 << slsl[0].it_num << " " << previous_s_R << " " << previous_s_z << " " << previous_A_0 << " " << previous_rho_prob << " " << last_asis_prob << " " << slsl[0].previous_norm_prob << " " << gal_update/slsl[0].it_num << " " << gal_update2/slsl[0].it_num << " " << slsl[0].running_A_mean[149].last_mean_A << " " << slsl[0].global_previous_prob << " " << slsl[0].star_prob << " " << last_theta_prior_prob <<endl ;// " " << slsl[0].global_previous_prob << endl ;
		}
	}		

		trace1.close();	
	cout << "total time: " << (time(NULL)-start) <<"s\n";
   	
// Write results to file

	mean_intervals();
	acl_calc();
	for(int it_conf=0; it_conf<config_file.size(); it_conf++)
	{
		slsl[it_conf].output_write(s_R_mean, s_z_mean);
	}


//	MPI_Finalize();
	return 0;
}



//-------------------------------------------------------------------
//-------------------------------------------------------------------
//       HELPER FUNCTIONS
//-------------------------------------------------------------------
//-------------------------------------------------------------------



void hyperprior_update_all(vector <LF> &LFs)
{
	test_rho_prob=0.;
	previous_rho_prob=0.;	
	last_theta_prior_prob=0;//slsl.size()*150*(-log(previous_s_R))/3.;

	for (int it=0; it<slsl.size(); it++){previous_rho_prob+=slsl[it].get_rho_last_prob_higher();}// cout << it << " " << slsl[it].get_rho_last_prob_higher() << endl;}

	test_s_R=previous_s_R+gsl_ran_gaussian_ziggurat(rng_handle,4.);
	test_s_z=previous_s_z;//+gsl_ran_gaussian_ziggurat(rng_handle,1.);
	test_A_0=previous_A_0+gsl_ran_gaussian_ziggurat(rng_handle,0.001)-(test_s_R-previous_s_R)*0.000212;

	test_theta_prior_prob=0;//slsl.size()*150*(-log(test_s_R))/3.;

	//cout << log(previous_s_R) << " " << log(test_s_R) << " " << last_theta_prior_prob << " " << test_theta_prior_prob << " " << slsl.size() << endl;

	//cout << test_s_R << " "<< test_s_z << " " << test_A_0 << endl;
	for (int it=slsl.size()-1; it>-1; it--)
	{
		slsl[it].make_new_test_m_vec(test_s_R, test_s_z, test_A_0);
		test_rho_prob+=slsl[it].get_rho_test_prob_higher();
	//	cout << it << " " << slsl[it].get_rho_test_prob_higher() << " " << slsl[it].get_rho_last_prob_higher() << endl;
	}
	
	if (test_rho_prob+test_theta_prior_prob > previous_rho_prob+last_theta_prior_prob)
	{
		previous_s_R=test_s_R;
		previous_s_z=test_s_z;
		previous_A_0=test_A_0;
		for (int it=0; it<slsl.size(); it++){slsl[it].last_m_vec=slsl[it].test_m_vec;}
		gal_update++;
	//	cout << "pass1 " << test_rho_prob << " " << previous_rho_prob << " " << slsl[0].test_m_vec[50]<< " " << slsl[0].last_m_vec[50] << endl;
		last_theta_prior_prob =test_theta_prior_prob ;
	}
	else if (exp(test_rho_prob+test_theta_prior_prob - (previous_rho_prob+last_theta_prior_prob))>gsl_ran_flat(rng_handle, 0., 1.) )
	{
		previous_s_R=test_s_R;
		previous_s_z=test_s_z;
		previous_A_0=test_A_0;
		for (int it=0; it<slsl.size(); it++){slsl[it].last_m_vec=slsl[it].test_m_vec;}
		gal_update++;
		last_theta_prior_prob =test_theta_prior_prob ;
	//	cout << "pass2 " << test_rho_prob << " " << previous_rho_prob << " " << slsl[0].test_m_vec[50]<< " " << slsl[0].last_m_vec[50] << endl;
	}
	//else {cout << "fail " << test_rho_prob << " " << previous_rho_prob << " " << slsl[0].test_m_vec[50]<< " " << slsl[0].last_m_vec[50] << endl;}


// -------------------------------------------------------------------------------------------------------------------------------	
// Now find z_dash|rho, Theta

//	last_asis_prob=0; test_asis_prob=0; last_norm_prob_sum=0; test_norm_prob_sum=0;
//	last_log_sum=0; test_log_sum=0;

//	last_theta_prior_prob=slsl.size()*150*(-log(previous_s_R));

//	for (int it=0; it<slsl.size(); it++)
//	{
//	//	last_asis_prob+=slsl[it].get_A_mean_last_prob();
//		for (int it2=0; it2<slsl[it].star_cat.size(); it2++)
//		{
//			last_asis_prob+=slsl[it].star_cat[it2].likelihood_eval(slsl[it].star_cat[it2].last_iso, slsl[it].star_cat[it2].last_A, slsl[it].star_cat[it2].last_dist_mod);
//			slsl[it].star_cat[it2].get_last_z();		
//		}
//		
//		last_norm_prob_sum+=slsl[it].previous_norm_prob;
//		slsl[it].last_z_dash=slsl[it].get_last_z_dash();
//	}

//// Update Theta|z_dash

//	if(slsl[0].it_num>0000)
//	{
//		test_s_R=previous_s_R+gsl_ran_gaussian_ziggurat(rng_handle,50.);
//		test_s_z=previous_s_z;//+gsl_ran_gaussian_ziggurat(rng_handle,1.);
//		test_A_0=previous_A_0;//+gsl_ran_gaussian_ziggurat(rng_handle,0.01);
//	}
//	else
//	{
//		test_s_R=previous_s_R;//+gsl_ran_gaussian_ziggurat(rng_handle,10.);
//		test_s_z=previous_s_z;//+gsl_ran_gaussian_ziggurat(rng_handle,1.);
//		test_A_0=previous_A_0;//+gsl_ran_gaussian_ziggurat(rng_handle,0.01)-(test_s_R-previous_s_R)*0.000212;
//	}

//	test_theta_prior_prob=slsl.size()*150*(-log(test_s_R));

//	for (int it=slsl.size()-1; it>-1; it--)	{slsl[it].make_new_test_m_vec(test_s_R, test_s_z, test_A_0);}


//	for (int it=slsl.size()-1; it>-1; it--)	
//	{
//		slsl[it].mvn_gen_internal_rel_from_z_dash();
//	//	test_asis_prob+=slsl[it].get_A_mean_test_prob();

//		for (int it2=0; it2<slsl[it].star_cat.size(); it2++)
//		{
//			slsl[it].star_cat[it2].set_test_A_from_z();		
//			test_asis_prob+=slsl[it].star_cat[it2].likelihood_eval(slsl[it].star_cat[it2].last_iso, slsl[it].star_cat[it2].last_A, slsl[it].star_cat[it2].last_dist_mod);
//		}
//	
//	// Normalisation term
//		slsl[it].current_norm_prob=0;
//		for (int it_LF=0; it_LF<LFs.size(); it_LF++)
//		{
//			slsl[it].current_norm_prob+=-LFs[it_LF].LF_prob_test(slsl[it].running_A_mean, slsl[it].J_min, slsl[it].J_max)*(slsl[it].star_cat.size()+1);
//		}
//		test_norm_prob_sum+=slsl[it].current_norm_prob;
//	}


//	if (test_asis_prob+test_norm_prob_sum+test_theta_prior_prob> last_asis_prob+last_norm_prob_sum+last_theta_prior_prob)
//	{
////	cout << "1 " << slsl[0].running_A_mean[0].A_chain.size() << " " << (test_asis_prob+test_norm_prob_sum - last_asis_prob-last_norm_prob_sum) << " " << last_asis_prob << " " << test_asis_prob  << " " << slsl[0].running_A_mean[50].last_sd_A << " " << slsl[0].running_A_mean[50].test_sd_A <<  " " << slsl[0].last_m_vec[50] << " " << slsl[0].test_m_vec[50] << " " << last_norm_prob_sum << " " << test_norm_prob_sum << " " << previous_s_R << " " << test_s_R  << endl;
//		previous_s_R=test_s_R;
//		previous_s_z=test_s_z;
//		previous_A_0=test_A_0;
//		for (int it=0; it<slsl.size(); it++)
//		{
//			slsl[it].last_m_vec=slsl[it].test_m_vec;
//			slsl[it].global_previous_prob==slsl[it].get_A_mean_test_prob();
//			slsl[it].previous_norm_prob=slsl[it].current_norm_prob;
//			for (int it2=0; it2<slsl[it].running_A_mean.size(); it2++){slsl[it].running_A_mean[it2].accept();}
//		}
//		gal_update2++;
//		last_theta_prior_prob =test_theta_prior_prob ;
//	}
//	else if (exp(test_asis_prob+test_norm_prob_sum+test_theta_prior_prob - last_asis_prob-last_norm_prob_sum-last_theta_prior_prob)>gsl_ran_flat(rng_handle, 0., 1.) )
//	{
////	cout << "1 " << slsl[0].running_A_mean[0].A_chain.size() << " " << (test_asis_prob+test_norm_prob_sum - last_asis_prob-last_norm_prob_sum) << " " << last_asis_prob << " " << test_asis_prob  << " " << slsl[0].running_A_mean[50].last_sd_A << " " << slsl[0].running_A_mean[50].test_sd_A <<  " " << slsl[0].last_m_vec[50] << " " << slsl[0].test_m_vec[50] << " " << last_norm_prob_sum << " " << test_norm_prob_sum << " " << previous_s_R << " " << test_s_R  << endl;
//		previous_s_R=test_s_R;
//		previous_s_z=test_s_z;
//		previous_A_0=test_A_0;
//		for (int it=0; it<slsl.size(); it++)
//		{
//			slsl[it].last_m_vec=slsl[it].test_m_vec;
//			slsl[it].global_previous_prob==slsl[it].get_A_mean_test_prob();
//			slsl[it].previous_norm_prob=slsl[it].current_norm_prob;
//			for (int it2=0; it2<slsl[it].running_A_mean.size(); it2++){slsl[it].running_A_mean[it2].accept();}
//		}
//		gal_update2++;

//	}
//	else
//	{	
////	cout << "0 " << slsl[0].running_A_mean[0].A_chain.size() << " " << (test_asis_prob+test_norm_prob_sum - last_asis_prob-last_norm_prob_sum) << " " << last_asis_prob << " " << test_asis_prob  << " " << slsl[0].running_A_mean[50].last_sd_A << " " << slsl[0].running_A_mean[50].test_sd_A <<  " " << slsl[0].last_m_vec[50] << " " << slsl[0].test_m_vec[50] << " " << last_norm_prob_sum << " " << test_norm_prob_sum << " " << previous_s_R << " " << test_s_R  << endl;
//		test_s_R=previous_s_R;
//		test_s_z=previous_s_z;
//		test_A_0=previous_A_0;

//		for (int it=0; it<slsl.size(); it++)
//		{
//			slsl[it].test_m_vec=slsl[it].last_m_vec;
//			slsl[it].global_previous_prob==slsl[it].get_A_mean_last_prob();
//			slsl[it].current_norm_prob=slsl[it].previous_norm_prob;
//			for (int it2=0; it2<slsl[it].running_A_mean.size(); it2++){slsl[it].running_A_mean[it2].reject();}
//		}
//	}

//	if (floor(slsl[0].it_num/100.)==slsl[0].it_num/100)
//	{
//		s_R_chain.push_back(previous_s_R);
//		s_z_chain.push_back(previous_s_z);
//		A_0_chain.push_back(previous_A_0);
//	}
	

}

void mean_intervals(void)
{
	for (int it=0; it<slsl.size(); it++){slsl[it].mean_intervals();}

	float s_R_sum=0., s_z_sum=0, A_0_sum=0.;
	float s_R_sum2=0., s_z_sum2=0, A_0_sum2=0.;
	for (int it=floor(0.5*s_R_chain.size()); it<s_R_chain.size(); it++)
	{
		s_R_sum+=s_R_chain[it];
		s_z_sum+=s_z_chain[it];
		A_0_sum+=A_0_chain[it];
		s_R_sum2+=pow(s_R_chain[it],2);
		s_z_sum2+=pow(s_z_chain[it],2);
		A_0_sum2+=pow(A_0_chain[it],2);
	}
	s_R_mean=s_R_sum/ceil(0.5*s_R_chain.size());
	s_z_mean=s_z_sum/ceil(0.5*s_z_chain.size());
	A_0_mean=A_0_sum/ceil(0.5*A_0_chain.size());

	s_R_sd=sqrt( s_R_sum2/ceil(0.5*s_R_chain.size()) - pow(s_R_mean,2));
	s_z_sd=sqrt( s_z_sum2/ceil(0.5*s_z_chain.size()) - pow(s_z_mean,2));
	A_0_sd=sqrt( A_0_sum2/ceil(0.5*A_0_chain.size()) - pow(A_0_mean,2));

	cout << s_R_mean << " " << s_z_mean << " "<< A_0_mean << " " << s_R_sd << " " << s_z_sd << " " << A_0_sd << " " << endl;
}

void acl_calc(void)
{
	ofstream acl_out;
	string dummy_string;
	dummy_string="disc.acl";
	acl_out.open(dummy_string.c_str(), ios::trunc);
	acl_out << "# lag acf_s_R acf_s_z acf_A_0" <<endl; 


	vector <float> acl;
	//vector <float> new_acl;
	vector <float> new_acl(int(ceil(0.5*s_R_chain.size())), 0.);
	vector <float> new_acl2(int(ceil(0.5*s_z_chain.size())), 0.);
	vector <float> new_acl3(int(ceil(0.5*A_0_chain.size())), 0.);

	for (int it1=floor(0.5*s_R_chain.size()); it1<s_R_chain.size(); it1++)
	{
		for (int lag=0; lag<it1-ceil(0.5*s_R_chain.size()); lag++)
		{
			new_acl[lag]+=(s_R_chain[it1]-s_R_mean)*(s_R_chain[it1-lag]-s_R_mean);
			new_acl2[lag]+=(s_z_chain[it1]-s_z_mean)*(s_z_chain[it1-lag]-s_z_mean);
			new_acl3[lag]+=(A_0_chain[it1]-A_0_mean)*(A_0_chain[it1-lag]-A_0_mean);
		}

	}

	for (int it1=0; it1<new_acl.size(); it1++)
	{
		new_acl[it1]/=(new_acl.size());
		new_acl2[it1]/=(new_acl2.size());
		new_acl3[it1]/=(new_acl3.size());
	}

	for (int it=0; it<new_acl.size(); it++)
	{
		acl_out << it << " " << new_acl[it]<< " " << new_acl2[it]<< " " << new_acl3[it] << endl;
	}
	acl_out.close();

	size_t const half_size = s_R_chain.size() / 2;
	vector <float> unburnt_s_R(s_R_chain.begin()+half_size, s_R_chain.end());
	vector <float> unburnt_s_z(s_z_chain.begin()+half_size, s_z_chain.end());
	vector <float> unburnt_A_0(A_0_chain.begin()+half_size, A_0_chain.end());

	cout << acl_block(unburnt_s_R) << " " << acl_block(unburnt_s_z) << " " << acl_block(unburnt_A_0) << " " << endl; 
}

void neighbour_find(vector<sl_obj>  &sl_list)
{
	for (int i=0; i<sl_list.size()-1; i++)
	{
		for (int j=i+1; j<sl_list.size(); j++)
		{
			if ( ( abs(sl_list[i].l-sl_list[j].l)<0.01 && abs(sl_list[i].b-sl_list[j].b)<0.27 )
				|| ( abs(sl_list[i].l-sl_list[j].l)<0.27 && abs(sl_list[i].b-sl_list[j].b)<0.01 ) )
			{
				sl_list[i].neighbour_set(&sl_list[j]) ;
				sl_list[j].neighbour_set(&sl_list[i]) ;
				sl_list[i].higher_neighbour_set(&sl_list[j]) ;
			}
		}
	}
}


double int_lookup(double A_max, double A_mean, double sd)
{
	if (A_mean>=19.9){A_mean=19.89;}
	if (A_mean<0.1){A_mean=0.11;}
	if (A_max>=19.9){A_max=19.89;}
	if (A_max<-1.9){A_max=-1.89;}
	if (sd>=9.9){sd=9.89;}
	if (sd<0.1){sd=0.11;}
//	cout << A_mean << " " << A_max << " " << sd << " " << int(floor((A_max+2)*5.-0.5))<< " " << int(floor(A_mean*5.-0.5)) << " " << int(floor(sd*5.-0.5)) << endl;
	return lookup_table[int(floor((A_max+2)*5.-0.5))][int(floor(A_mean*5.-0.5))][int(floor(sd*5.-0.5))];
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
	vector <vector <vector <double> > > dummy_table(220, vector <vector <double> > (200, vector <double> (50, 0)));

	struct params_struct {double A_max; double A_mean; double sigma;};

	double res=1, err;

	double low[1]={0.};
	double hi[1]={30.};

	params_struct int_params;

	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (1);


	for (int A_max_it=0; A_max_it<110; A_max_it++)
	{
		int_params.A_max=A_max_it*0.2+0.2-2.;
	//	cout << A_max_it << endl;

		for (int A_mean_it=0; A_mean_it<100; A_mean_it++)
		{
			int_params.A_mean=A_mean_it*0.2+0.2;

			for (int sd_it=0; sd_it<50; sd_it++)
			{
				int_params.sigma=sd_it*0.2+0.2;
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

