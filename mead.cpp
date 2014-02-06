#include <cstdlib> 
#include <ctime>
#include <unistd.h>
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

vector<iso_obj> iso_read(const string &filename);
vector<iso_obj> iso_read_Tg(const string &filename);


	


double r_min, i_min, ha_min, r_max, i_max, ha_max;
double J_min, H_min, K_min, J_max, H_max, K_max;

string config_dir;
string hmd_dir;
gsl_rng* rng_handle;

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
	else if (hostname_s=="auriga")
	{
		cout << "hostname is " << hostname_s << endl;
		config_dir="/home/sale/work-oxford/distance_red/config/";
	}
	else
	{
		cout << "hostname is " << hostname_s << endl;	
		config_dir="/usersVol1/sale/distance_red/config/";
	}

	char buff[PATH_MAX];
	getcwd( buff, PATH_MAX );
	std::string cwd( buff );

//	// Initialise MPI

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
	vector<iso_obj> isochrones=iso_read_Tg(config_dir+"padova-iso_tefflogg_0214.txt");

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

	vector <LF> lfs;

	LF lfzero(config_dir+"iphas_LFs/lfp0000_r.dat");
	lfs.push_back(lfzero);


// Read in config file

	vector <vector <string> > config_file;
	string hmd_file_path=argv[1];
	config_file=config_read(hmd_file_path);
	
	unsigned found=hmd_file_path.find_last_of("/");
	if (found==string::npos){hmd_dir=cwd+"/";}
	else{hmd_dir=cwd+"/"+hmd_file_path.substr(0,found+1);}

	cout << "hmd_dir: " << hmd_dir  << endl;

//	vector <sl_obj> slsl(config_file.size());

	for(int it_conf=0; it_conf<config_file.size(); it_conf++)
	{

   // Read in data

		sl_obj sl1( hmd_dir+config_file[it_conf][0], atof(config_file[it_conf][1].c_str()), atof(config_file[it_conf][2].c_str()),
			config_file[it_conf][3], atof(config_file[it_conf][4].c_str()), atof(config_file[it_conf][5].c_str()), atof(config_file[it_conf][6].c_str()) );

		for (int lf_it=0; lf_it<lfs.size(); lf_it++){lfs[lf_it].precompute_Aminmax();}

		cout << it_conf << " " << config_file[it_conf][0] ;
		cout << " r_max: " << sl1.r_max << " i_max: " << sl1.i_max << " ha_max: " << sl1.ha_max << endl;

		sl1.initial_guess(isochrones, guess_set, lfs);

		clock_t start;
		start=time(NULL);

		while (sl1.it_num<150000)
		{
				sl1.update(isochrones, lfs);
		}		
	
		cout << "total time: " << (time(NULL)-start) <<"s\n";
   	
// Write results to file

		sl1.mean_intervals();
		sl1.output_write();
	}


//	MPI_Finalize();
	return 0;
}



//-------------------------------------------------------------------
//-------------------------------------------------------------------
//       HELPER FUNCTIONS
//-------------------------------------------------------------------
//-------------------------------------------------------------------


vector<iso_obj> iso_read(const string &filename)		// Function to read in calibration data
{
	ifstream input1;
	input1.open(filename.c_str());
	if(!input1) { //output file couldn't be opened
		cerr << "Error: Isochrone file could not be opened" << endl;
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
		cerr << "Error: Isochrone_Tg file could not be opened" << endl;
		exit(1);
	}

	vector<iso_obj> totalfile;		// making a vector to store the calibration_obj
	
	while (!input1.eof())				// Running down file - reading it in
	{
//		double jac;
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
			if (fromfile.size()==11)		// check there's something in the line
			{
         			iso_obj objnew(fromfile[0], fromfile[1], fromfile[2], fromfile[3], fromfile[4], fromfile[6], fromfile[7], fromfile[8], fromfile[5]);
				totalfile.push_back(objnew);
			}		
		}
	}
	return totalfile;
	input1.close();
}

//double int_lookup(double A_max, double A_mean, double sd)
//{
//	if (A_mean>=11.95){A_mean=11.94;}
//	if (A_mean<0.05){A_mean=0.05;}
//	if (A_max>=14.95){A_max=14.94;}
//	if (A_max<-1.95){A_max=-1.95;}
//	if (sd>=3.95){sd=3.94;}
//	if (sd<0.05){sd=0.05;}
//	return lookup_table[int(floor((A_max+2)*10.-0.5))][int(floor(A_mean*10.-0.5))][int(floor(sd*10.-0.5))];
//}

//double integral_func (double *A_test, size_t dim, void *params)
//{
//	params_struct *p;
//	p=(params_struct *)params;
//	double xi=sqrt(log(1+pow(p->sigma/(p->A_mean),2)));
//	double mu= log(p->A_mean)-xi/2;
//	//double cdf= gsl_cdf_lognormal_P(p->A_max, (mu), xi);
//	return gsl_ran_lognormal_pdf(*A_test, (mu),  xi) / (1 + exp(4*(*A_test-p->A_max)))  ;
//}


//vector <vector <vector <double> > > lookup_creator(void)
//{
//	vector <vector <vector <double> > > dummy_table(170, vector <vector <double> > (120, vector <double> (40, 0)));

//	struct params_struct {double A_max; double A_mean; double sigma;};

//	double res=1, err;

//	double low[1]={0.};
//	double hi[1]={10.};

//	params_struct int_params;

//	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (1);

//	for (int A_max_it=0; A_max_it<170; A_max_it++)
//	{
//		int_params.A_max=A_max_it*0.1+0.1-2.;
//	//	cout << A_max_it << endl;

//		for (int A_mean_it=0; A_mean_it<120; A_mean_it++)
//		{
//			int_params.A_mean=A_mean_it*0.1+0.1;

//			for (int sd_it=0; sd_it<40; sd_it++)
//			{
//				int_params.sigma=sd_it*0.1+0.1;
//				gsl_monte_function F={&integral_func, 1, &int_params};
//				gsl_monte_vegas_init(s);

//				gsl_monte_vegas_integrate (&F, low, hi, 1, 100, rng_handle, s, &res, &err);
//			
//				dummy_table[A_max_it][A_mean_it][sd_it]=res;
//			}
//		}
//	}

//	gsl_monte_vegas_free(s);

//	return dummy_table;
//}

