#include <time.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

//VariLogNormal VLN;

int main (void)
{

// Setting up random number generators

//	srand((unsigned)time(0));
//	double s=double(rand())/(double(RAND_MAX)+1.0);	//Generating a random seed for the MT generator in range 0-1
	//MT urng(s);
//	MultWithCarry urng(s);
//	Random::Set(urng);
//	srand ( time(NULL) );

// testing newran03

/*	clock_t start = clock();

	for (int i=0; i<100000; i++){cout << VLN.Next(exp(1.5),sqrt((exp(1)-1)*exp(3)))<<endl;}

	//cout << "MCMC_bin time = " << (clock()-start)/((double)CLOCKS_PER_SEC) << endl;

*/
	//initialize up random number generator
    	gsl_rng_env_setup ();
	gsl_rng* rng_handle = gsl_rng_alloc (gsl_rng_taus2);

	//seed the random no generator
	gsl_rng_set(rng_handle, time(0));
//	cout << s << " " << " " << int(s*pow(2,30)) << " " << gsl_ran_lognormal(rng_handle,1 ,1) <<endl;

	//alas, here is a number output!
//	for (int i=0; i<100000; i++){cout << gsl_ran_lognormal(rng_handle,1,1)<<endl;}

//	cout << "MCMC_bin time = " << (clock()-start)/((double)CLOCKS_PER_SEC) << endl;


	for (int j=0; j<1000; j++)
	{
		double sum=0, num=0;
		double mean=0.01, new_mean, theta=0.0002, k=100, sigma=0.02, sd=0.002;
		double trans_prob, prior=0;
		while(num<10000)
		{
		//	new_mean=gsl_ran_gamma(rng_handle,mean/theta,theta);
		//	new_mean=gsl_ran_gamma(rng_handle, k, mean/k);
			new_mean=gsl_ran_lognormal(rng_handle,log(mean)-pow(sigma,2)/2,sigma);
			if (new_mean!=0)
			{
		//	cout << mean << "  ttt " << new_mean << endl;
			// new to old
		//	trans_prob=log(gsl_ran_gamma_pdf(mean, new_mean/theta, theta));
		//	trans_prob=log(gsl_ran_gamma_pdf(mean,k, new_mean/k));
			trans_prob=log(gsl_ran_lognormal_pdf(mean, log(new_mean)-pow(sigma,2)/2 ,sigma));

			// old to new
		//	trans_prob-=log(gsl_ran_gamma_pdf(new_mean, mean/theta, theta));
		//	trans_prob-=log(gsl_ran_gamma_pdf(new_mean,k, mean/k));
			trans_prob-=log(gsl_ran_lognormal_pdf(new_mean, log(mean)-pow(sigma,2)/2 ,sigma));

			prior=2*log(mean/new_mean);//+(pow(log(mean),2)+pow(log(new_mean),2))/(2*sigma*sigma);

		//	cout << mean << "ttt   " << trans_prob << " " << prior << " " << log(gsl_ran_lognormal_pdf(mean, 0 ,sigma)) << endl;

			if (gsl_ran_flat(rng_handle, 0., 1.)<exp(trans_prob+prior)){mean=new_mean;sum+=mean;num++;}

		//	mean=new_mean;
			}
		}
		cout << sum/num << " " << mean <<  " " << num << endl;
	}

	

	return 0;
}
