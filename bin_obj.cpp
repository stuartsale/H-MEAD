#include "bin_obj.h"
#include "iphas_obj.h"


bin_obj::bin_obj(void)
{
	last_mean_A=0.;
	last_sd_A=0.;
	last_mu=0.;
	last_sigma=0.;

	last_mean_rho=0.;

	last_n=0;
	last_lnA_sum=0.;
	last_lnA2_sum=0.;
	test_lnA_sum=0.; 
	test_lnA2_sum=0.;

}

void bin_obj::initial_add(iphas_obj * star)
{
	last_n++;
	last_lnA_sum+=log(star->last_A);
	last_lnA2_sum+=pow(log(star->last_A),2.);
}

void bin_obj::try_add(iphas_obj * star)
{
	test_n=last_n+1;
//	if (dist_bin==130){cout << "+ " << last_n << " " << last_lnA_sum << " " << log(star->test_A) << endl;}
	test_lnA_sum=last_lnA_sum+log(star->test_A);
	test_lnA2_sum=last_lnA2_sum+pow(log(star->test_A),2.);

	set_test_prob();
}

void bin_obj::try_update(iphas_obj * star)
{
//	if (dist_bin==130){cout << "x " << last_n << " " << last_lnA_sum << " " << log(star->last_A)-log(star->test_A) << endl;}
	test_lnA_sum=last_lnA_sum+log(star->test_A) - log(star->last_A); 
	test_lnA2_sum=last_lnA2_sum+pow(log(star->test_A),2.) - pow(log(star->last_A),2.);

	set_test_prob();
}

void bin_obj::try_remove(iphas_obj * star)
{
	test_n=last_n-1;
//	if (dist_bin==130){cout << "- " << last_n << " " << last_lnA_sum << " " << log(star->last_A) << endl;}
	test_lnA_sum=last_lnA_sum-log(star->last_A);
	test_lnA2_sum=last_lnA2_sum-pow(log(star->last_A),2.);

	set_test_prob();
}

void bin_obj::initial_rho_to_A(float near_A)
{
	last_mean_A=last_mean_rho+near_A;
	last_sd_A=0.4;//*last_mean_A;
	last_sigma=sqrt(log(1+pow(last_sd_A/last_mean_A,2.) ));
	last_mu=log(last_mean_A)-pow(last_sigma,2);
}

void bin_obj::rho_to_A(float near_A)
{
	test_mean_A=test_mean_rho+near_A;
	test_sd_A=0.4;//*test_mean_A;
	test_sigma=sqrt(log(1+pow(test_sd_A/test_mean_A,2.) ));
	test_mu=log(test_mean_A)-pow(test_sigma,2)/2.;
}

void bin_obj::set_last_prob(void)
{
	last_prob=-last_n*log(last_sigma) - last_lnA_sum -(last_lnA2_sum - 2*last_mu*last_lnA_sum + last_n*pow(last_mu,2) )/(2*pow(last_sigma,2) ) ;
}

void bin_obj::set_test_prob(void)
{
	test_prob=-test_n*log(test_sigma) - test_lnA_sum -(test_lnA2_sum - 2*test_mu*test_lnA_sum + test_n*pow(test_mu,2) )/(2*pow(test_sigma,2) ) ;
}

void bin_obj::accept(void)
{
	last_n=test_n;
	last_lnA_sum=test_lnA_sum;
	last_lnA2_sum=test_lnA2_sum;

	last_mean_A=test_mean_A;
	last_sd_A=test_sd_A;
	last_mu=test_mu;
	last_sigma=test_sigma;
	last_mean_rho=test_mean_rho;
	
	last_prob=test_prob;
}

void bin_obj::reject(void)
{
	test_n=last_n;
	test_lnA_sum=last_lnA_sum;
	test_lnA2_sum=last_lnA2_sum;

	test_mean_A=last_mean_A;
	test_sd_A=last_sd_A;
	test_mu=last_mu;
	test_sigma=last_sigma;
	test_mean_rho=last_mean_rho;
	
	test_prob=last_prob;
}

void bin_obj::chain_push_back(void)
{
	A_chain.push_back(last_mean_A);
	sd_chain.push_back(last_sd_A);
	rho_chain.push_back(last_mean_rho);
}

void bin_obj::mean_intervals(void)
{
	float A_sum=0., sd_sum=0., rho_sum=0.;
	for (int m=floor(0.50*A_chain.size()); m<A_chain.size(); m++)
	{
		A_sum+=A_chain[m];
		sd_sum+=sd_chain[m];
		rho_sum+=rho_chain[m];
	}
	final_A=A_sum/ceil(0.50*A_chain.size());
	final_sd=sd_sum/ceil(0.50*sd_chain.size());
	final_rho=rho_sum/ceil(0.50*rho_chain.size());

	vector <float> A_diffs, sd_diffs, rho_diffs;
	for (int m=floor(0.50*A_chain.size()); m<A_chain.size(); m++)
	{
		A_diffs.push_back(abs(A_chain[m]-final_A));
		sd_diffs.push_back(abs(sd_chain[m]-final_sd));
		rho_diffs.push_back(abs(rho_chain[m]-final_rho));
	}
	sort(A_diffs.begin(), A_diffs.end());
	sort(sd_diffs.begin(), sd_diffs.end());
	sort(rho_diffs.begin(), rho_diffs.end());

	final_dA=A_diffs[int(0.682*A_diffs.size())];
	final_dsd=sd_diffs[int(0.682*sd_diffs.size())];
	final_drho=rho_diffs[int(0.682*rho_diffs.size())];

	//A_mean[it].error_measure=sqrt(pow(A_mean[it].d_mean/A_mean[it].mean_A,2)+pow(A_mean[it].d_sigma/A_mean[it].sigma,2));
}

