#include "bin_obj.h"
#include "iphas_obj.h"


bin_obj::bin_obj(void)
{
	mean_A=0;
	sigma=0;

	last_n=0;
	last_A_sum=0.;
	last_A2_sum=0.;
	test_A_sum=0.; 
	test_A2_sum=0.;

}

void bin_obj::initial_add(iphas_obj * star)
{
	last_n++;
	last_A_sum+=star->last_A;
	last_A2_sum+=pow(star->last_A,2.);
}

void bin_obj::try_add(iphas_obj * star)
{
	test_n=last_n+1;
	test_A_sum=last_A_sum+star->test_A;
	test_A2_sum=last_A2_sum+pow(star->test_A,2.);
}

void bin_obj::try_remove(iphas_obj * star)
{
	test_n=last_n-1;
	test_A_sum=last_A_sum-star->test_A;
	test_A2_sum=last_A2_sum-pow(star->test_A,2.);
}

void bin_obj::accept(void)
{
	last_n=test_n;
	last_A_sum=test_A_sum;
	last_A2_sum=test_A2_sum;
}


//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

bin_obj2::bin_obj2(void)
{
	mean_A=0;
	d_mean=0;
	sigma=0;
	d_sigma=0;
	size=0;
	error_measure=0;
	sum=0;
	diff=0;
	d_diff=0;
}

