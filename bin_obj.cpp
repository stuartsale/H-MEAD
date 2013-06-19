#include "bin_obj.h"


bin_obj::bin_obj(void)
{
	mean_A=0;
	sigma=0;

	last_n=0;
	
	test_star_list.reserve(1000);
	last_star_list.reserve(1000);
}

void bin_obj::initial_add(iphas_obj * star)
{
	last_star_list.push_back(star);
	last_n++;
}

void bin_obj::try_add(iphas_obj * star)
{
	test_star_list=last_star_list;
	test_star_list.push_back(star);
	test_n=last_n+1;
}

void bin_obj::try_remove(iphas_obj * star)
{
	test_star_list=last_star_list;
	it_found=find (test_star_list.begin(), test_star_list.end(), star);
	test_star_list.erase(it_found);
	test_n=last_n-1;
}

void bin_obj::accept(void)
{
	//last_star_list=test_star_list;
	last_n=test_n;
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

