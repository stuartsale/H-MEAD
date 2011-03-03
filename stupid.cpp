#include big.h

vector <bin_obj2> big_func(vector <iphas_obj> stars, vector <bin_obj2> startA, vector<iso_obj> &isochrones, vector<iso_obj> &guess_set, double l, double b)
{
	int t=0;
	while (t<stars.size())
	{
		if (stars[t].r-stars[t].ha>guess_set[0].redline(stars[t].r-stars[t].i) || stars[t].r-stars[t].ha<guess_set[guess_set.size()-1].redline(stars[t].r-stars[t].i)){stars.erase(t);}
		else {t++;}
	}

// Make initial guess

	t=0;
	
	while (t<stars.size())
	{
		for (int t=0; t<guess_set.size();t++)
		{
			if (stars[t].r-stars[t].ha>guess_set[it].redline(stars[t].r-stars[t].i))
			{
				test_Mi=(((stars[t].r-stars[t].ha)-guess_set[it].redline(stars[t].r-stars[t].i))*guess_set[it-1].Mi + (guess_set[it-1].redline(stars[t]r-stars[t]i)-(stars[t].r-stars[t].ha))*guess_set[it].Mi)/(guess_set[it-1].redline(stars[t].r-stars[t].i)-guess_set[it].redline(stars[t].r-stars[t].i));
				try{test_logAge=max_age(test_Mi, isochrones)-log(2);} catch (int e){test_logAge=8.5;}
				test_iso=iso_get(0.,test_Mi, test_logAge, isochrones);
				Mi_sd = 0.025*sqrt(stars[t].d_r*stars[t].d_r+stars[t].d_ha*stars[t].d_ha)*(guess_set[it].Mi-guess_set[it-1].Mi)/(guess_set[it-1].redline(stars[t].r-stars[t].i)-guess_set[it].redline(stars[t].r-stars[t].i));
				break;
			}
		}
		test_A=quadratic(test_iso.u-test_iso.u_i, test_iso.v-test_iso.v_i, (test_iso.w-test_iso.w_i)+(test_iso.r0-test_iso.i0)-(stars[t].r-stars[t].i), +1);
		if (test_A<0){stars.erase(t);}
		else
		{
			test_dist_mod=stars[t].r-(test_iso.u*pow(test_A,2)+test_iso.v*test_A+test_iso.w)-test_iso.r0;
			test_dist=pow(10,test_dist_mod/5+1);
	
			stars[t].iso_obj_chain.push_back(test_iso);
			stars[t].dist_mod_chain.push_back(test_dist_mod);
			stars[t].A_chain.push_back(test_A);
		}
	}
