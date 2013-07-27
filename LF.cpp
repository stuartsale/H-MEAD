#include "LF.h"


LF::LF(string filename)
{

	ifstream input1;
	input1.open(filename.c_str());
	if(!input1) { //output file couldn't be opened
		cerr << "Error: file " << filename << " could not be opened" << endl;
	}

	while (!input1.eof())				// Running down file - reading it in
	{
		string str; 	
		getline(input1, str);			// grab a line
		string temp;
		stringstream sstest(str);
		sstest>>temp;
		if (temp!="#")				// check the line isn't commented out
		{

			float buf;
			stringstream ss(str);		// turn that line into a stringstream
		
			vector<float> fromfile;	//vector to put contents of line into
		
			while (ss>>buf)
			{			// Includes implicit conversion from string to float
				fromfile.push_back(buf);	
			}
			if (fromfile.size()==2)		// check there's something in the line
			{
				LF_vec.push_back(fromfile);
			}
		}
	}

	metal_prob=1;
	
	norm_lf=0;
	for (int it=0; it<150; it++)
	{
		prior_lf.push_back(exp(log_prior(5*log10(it*100.+50.)-5., feh, PI, 0.)));
		norm_lf+=prior_lf[it];
	}
	cout << prior_lf.size() << " " << prior_lf[10] << " " << norm_lf << endl;
	

}

float LF::LF_prob_last(vector <bin_obj> A_rel)
{
	float prob=0;
	//float norm=0;
	//float prior;
	float A_max;

	float dist_mod_lf;

//	for (int it=0; it<A_rel.size(); it++)	// run though A(d)
//	{
//		dist_mod_lf=5*log10(it*100.+50.)-5;
//		//prior=exp(log_prior(5*log10(it*100.+50.)-5., feh, PI, 0.));
//		//norm+=prior;
//		for (int it2=0; it2<LF_vec.size(); it2++)
//		{
//			if (LF_vec[it2][0]+dist_mod_lf+1.02*A_rel[it][0]>=r_min)
//			{
//				A_max=(r_max-dist_mod_lf-LF_vec[it2][0])/1.02;
//				//if (LF_vec[it2][0]+5*log10(it*100.+50.)-5.+1.02*A_rel[it][0]<r_max)
//				if (A_max>10)
//				{
//					prob+=prior_lf[it]/*lookup_table[0][0][0]*/*LF_vec[it2][1];
//				}
//				else if (A_max>0 && A_rel[it][0]<10)
//				{
//					prob+=prior_lf[it]*lookup_table[int(A_max*10.)][int(A_rel[it][0]*10.)][int(A_rel[it][1]*10.)]*LF_vec[it2][1];
//				}
//				else {break;}
//			}
//	
//		}
//	}
//	//cout << log(prob/norm_lf) << " " << prob << " " << norm_lf << " " << prior_lf[10] << " " << A_rel.size() << endl;

	for (int it=0; it<A_rel.size(); it++)	// run though A(d)
	{
		dist_mod_lf=5*log10(it*100.+50.)-5;
		//prior=exp(log_prior(5*log10(it*100.+50.)-5., feh, PI, 0.));
		//norm+=prior;
		for (int it2=0; it2<LF_vec.size(); it2++)
		{
			if (LF_vec[it2][0]+dist_mod_lf+0.276*A_rel[it].last_mean_A>=J_min)
			{
				A_max=(r_max-dist_mod_lf-LF_vec[it2][0])/0.276;
				//if (LF_vec[it2][0]+5*log10(it*100.+50.)-5.+1.02*A_rel[it].last_mean_A<r_max)
				if (A_max>10)
				{
					prob+=prior_lf[it]/*lookup_table[0][0][0]*/*LF_vec[it2][1];
				}
				else if (A_max>0 && A_rel[it].last_mean_A<10)
				{
					prob+=prior_lf[it]*lookup_table[int(A_max*10.)][int(A_rel[it].last_mean_A*10.)][int(A_rel[it].last_sd_A*10.)]*LF_vec[it2][1];
				}
				else {break;}
			}
	
		}
	}

	return log(prob/norm_lf);


}	

float LF::LF_prob_test(vector <bin_obj> A_rel)
{
	float prob=0;

	float A_max;

	float dist_mod_lf;


	for (int it=0; it<A_rel.size(); it++)	// run though A(d)
	{
		dist_mod_lf=5*log10(it*100.+50.)-5;
		//prior=exp(log_prior(5*log10(it*100.+50.)-5., feh, PI, 0.));
		//norm+=prior;
		for (int it2=0; it2<LF_vec.size(); it2++)
		{
			if (LF_vec[it2][0]+dist_mod_lf+0.276*A_rel[it].test_mean_A>=J_min)
			{
				A_max=(r_max-dist_mod_lf-LF_vec[it2][0])/0.276;
				//if (LF_vec[it2][0]+5*log10(it*100.+50.)-5.+1.02*A_rel[it].test_mean_A<r_max)
				if (A_max>10)
				{
					prob+=prior_lf[it]/*lookup_table[0][0][0]*/*LF_vec[it2][1];
				}
				else if (A_max>0 && A_rel[it].test_mean_A<10)
				{
					prob+=prior_lf[it]*lookup_table[int(A_max*10.)][int(A_rel[it].test_mean_A*10.)][int(A_rel[it].test_sd_A*10.)]*LF_vec[it2][1];
				}
				else {break;}
			}
	
		}
	}

	return log(prob/norm_lf);


}
