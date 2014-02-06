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
	feh=0;
	
	norm_lf=0;
	for (int it=0; it<150; it++)
	{
		prior_lf.push_back(exp(log_prior(5*log10(it*100.+50.)-5., feh, PI, 0.)));
		norm_lf+=prior_lf[it];
	}
}

void LF::precompute_Aminmax(void)
{
	float A_max1, A_min1;
	float dist_mod_lf;

	vector <vector <float> > A_max2, A_min2;

	A_max.resize(LF_vec.size(), vector <float> (150));
	A_min.resize(LF_vec.size(), vector <float> (150));
	A_max2.resize(LF_vec.size(), vector <float> (150));
	A_min2.resize(LF_vec.size(), vector <float> (150));

	for (int it=0; it<150; it++)	// run though A(d)
	{	
		dist_mod_lf=5*log10(it*100.+50.)-5;
		for (int it2=0; it2<LF_vec.size(); it2+=1)
		{
			A_max1=(r_max-dist_mod_lf-LF_vec[it2][0])/0.838;
			if (A_max1>0){A_max2[it2][it]=A_max1;}
			else {A_max2[it2][it]=1E-9;}
			
			A_min1=(r_min-dist_mod_lf-LF_vec[it2][0])/0.838;
			if (A_min1>0){A_min[it2][it]=A_min1;}
			else {A_min2[it2][it]=1E-9;}
		}
	}
	vmlSetMode( VML_EP );
	for (int it2=0; it2<LF_vec.size(); it2+=1)
	{
		vsLn(150, &A_min2[it2][0], &A_min[it2][0]);
		vsLn(150, &A_max2[it2][0], &A_max[it2][0]);
	}
}


float LF::LF_prob2(vector < vector <float> > A_rel)
{
	float prob=0;
	float prob1=0;
	float prob2=0;

	size_t rel_length=A_rel[0].size();
	float prob_array[rel_length];
	float part_array1_max[rel_length];
	float part_array2_max[rel_length];
	float part_array3_max[rel_length];
	float part_array4_max[rel_length];
	float part_array1_min[rel_length];
	float part_array2_min[rel_length];
	float part_array3_min[rel_length];
	float part_array4_min[rel_length];
	float part_array5[rel_length];

	float dist_mod_lf;
//	float dA[150];
//	float aa[150];
//	for (int i=0; i<150; i++){aa[i]=A_rel[0][i]-10.;}

//	vsLn(150, aa, dA);
//	cout << dA[0] << endl;

//	for (int it=0; it<A_rel[0].size(); it++)	// run though A(d)
//	{
//		for (int it2=0; it2<LF_vec.size(); it2+=1)
//		{
//			if (A_max[it2][it]>0 && A_min[it2][it]>0)
//			{
//				prob+=1*prior_lf[it]*(gsl_cdf_lognormal_P(A_max[it2][it], A_rel[2][it], A_rel[3][it])-gsl_cdf_lognormal_P(A_min[it2][it], A_rel[2][it], A_rel[3][it]))*LF_vec[it2][1];
//			}
//			else if (A_max[it2][it]>0)
//			{
//				prob+=1*prior_lf[it]*gsl_cdf_lognormal_P(A_max[it2][it], A_rel[2][it], A_rel[3][it])*LF_vec[it2][1];
//			}
//			else {break;}
//	
//		}
//	}


	for (int it2=0; it2<LF_vec.size(); it2+=1)
	{
		vsSub(rel_length, &A_max[it2][0], &A_rel[2][0], part_array2_max);
		vsDiv(rel_length, part_array2_max, &A_rel[3][0], part_array3_max);
		vsCdfNorm(rel_length, part_array3_max, part_array4_max);

		vsSub(rel_length, &A_min[it2][0], &A_rel[2][0], part_array2_min);
		vsDiv(rel_length, part_array2_min, &A_rel[3][0], part_array3_min);
		vsCdfNorm(rel_length, part_array3_min, part_array4_min);

		vsSub(rel_length, part_array4_max, part_array4_min, part_array5);
		vsMul(rel_length, part_array5, &prior_lf[0], prob_array);

		for (int i =0; i<rel_length; i++){prob2+=prob_array[i];}
		prob2*=LF_vec[it2][1];
		prob1+=prob2;
	}	


	return log(prob1/norm_lf);


}	
