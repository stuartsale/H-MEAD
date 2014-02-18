#include "LF.h"


LF::LF(string filename, float feh_in)
{
	beta=1.;

	feh=feh_in;
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
	l=PI;
	b=0;
	metal_prob=1;

	gauss_table[0]=0.;
	gauss_table[800]=1.;
	for (int it_g=1; it_g<800;it_g++)
	{
		gauss_table[it_g]=gsl_cdf_gaussian_P(it_g/100.-4.,1.);
	}
}

void LF::set_prior_lf(float l_in, float b_in, float dl, float db)
{
	l=l_in;
	b=b_in;

	cout << l << " " << b << endl;
	norm_lf=0;
	prior_lf.clear();
	for (int it=0; it<150; it++)
	{
		prior_lf.push_back(exp(log_prior_LF(5*log10(it*100.+50.)-5., feh, l*PI/180, b*PI/180))*100*sin(dl*PI/180.)*sin(db*PI/180.)*0.00679);	// 0.00679 = K$ or earlier stars pc^-3
		norm_lf+=prior_lf[it];
	}
	
	alpha=norm_lf/beta;
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
	float prob1=0;

	for (int it2=0; it2<LF_vec.size(); it2+=1)
	{
		for (int it=0; it<A_rel[0].size(); it++)	// run though A(d)
		{
			prob1+=(gauss_table[max(min(800, int(((A_max[it2][it]-A_rel[2][it])/A_rel[3][it] +4)*100)),0)] 
				- gauss_table[max(min(800, int(((A_min[it2][it]-A_rel[2][it])/A_rel[3][it] +4)*100)),0)])*prior_lf[it]*LF_vec[it2][1] ;
		}
	}	


	return prob1/norm_lf;


}	
