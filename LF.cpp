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
}

float LF::LF_prob(vector < vector <float> > A_rel)
{
	float prob=0;
	float norm=0;
	float prior;

	for (int it=0; it<A_rel.size(); it++)	// run though A(d)
	{
		prior=exp(log_prior(it*100.+50., feh, 180., 0.));
		norm+=prior;
		for (int it2=0; it2<LF_vec.size(); it2++)
		{
			prob+=prior*lookup_table[0][0][0]*LF_vec[it2][1];
		}
	}
	return log(prob/norm);


}	
