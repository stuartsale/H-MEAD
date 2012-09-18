#include "iso_obj.h"

// constructor
iso_obj::iso_obj(double feh_in, double Mi_in, double age_in, double logT_in, double logg_in, double r0_in, double i0_in, double ha0_in, double Jac_in)
{
	Mi=Mi_in;
	logAge=age_in;
	logT=logT_in;
	logg=logg_in;
	feh=feh_in;

	r0=r0_in;
	i0=i0_in;
	ha0=ha0_in;

	Jac=Jac_in;
	
   
/*	a=((0.00051868526*(r0-i0)-0.0029977593)-(0.0001234*(r0-i0)-0.001616)) /1.08573;			// 2nd order on r-i	}
	b=((-0.03202578*(r0-i0)+1.1109227)-(-0.014138594*(r0-i0)+0.8480197)) /1.08573;			// 1st order on r-i	}A(r)-A(i)
	c=((0.0028546092*(r0-i0)+0.0048040543)-(0.000120494*(r0-i0)+0.0016170746)) /1.08573;		// 0th order on r-i	}

	l=((0.00051868526*(r0-i0)-0.0029977593)-(0.000021304992*(r0-i0)-0.000029709503)) /1.08573;	// 2nd order on r-ha	}
	m=((-0.03202578*(r0-i0)+1.1109227)-(0.00016239351*(r0-i0)+1.0192878)) /1.08573;			// 1st order on r-ha	}A(r)-A(Ha)
	n=((0.0028546092*(r0-i0)+0.0048040543)-(-0.0007983722*(r0-i0)+0.0004974765)) /1.08573;		// 0th order on r-ha	}
*/

//R=3.1

	u=(0.0007075*(r0-i0)-0.00504174) /1.08573;						// 2nd order on r 				
	v=(-0.0396596*(r0-i0)+1.1106105) /1.08573;							// 1st order on r
	w=(0.001023008*(r0-i0)+0.0008454) /1.08573;							// 0th order on r	

	u_i=(0.00004671*(r0-i0)-0.0025414) /1.08573;
	v_i=(-0.0169457*(r0-i0)+0.7929863) /1.08573;
	w_i=(0.0011835*(r0-i0)-0.0037645) /1.08573;

	u_ha=(0.00007665*(r0-i0)-0.00007222) /1.08573;
	v_ha=(-0.00054454*(r0-i0)+1.003924) /1.08573;
	w_ha=(-0.0004797*(r0-i0)+0.00028745) /1.08573;


/*//R=2.9

	u=(0.000825*(r0-i0)-0.005582) /1.08573;							// 2nd order on r 				
	v=(-0.04209*(r0-i0)+1.11201874) /1.08573;							// 1st order on r
	w=(0.00115325*(r0-i0)+0.0010913316) /1.08573;							// 0th order on r	

	u_i=(-0.0000659*(r0-i0)-0.002717) /1.08573;
	v_i=(-0.017621*(r0-i0)+0.77944) /1.08573;
	w_i=(0.001279*(r0-i0)-0.0003786) /1.08573;

	u_ha=(0.00007835*(r0-i0)-0.00007595) /1.08573;
	v_ha=(-0.0005294*(r0-i0)+0.99945) /1.08573;
	w_ha=(-0.0005279*(r0-i0)+0.0003166) /1.08573; */



}

iso_obj::iso_obj(void)
{
	Mi=1;
	logAge=9;
	logT=3.75;
	logg=4.5;
	feh=0;

	r0=0;
	i0=0;
	ha0=0;

	Jac=0;
}

double iso_obj::IMF(void)
{
	return pow(Mi, -2.7);
}

double iso_obj::redline(double r_i1)
{
	double A_int;
	A_int=0.;// (-(v-v_i)+sqrt(pow(v-v_i,2)-4*(u-u_i)*((w-w_i)+(r0-i0)-r_i1)))/(2*(u-u_i));  //(u-u_i, v-v_i, (w-w_i)+(r0-i0)-r_i1, +1);
	A_int=(-(r0-i0)+r_i1)/(v-v_i);
	return (u-u_ha)*pow(A_int,2) + (v-v_ha)*A_int + (w-w_ha)+(r0-ha0);
}







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

	for (int it=0; it<A_rel.size(); it++)	// run though A(d)
	{
		norm+=
		for (int it2=0; it2<LF_vec.size(); it2++)
		{
			prob+=


}	

	
