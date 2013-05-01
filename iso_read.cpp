#include "iso_read.h"


vector<iso_obj> iso_read(const string &filename)		// Function to read in calibration data
{
	ifstream input1;
	input1.open(filename.c_str());
	if(!input1) { //output file couldn't be opened
		cerr << "Error: file could not be opened" << endl;
		exit(1);
	}

	vector<iso_obj> totalfile;		// making a vector to store the calibration_obj
	
	while (!input1.eof())				// Running down file - reading it in
	{
		string str; 	
		getline(input1, str);			// grab a line
		string temp;
		stringstream sstest(str);
		sstest>>temp;
		if (temp!="#")				// check the line isn't commented out
		{		
			double buf;
			stringstream ss(str);		// turn that line into a stringstream
		
			vector<double> fromfile;	//vector to put contents of line into
		
			while (ss>>buf){			// Includes implicit conversion from string to double
				fromfile.push_back(buf);	
			}
			if (fromfile.size()==8)		// check there's something in the line
			{
         			iso_obj objnew(fromfile[0], fromfile[1], fromfile[2], fromfile[3], fromfile[4], fromfile[5], fromfile[6], fromfile[7],1, "IPHAS");
				totalfile.push_back(objnew);
			}		
		}
	}
	return totalfile;
	input1.close();
}

vector<iso_obj> iso_read_Tg(const string &filename)		// Function to read in calibration data
{
	ifstream input1;
	input1.open(filename.c_str());
	if(!input1) { //output file couldn't be opened
		cerr << "Error: file could not be opened" << endl;
		exit(1);
	}

	vector<iso_obj> totalfile;		// making a vector to store the calibration_obj
	
	while (!input1.eof())				// Running down file - reading it in
	{
		double jac;
		string str; 	
		getline(input1, str);			// grab a line
		string temp;
		stringstream sstest(str);
		sstest>>temp;
		if (temp!="#")				// check the line isn't commented out
		{		
			double buf;
			stringstream ss(str);		// turn that line into a stringstream
		
			vector<double> fromfile;	//vector to put contents of line into
		
			while (ss>>buf){			// Includes implicit conversion from string to double
				fromfile.push_back(buf);	
			}
			if (fromfile.size()==14)		// check there's something in the line
			{
				if (fromfile[5]!=0 && fromfile[6]!=0 && fromfile[7]!=0 && fromfile[8]!=0)
				{
					jac=abs(fromfile[5]*fromfile[8]-fromfile[6]*fromfile[7]);
				}
				else {jac=0;}
         			iso_obj objnew(fromfile[0], fromfile[1], fromfile[2], fromfile[3], fromfile[4], fromfile[9], fromfile[10], fromfile[11], jac, "IPHAS");
				totalfile.push_back(objnew);
			}		
		}
	}
	return totalfile;
	input1.close();
}

vector<iso_obj> iso_read_Tg_2MASS(const string &filename)		// Function to read in calibration data
{
	ifstream input1;
	input1.open(filename.c_str());
	if(!input1) { //output file couldn't be opened
		cerr << "Error: file could not be opened" << endl;
		exit(1);
	}

	vector<iso_obj> totalfile;		// making a vector to store the calibration_obj
	
	while (!input1.eof())				// Running down file - reading it in
	{
		double jac;
		string str; 	
		getline(input1, str);			// grab a line
		string temp;
		stringstream sstest(str);
		sstest>>temp;
		if (temp!="#")				// check the line isn't commented out
		{		
			double buf;
			stringstream ss(str);		// turn that line into a stringstream
		
			vector<double> fromfile;	//vector to put contents of line into
		
			while (ss>>buf){			// Includes implicit conversion from string to double
				fromfile.push_back(buf);	
			}
			if (fromfile.size()==11)		// check there's something in the line
			{
				iso_obj objnew(fromfile[0], fromfile[1], fromfile[2], fromfile[3], fromfile[4], fromfile[6], fromfile[7], fromfile[8], fromfile[5], "2MASS");
				totalfile.push_back(objnew);
			}
		}
	}
	cout << "iso length:" << totalfile.size() << endl;
	return totalfile;
	input1.close();
}
