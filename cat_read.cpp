#include "cat_read.h"

// function to read in IPHAS data, works in much the same manner as calibration_read
vector<iphas_obj> iphas_read(string filename,float &r_min1,float &i_min1,float &ha_min1,float &r_max1, float &i_max1, float &ha_max1)		
{						
   //-----------------------
   // Pre-cond: filename is a string containing the filename of the IPHAS 
   // catalogue data for the region examined.
   //--------------------
   // Post-cond: iphas_colours is returned. As a vector of iphas_obj, it contains all sources in the IPHAS catalogue
   // iff they match these conditions: (a) Classified as stellar or _probably_ stellar  in all three bands (Ha, r' and i')
   //----------------------------------------------

	vector<iphas_obj> iphas_colours;
	ifstream iphas_data;
	iphas_data.open(filename.c_str());
	if(!iphas_data) { //output file couldn't be opened
		cerr << "Error: IPHAS catalogue file could not be opened \n";
		exit(1);
	}	
	while (!iphas_data.eof())				// Running down file
	{
		string str1; 	
		getline(iphas_data, str1);
		string temp;
		stringstream sstest(str1);
		sstest>>temp;
		if (temp!="#")
		{		
			float buffer;
			stringstream ss1(str1);
		
			vector<float> infromfile;
		
			while (ss1>>buffer){			// Includes implicit conversion from string to float
				infromfile.push_back(buffer);
			}// correct length		-r_class is stellar or prob stellar-	--i_class is stellar or prob stellar-----	-ha class is stellar or prob stellar---		--------r,i,Ha photometry non-zero----------------		-----r,i,Ha photometric errors non-zero-------------------		---------small RA & DEC offsets between r and i, and r and Ha-------------------
			if (infromfile.size()==30 && (infromfile[6]==-1 || infromfile[6]==-2) && (infromfile[11]==-1 || infromfile[11]==-2) && (infromfile[16]==-1 || infromfile[16]==-2) && (abs(infromfile[4])>0.1) && (abs(infromfile[9])>0.1) && (abs(infromfile[14])>0.1) && (infromfile[5]!=0) && (infromfile[10]!=0) && (infromfile[15]!=0) && infromfile[18]<=1.0 && infromfile[19]<=1.0 && infromfile[20]<=1.0 && infromfile[21]<=1.0) 	//selecting only stellar or probably stellar objects and those with small RA & DEC offsets
			{

//   		        	if(infromfile[4] > r_max1) { r_max1 = infromfile[4];}
//            			if(infromfile[9] > i_max1) { i_max1 = infromfile[9];}
//            			if(infromfile[14] > ha_max1) { ha_max1 = infromfile[14];}
			        iphas_obj next_obj(infromfile[4], infromfile[9], infromfile[14], sqrt(pow(infromfile[5],2)+pow(0.016165105,2)), sqrt(pow(infromfile[10],2)+pow(0.016165105,2)), sqrt(pow(infromfile[15],2)+pow(0.016165105,2)), infromfile[0], infromfile[1], infromfile[22], infromfile[23]*1.0857, infromfile[28], infromfile[27], infromfile[29]);	// making the new iphas_obj

				iphas_colours.push_back(next_obj);																											// pushing it into the vetor
			}

			else if (infromfile.size()==22 && (infromfile[6]==-1 || infromfile[6]==-2) && (infromfile[11]==-1 || infromfile[11]==-2) && (infromfile[16]==-1 || infromfile[16]==-2) && (infromfile[4]!=0) && (infromfile[9]!=0) && (infromfile[14]!=0) && (infromfile[5]!=0) && (infromfile[10]!=0) && (infromfile[15]!=0) && infromfile[18]<=1.0 && infromfile[19]<=1.0 && infromfile[20]<=1.0 && infromfile[21]<=1.0) 	//selecting only stellar or probably stellar objects and those with small RA & DEC offsets
			{

//   		        	if(infromfile[4] > r_max1) { r_max1 = infromfile[4];}
//            			if(infromfile[9] > i_max1) { i_max1 = infromfile[9];}
//            			if(infromfile[14] > ha_max1) { ha_max1 = infromfile[14];}
			        iphas_obj next_obj(infromfile[4], infromfile[9], infromfile[14], sqrt(pow(infromfile[5],2)+pow(0.016165105,2)), sqrt(pow(infromfile[10],2)+pow(0.016165105,2)), sqrt(pow(infromfile[15],2)+pow(0.016165105,2)), infromfile[0], infromfile[1]);	// making the new iphas_obj

				iphas_colours.push_back(next_obj);																											// pushing it into the vetor
			}

			else if (infromfile.size()==8 && infromfile[2]>0 && infromfile[3]>0 && infromfile[4]>0 && infromfile[5]<0.16 && infromfile[6]<0.16 && infromfile[7]<0.16)
			{
			        iphas_obj next_obj(infromfile[2], infromfile[3], infromfile[4], sqrt(pow(infromfile[5],2)+pow(0.016165105,2)), sqrt(pow(infromfile[6],2)+pow(0.016165105,2)), sqrt(pow(infromfile[7],2)+pow(0.016165105,2)), infromfile[0], infromfile[1]);	// making the new iphas_obj

				iphas_colours.push_back(next_obj);
			} 

			else if (infromfile.size()==22)		// if sources are saturated alter bright limits
			{
				if (infromfile[6]==-9 && infromfile[4]>r_min1){r_min1=infromfile[4];}		//r'
				if (infromfile[11]==-9 && infromfile[9]>i_min1){i_min1=infromfile[9];}		//i'
				if (infromfile[16]==-9 && infromfile[14]>ha_min1){ha_min1=infromfile[14];} 	//Ha
			}
		}	
	}
	iphas_data.close();
	return  iphas_colours;
}


// Function to read in 2MASS data, based on IPHAS read
vector<iphas_obj> TWOMASS_read(string filename,float &J_min1,float &H_min1,float &K_min1,float &J_max1, float &H_max1, float &K_max1)
{
	vector<iphas_obj> TWOMASS_colours;
	ifstream TWOMASS_data;
	TWOMASS_data.open(filename.c_str());
	if(!TWOMASS_data) { //output file couldn't be opened
		cerr << "Error: 2MASS catalogue file could not be opened \n";
		exit(1);
	}
	while (!TWOMASS_data.eof())				// Running down file
	{
		string str1; 	
		getline(TWOMASS_data, str1);
		string temp;
		stringstream sstest(str1);
		sstest>>temp;
		if (temp!="#")
		{		
			float buffer;
			stringstream ss1(str1);
		
			vector<float> infromfile;
		
			while (ss1>>buffer){			// Includes implicit conversion from string to float
				infromfile.push_back(buffer);	
				if (infromfile.size()==9)
				{
					iphas_obj new_obj(infromfile[3], infromfile[5], infromfile[7], infromfile[4], infromfile[6], infromfile[7], infromfile[0], infromfile[1], "2MASS");
					TWOMASS_colours.push_back(new_obj);
				}
			}

		}
	}	
	TWOMASS_data.close();
	return  TWOMASS_colours;
}
