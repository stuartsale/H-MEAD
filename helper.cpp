#include "helper.h"

// solves a quadratic eqn, sign determines which of the two solutions is used
double quadratic(double a, double b, double c, int sign)	
{
	double det;
	det=pow(b,2) - 4 * a * c;
	if (det>=0)
	{
		if (sign>0){return (-b+sqrt(det))/(2 * a);}
		else if (sign<0){return (-b-sqrt(det))/(2 * a);}
	}
	else {throw 99;}
}


iso_obj iso_get(double targ_feh, double targ_Mi, double targ_logAge, vector<iso_obj> &isochrones)
{
	int feh_line1, feh_line2;
	int age_line1, age_line2;
	int mass_line;
	double mass_weight, feh_weight, age_weight;
	double r1=0, i1=0, ha1=0;

	if (targ_feh<-0.914 || targ_feh>0.417 || targ_Mi<0.15 || targ_Mi>19.14 || targ_logAge<7.05 || targ_logAge>10.25)
	{
		throw (5);
	}


	if (targ_feh>=-0.914 && targ_feh<-0.778){feh_line1=0; feh_line2=39000;feh_weight=(targ_feh+0.914)/0.136;}
	else if (targ_feh>=-0.778 && targ_feh<-0.652){feh_line1=39000; feh_line2=78000;feh_weight=(targ_feh+0.778)/0.126;}
	else if (targ_feh>=-0.652 && targ_feh<-0.472){feh_line1=78000; feh_line2=117000;feh_weight=(targ_feh+0.652)/0.180;}
	else if (targ_feh>=-0.472 && targ_feh<-0.343){feh_line1=117000; feh_line2=156000;feh_weight=(targ_feh+0.472)/0.129;}
	else if (targ_feh>=-0.343 && targ_feh<-0.243){feh_line1=156000; feh_line2=195000;feh_weight=(targ_feh+0.343)/0.100;}
	else if (targ_feh>=-0.243 && targ_feh<-0.160){feh_line1=195000; feh_line2=234000;feh_weight=(targ_feh+0.243)/0.083;}
	else if (targ_feh>=-0.160 && targ_feh<-0.090){feh_line1=234000; feh_line2=273000;feh_weight=(targ_feh+0.160)/0.070;}
	else if (targ_feh>=-0.090 && targ_feh<0.000){feh_line1=273000; feh_line2=312000;feh_weight=(targ_feh+0.090)/0.090;}
	else if (targ_feh>=0.000 && targ_feh<0.077){feh_line1=312000; feh_line2=351000;feh_weight=(targ_feh)/0.077;}
	else if (targ_feh>=0.077 && targ_feh<0.202){feh_line1=351000; feh_line2=390000;feh_weight=(targ_feh-0.077)/0.125;}
	else if (targ_feh>=0.202 && targ_feh<0.363){feh_line1=390000; feh_line2=429000;feh_weight=(targ_feh-0.202)/0.161;}
	else if (targ_feh>=0.363 && targ_feh<0.417){feh_line1=429000; feh_line2=468000;feh_weight=(targ_feh-0.363)/0.054;}
	else {throw(5);}



	age_line1=ceil((10.25-targ_logAge)/0.05)*600;
	age_line2=age_line1-600;
	age_weight=(targ_logAge-0.05*floor(targ_logAge/0.05))/0.05;

	if (isochrones[feh_line1+age_line1].Mi>targ_Mi || isochrones[feh_line1+age_line1+599].Mi<targ_Mi || isochrones[feh_line1+age_line2].Mi>targ_Mi || isochrones[feh_line1+age_line2+599].Mi<targ_Mi ||
	   isochrones[feh_line2+age_line1].Mi>targ_Mi || isochrones[feh_line2+age_line1+599].Mi<targ_Mi || isochrones[feh_line2+age_line2].Mi>targ_Mi || isochrones[feh_line2+age_line2+599].Mi<targ_Mi )
	{
		throw (6);
	}

	// 1,1
	if (targ_Mi<0.92*(isochrones[feh_line1+age_line1+599].Mi - isochrones[feh_line1+age_line1].Mi)+isochrones[feh_line1+age_line1].Mi)
	{
		mass_line=floor(300*(targ_Mi - isochrones[feh_line1+age_line1].Mi)/(0.92*(isochrones[feh_line1+age_line1+599].Mi - isochrones[feh_line1+age_line1].Mi)));
		mass_weight=300*(targ_Mi - isochrones[feh_line1+age_line1].Mi)/(0.92*(isochrones[feh_line1+age_line1+599].Mi - isochrones[feh_line1+age_line1].Mi)) - mass_line;
	}
	else
	{
		mass_line=300 + floor(300*(targ_Mi - isochrones[feh_line1+age_line1+300].Mi)/(0.08*(isochrones[feh_line1+age_line1+599].Mi - isochrones[feh_line1+age_line1].Mi)));
		mass_weight=300 + 300*(targ_Mi - isochrones[feh_line1+age_line1+300].Mi)/(0.08*(isochrones[feh_line1+age_line1+599].Mi - isochrones[feh_line1+age_line1].Mi)) - mass_line;
	}

	r1+=(1-feh_weight)*(1-age_weight)*(1-mass_weight) * isochrones[feh_line1+age_line1+mass_line].r0
	  + (1-feh_weight)*(1-age_weight)*(mass_weight) * isochrones[feh_line1+age_line1+mass_line+1].r0;
	i1+=(1-feh_weight)*(1-age_weight)*(1-mass_weight) * isochrones[feh_line1+age_line1+mass_line].i0
	  + (1-feh_weight)*(1-age_weight)*(mass_weight) * isochrones[feh_line1+age_line1+mass_line+1].i0;
	ha1+=(1-feh_weight)*(1-age_weight)*(1-mass_weight) * isochrones[feh_line1+age_line1+mass_line].ha0
	  + (1-feh_weight)*(1-age_weight)*(mass_weight) * isochrones[feh_line1+age_line1+mass_line+1].ha0;


	// 1,2
	if (targ_Mi<0.92*(isochrones[feh_line1+age_line2+599].Mi - isochrones[feh_line1+age_line2].Mi)+isochrones[feh_line1+age_line2].Mi)
	{
		mass_line=floor(300*(targ_Mi - isochrones[feh_line1+age_line2].Mi)/(0.92*(isochrones[feh_line1+age_line2+599].Mi - isochrones[feh_line1+age_line2].Mi)));
		mass_weight=300*(targ_Mi - isochrones[feh_line1+age_line2].Mi)/(0.92*(isochrones[feh_line1+age_line2+599].Mi - isochrones[feh_line1+age_line2].Mi)) - mass_line;
	}
	else
	{
		mass_line=300 + floor(300*(targ_Mi - isochrones[feh_line1+age_line2+300].Mi)/(0.08*(isochrones[feh_line1+age_line2+599].Mi - isochrones[feh_line1+age_line2].Mi)));
		mass_weight=300 + 300*(targ_Mi - isochrones[feh_line1+age_line2+300].Mi)/(0.08*(isochrones[feh_line1+age_line2+599].Mi - isochrones[feh_line1+age_line2].Mi)) - mass_line;
	}

	r1+=(1-feh_weight)*(age_weight)*(1-mass_weight) * isochrones[feh_line1+age_line2+mass_line].r0
	  + (1-feh_weight)*(age_weight)*(mass_weight) * isochrones[feh_line1+age_line2+mass_line+1].r0;
	i1+=(1-feh_weight)*(age_weight)*(1-mass_weight) * isochrones[feh_line1+age_line2+mass_line].i0
	  + (1-feh_weight)*(age_weight)*(mass_weight) * isochrones[feh_line1+age_line2+mass_line+1].i0;
	ha1+=(1-feh_weight)*(age_weight)*(1-mass_weight) * isochrones[feh_line1+age_line2+mass_line].ha0
	  + (1-feh_weight)*(age_weight)*(mass_weight) * isochrones[feh_line1+age_line2+mass_line+1].ha0;

	// 2,1
	if (targ_Mi<0.92*(isochrones[feh_line2+age_line1+599].Mi - isochrones[feh_line2+age_line1].Mi)+isochrones[feh_line2+age_line1].Mi)
	{
		mass_line=floor(300*(targ_Mi - isochrones[feh_line2+age_line1].Mi)/(0.92*(isochrones[feh_line2+age_line1+599].Mi - isochrones[feh_line2+age_line1].Mi)));
		mass_weight=300*(targ_Mi - isochrones[feh_line2+age_line1].Mi)/(0.92*(isochrones[feh_line2+age_line1+599].Mi - isochrones[feh_line2+age_line1].Mi)) - mass_line;
	}
	else
	{
		mass_line=300 + floor(300*(targ_Mi - isochrones[feh_line2+age_line1+300].Mi)/(0.08*(isochrones[feh_line2+age_line1+599].Mi - isochrones[feh_line2+age_line1].Mi)));
		mass_weight=300 + 300*(targ_Mi - isochrones[feh_line2+age_line1+300].Mi)/(0.08*(isochrones[feh_line2+age_line1+599].Mi - isochrones[feh_line2+age_line1].Mi)) - mass_line;
	}

	r1+=(feh_weight)*(1-age_weight)*(1-mass_weight) * isochrones[feh_line2+age_line1+mass_line].r0
	  + (feh_weight)*(1-age_weight)*(mass_weight) * isochrones[feh_line2+age_line1+mass_line+1].r0;
	i1+=(feh_weight)*(1-age_weight)*(1-mass_weight) * isochrones[feh_line2+age_line1+mass_line].i0
	  + (feh_weight)*(1-age_weight)*(mass_weight) * isochrones[feh_line2+age_line1+mass_line+1].i0;
	ha1+=(feh_weight)*(1-age_weight)*(1-mass_weight) * isochrones[feh_line2+age_line1+mass_line].ha0
	  + (feh_weight)*(1-age_weight)*(mass_weight) * isochrones[feh_line2+age_line1+mass_line+1].ha0;


	// 2,2
	if (targ_Mi<0.92*(isochrones[feh_line2+age_line2+599].Mi - isochrones[feh_line2+age_line2].Mi)+isochrones[feh_line2+age_line2].Mi)
	{
		mass_line=floor(300*(targ_Mi - isochrones[feh_line2+age_line2].Mi)/(0.92*(isochrones[feh_line2+age_line2+599].Mi - isochrones[feh_line2+age_line2].Mi)));
		mass_weight=300*(targ_Mi - isochrones[feh_line2+age_line2].Mi)/(0.92*(isochrones[feh_line2+age_line2+599].Mi - isochrones[feh_line2+age_line2].Mi)) - mass_line;
	}
	else
	{
		mass_line=300 + floor(300*(targ_Mi - isochrones[feh_line2+age_line2+300].Mi)/(0.08*(isochrones[feh_line2+age_line2+599].Mi - isochrones[feh_line2+age_line2].Mi)));
		mass_weight=300 + 300*(targ_Mi - isochrones[feh_line2+age_line2+300].Mi)/(0.08*(isochrones[feh_line2+age_line2+599].Mi - isochrones[feh_line2+age_line2].Mi)) - mass_line;
	}

	r1+=(feh_weight)*(age_weight)*(1-mass_weight) * isochrones[feh_line2+age_line2+mass_line].r0
	  + (feh_weight)*(age_weight)*(mass_weight) * isochrones[feh_line2+age_line2+mass_line+1].r0;
	i1+=(feh_weight)*(age_weight)*(1-mass_weight) * isochrones[feh_line2+age_line2+mass_line].i0
	  + (feh_weight)*(age_weight)*(mass_weight) * isochrones[feh_line2+age_line2+mass_line+1].i0;
	ha1+=(feh_weight)*(age_weight)*(1-mass_weight) * isochrones[feh_line2+age_line2+mass_line].ha0
	  + (feh_weight)*(age_weight)*(mass_weight) * isochrones[feh_line2+age_line2+mass_line+1].ha0;

	if (r1<-10||i1<-10||ha1<-10){throw(7);}
	iso_obj new_iso(targ_feh, targ_Mi, targ_logAge, 0., 0., r1, i1, ha1,1);

//*/
	//if ((r1>15||i1>15||ha1>15)){cout <<targ_feh << " " << targ_Mi << " " << targ_logAge << " "<<feh_weight << " " <<(targ_feh+0.914)/0.136 << " " << age_weight  << " " << r1 << " " << i1 << " " << ha1 << endl;}
//	iso_obj new_iso(0., 1., 9., 3.76, 4.5, 4.5, 4.15, 4.28);	
	return new_iso;
}


double max_age(double targ_Mi, vector<iso_obj> &isochrones)
{
	double max_age1=99.;

	if (targ_Mi<0.15 || targ_Mi>19.14)
	{
		throw (5);
	}


	for (double age=10.25; age>=7.05; age-=0.05)
	{
		if (isochrones[312000+ceil((10.25-age)/0.05)*600+599].Mi>targ_Mi)
		{
			max_age1=age;
			break;
		}
	}


	if (max_age1<50){return max_age1;}
	else {throw (6);return max_age1;}
}


iso_obj iso_get_Tg(double targ_feh, double targ_logT, double targ_logg, vector<iso_obj> &isochrones)
{
	int feh_line1, feh_line2;
	int logT_line1, logT_line2;
	int logg_line;
	double logg_weight, feh_weight, logT_weight;
	double r1=0, i1=0, ha1=0, Jac1=0, Mi1=0, logAge1=0;

	if (targ_feh<-0.914 || targ_feh>0.870 || targ_logT<3.5 || targ_logT>4.495 || targ_logg<-0.5 || targ_logg>5.495)
	{
		throw (5);
	}


	if (targ_feh>=-0.914 && targ_feh<-0.778){feh_line1=0;}
	else if (targ_feh>=-0.778 && targ_feh<-0.646){feh_line1=60701;}
	else if (targ_feh>=-0.646 && targ_feh<-0.493){feh_line1=121402;}
	else if (targ_feh>=-0.493 && targ_feh<-0.340){feh_line1=182103;}
	else if (targ_feh>=-0.340 && targ_feh<-0.226){feh_line1=242804;}
	else if (targ_feh>=-0.226 && targ_feh<-0.134){feh_line1=303505;}
	else if (targ_feh>=-0.134 && targ_feh<-0.058){feh_line1=364206;}
	else if (targ_feh>=-0.058 && targ_feh<0.022){feh_line1=424907;}
	else if (targ_feh>=0.022 && targ_feh<0.104){feh_line1=485608;}
	else if (targ_feh>=0.104 && targ_feh<0.205){feh_line1=546309;}
	else if (targ_feh>=0.205 && targ_feh<0.348){feh_line1=607010;}
	else if (targ_feh>=0.348 && targ_feh<0.455){feh_line1=667711;}
	else if (targ_feh>=0.455 && targ_feh<0.512){feh_line1=728412;}
	else if (targ_feh>=0.512 && targ_feh<0.570){feh_line1=789113;}
	else if (targ_feh>=0.570 && targ_feh<0.688){feh_line1=849814;}
	else if (targ_feh>=0.688 && targ_feh<0.870){feh_line1=910515;}
	else {throw(5);}

	logT_line1=int((targ_logT-3.5)/0.01)*601;
	logg_line=int((targ_logg+0.5)/0.01);

	Mi1+=isochrones[feh_line1+logT_line1+logg_line].Mi;
	logAge1+=isochrones[feh_line1+logT_line1+logg_line].logAge;
	r1+=isochrones[feh_line1+logT_line1+logg_line].r0;
	i1+=isochrones[feh_line1+logT_line1+logg_line].i0;
	ha1+=isochrones[feh_line1+logT_line1+logg_line].ha0;
	Jac1+=isochrones[feh_line1+logT_line1+logg_line].Jac;


	if (r1<-10 || i1<-10 || ha1<-10 || Jac1==0. || Mi1<0.){throw(7);}

	iso_obj new_iso(targ_feh, Mi1, logAge1, targ_logT, targ_logg, r1, i1, ha1, Jac1);

	return new_iso;
}


// gets stdout following the use of a command - used to read in Schlegel Galactic reddening - code 'borrowed' from somewhere on the net
string getStdoutFromCommand(string cmd)				
{
	// setup
	string str;
	FILE *stream;
	char buffer[100];

	// do it
	stream = popen(cmd.c_str(), "r");
	while ( fgets(buffer, 100, stream) != NULL )
		str.append(buffer);
	pclose(stream);

	// exit
	return str;
}



vector<bin_obj2> backup_A_mean_find(double l_gal, double b_gal)
{
	double Sch_max, density_dust, A_6250;
	vector<bin_obj2> backup_A_mean (150);

	// retrieve Schlegel et al limit
	Sch_max=SFD_read(l_gal, b_gal)*3.1/1.215;
	cout << "Sch_max = " << Sch_max << endl;

	// integrate dust density to ~infinity, used to normalise the dust distribution so that at infinity it gives the Schlegel value
	double dust_inf=0;
	for (double d=0; d<50000; d+=10)
	{
		dust_inf+=exp(-sqrt(pow(8080.,2)+pow(d*cos(b_gal*PI/180.),2)-2.*8080.*d*cos(b_gal*PI/180.)*cos(l_gal*PI/180.))/2500 - fabs(d*sin(b_gal*PI/180.)+17)/125)*10;	// Dust scale height and lengh from Marshall et al 2006
	}

	double const_term=Sch_max/dust_inf;

	A_6250=0;
	for (double d=0.0; d<=15001.0; d+=10.0)
	{
		density_dust=exp(-sqrt(pow(8080.,2)+pow(d*cos(b_gal*PI/180.),2)-2.*8080.*d*cos(b_gal*PI/180.)*cos(l_gal*PI/180.))/2500 - fabs(d*sin(b_gal*PI/180.)+17)/125);
		A_6250+=const_term*10*density_dust;		// max/total_int * delta_d * rho(d)

		if (d/100!=int(d/100) && d/50==int(d/50))				
		{
         		//cout << "d/100=" << d/100 << " A=" << A_6250 << " " << backup_A_mean.size() << endl; 
										// also make backup_A_mean at this point
			backup_A_mean[int((d-50)/100)].mean_A=A_6250;
			backup_A_mean[int((d-50)/100)].sigma=0.1*A_6250;
			backup_A_mean[int((d-50)/100)].d_mean=0.1;
		}
	}
	return backup_A_mean;   
}

string stringify(double x)		
{  
	ostringstream o;
	if(!(o << x)){throw 66;}
	return o.str();
}

double StrToDbl(string s) 
{
     double d;
     stringstream ss(s); //turn the string into a stream
     ss >> d; //convert
     return d;
}

vector <vector <string> > config_read(string filename)
{
	vector <vector <string> > totalfile;

	ifstream input1;
	input1.open(filename.c_str());
	if(!input1) { //output file couldn't be opened
		cerr << "Error: file " << filename << " could not be opened" << endl;
		exit(1);
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

			string buf;
			stringstream ss(str);		// turn that line into a stringstream
		
			vector<string> fromfile;	//vector to put contents of line into
		
			while (ss>>buf)
			{			// Includes implicit conversion from string to string
				fromfile.push_back(buf);	
			}
			if (fromfile.size()==7)		// check there's something in the line
			{
				totalfile.push_back(fromfile);
			}
		}
	}

	return totalfile;
}

//vector <vector <double> > pack_A_mean(double input[])
//{
//	vector < vector <double> > output(150, vector<double> (4));
//	for (int it=0; it<600; it++)
//	{
//		output[int(floor(it/4))][it%2]=input[it];
//	}

//	return output;
//}

//void vector_send(vector < vector <double> > input, int tag, int dest)
//{
//	double input1[input.size()*input[0].size()];
//	for (int it1=0; it1<input.size(); it1++)
//	{
//		for (int it2=0; it2<input[it1].size(); it2++)
//		{
//			input1[it1*4+it2]=input[it1][it2];
//		}
//	}

//	MPI_Send(&input1, 600, MPI_DOUBLE, tag, dest, MPI_COMM_WORLD);
//}

//vector <vector <double> > vector_recv(int tag, int source)
//{
//	MPI_Status Stat;
//	vector <vector <double> > output;
//	double output1[600];

//	MPI_Recv(&output1, 600, MPI_DOUBLE, tag, source, MPI_COMM_WORLD, &Stat);

//	output=pack_A_mean(output1);

//	return output;
//}

