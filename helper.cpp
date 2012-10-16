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

	if (targ_feh<-0.914 || targ_feh>0.417 || targ_logT<3.5 || targ_logT>4.495 || targ_logg<-0.5 || targ_logg>5.495)
	{
		throw (5);
	}


	if (targ_feh>=-0.914 && targ_feh<-0.778){feh_line1=0; feh_line2=240000;feh_weight=(targ_feh+0.914)/0.136;}
	else if (targ_feh>=-0.778 && targ_feh<-0.652){feh_line1=240000; feh_line2=480000;feh_weight=(targ_feh+0.778)/0.126;}
	else if (targ_feh>=-0.652 && targ_feh<-0.472){feh_line1=480000; feh_line2=720000;feh_weight=(targ_feh+0.652)/0.180;}
	else if (targ_feh>=-0.472 && targ_feh<-0.343){feh_line1=720000; feh_line2=960000;feh_weight=(targ_feh+0.472)/0.129;}
	else if (targ_feh>=-0.343 && targ_feh<-0.243){feh_line1=1200000; feh_line2=1440000;feh_weight=(targ_feh+0.343)/0.100;}
	else if (targ_feh>=-0.243 && targ_feh<-0.160){feh_line1=1440000; feh_line2=1680000;feh_weight=(targ_feh+0.243)/0.083;}
	else if (targ_feh>=-0.160 && targ_feh<-0.090){feh_line1=1680000; feh_line2=1920000;feh_weight=(targ_feh+0.160)/0.070;}
	else if (targ_feh>=-0.090 && targ_feh<0.000){feh_line1=1920000; feh_line2=2160000;feh_weight=(targ_feh+0.090)/0.090;}
	else if (targ_feh>=0.000 && targ_feh<0.077){feh_line1=2160000; feh_line2=2400000;feh_weight=(targ_feh)/0.077;}
	else if (targ_feh>=0.077 && targ_feh<0.202){feh_line1=2400000; feh_line2=2640000;feh_weight=(targ_feh-0.077)/0.125;}
	else if (targ_feh>=0.202 && targ_feh<0.363){feh_line1=2640000; feh_line2=2880000;feh_weight=(targ_feh-0.202)/0.161;}
	else if (targ_feh>=0.363 && targ_feh<0.417){feh_line1=288000; feh_line2=3120000;feh_weight=(targ_feh-0.363)/0.054;}
	else {throw(5);}


	logT_line1=floor((targ_logT-3.5)/0.005)*1200;
	logT_line2=logT_line1+1200;
	logT_weight=(targ_logT-3.5)/0.005-logT_line1/1200;

	logg_line=floor((targ_logg+0.5)/0.005);
	logg_weight=(targ_logg+0.5)/0.005 - logg_line;

	// 1,1

	Mi1+=(1-feh_weight)*(1-logT_weight)*(1-logg_weight) * isochrones[feh_line1+logT_line1+logg_line].Mi
	  + (1-feh_weight)*(1-logT_weight)*(logg_weight) * isochrones[feh_line1+logT_line1+logg_line+1].Mi;
	logAge1+=(1-feh_weight)*(1-logT_weight)*(1-logg_weight) * isochrones[feh_line1+logT_line1+logg_line].logAge
	  + (1-feh_weight)*(1-logT_weight)*(logg_weight) * isochrones[feh_line1+logT_line1+logg_line+1].logAge;
	r1+=(1-feh_weight)*(1-logT_weight)*(1-logg_weight) * isochrones[feh_line1+logT_line1+logg_line].r0
	  + (1-feh_weight)*(1-logT_weight)*(logg_weight) * isochrones[feh_line1+logT_line1+logg_line+1].r0;
	i1+=(1-feh_weight)*(1-logT_weight)*(1-logg_weight) * isochrones[feh_line1+logT_line1+logg_line].i0
	  + (1-feh_weight)*(1-logT_weight)*(logg_weight) * isochrones[feh_line1+logT_line1+logg_line+1].i0;
	ha1+=(1-feh_weight)*(1-logT_weight)*(1-logg_weight) * isochrones[feh_line1+logT_line1+logg_line].ha0
	  + (1-feh_weight)*(1-logT_weight)*(logg_weight) * isochrones[feh_line1+logT_line1+logg_line+1].ha0;
	Jac1+=(1-feh_weight)*(1-logT_weight)*(1-logg_weight) * isochrones[feh_line1+logT_line1+logg_line].Jac
	  + (1-feh_weight)*(1-logT_weight)*(logg_weight) * isochrones[feh_line1+logT_line1+logg_line+1].Jac;

	if (r1<-10||i1<-10||ha1<-10||isochrones[feh_line1+logT_line1+logg_line].Jac==0.||isochrones[feh_line1+logT_line1+logg_line+1].Jac==0.){throw(7);}

	// 1,2

	Mi1+=(1-feh_weight)*(logT_weight)*(1-logg_weight) * isochrones[feh_line1+logT_line2+logg_line].Mi
	  + (1-feh_weight)*(logT_weight)*(logg_weight) * isochrones[feh_line1+logT_line2+logg_line+1].Mi;
	logAge1+=(1-feh_weight)*(logT_weight)*(1-logg_weight) * isochrones[feh_line1+logT_line2+logg_line].logAge
	  + (1-feh_weight)*(logT_weight)*(logg_weight) * isochrones[feh_line1+logT_line2+logg_line+1].logAge;
	r1+=(1-feh_weight)*(logT_weight)*(1-logg_weight) * isochrones[feh_line1+logT_line2+logg_line].r0
	  + (1-feh_weight)*(logT_weight)*(logg_weight) * isochrones[feh_line1+logT_line2+logg_line+1].r0;
	i1+=(1-feh_weight)*(logT_weight)*(1-logg_weight) * isochrones[feh_line1+logT_line2+logg_line].i0
	  + (1-feh_weight)*(logT_weight)*(logg_weight) * isochrones[feh_line1+logT_line2+logg_line+1].i0;
	ha1+=(1-feh_weight)*(logT_weight)*(1-logg_weight) * isochrones[feh_line1+logT_line2+logg_line].ha0
	  + (1-feh_weight)*(logT_weight)*(logg_weight) * isochrones[feh_line1+logT_line2+logg_line+1].ha0;
	Jac1+=(1-feh_weight)*(logT_weight)*(1-logg_weight) * isochrones[feh_line1+logT_line2+logg_line].Jac
	  + (1-feh_weight)*(logT_weight)*(logg_weight) * isochrones[feh_line1+logT_line2+logg_line+1].Jac;

	if (r1<-10||i1<-10||ha1<-10||isochrones[feh_line1+logT_line2+logg_line].Jac==0.||isochrones[feh_line1+logT_line2+logg_line+1].Jac==0.){throw(7);}

	// 2,1

	Mi1+=(feh_weight)*(1-logT_weight)*(1-logg_weight) * isochrones[feh_line2+logT_line1+logg_line].Mi
	  + (feh_weight)*(1-logT_weight)*(logg_weight) * isochrones[feh_line2+logT_line1+logg_line+1].Mi;
	logAge1+=(feh_weight)*(1-logT_weight)*(1-logg_weight) * isochrones[feh_line2+logT_line1+logg_line].logAge
	  + (feh_weight)*(1-logT_weight)*(logg_weight) * isochrones[feh_line2+logT_line1+logg_line+1].logAge;
	r1+=(feh_weight)*(1-logT_weight)*(1-logg_weight) * isochrones[feh_line2+logT_line1+logg_line].r0
	  + (feh_weight)*(1-logT_weight)*(logg_weight) * isochrones[feh_line2+logT_line1+logg_line+1].r0;
	i1+=(feh_weight)*(1-logT_weight)*(1-logg_weight) * isochrones[feh_line2+logT_line1+logg_line].i0
	  + (feh_weight)*(1-logT_weight)*(logg_weight) * isochrones[feh_line2+logT_line1+logg_line+1].i0;
	ha1+=(feh_weight)*(1-logT_weight)*(1-logg_weight) * isochrones[feh_line2+logT_line1+logg_line].ha0
	  + (feh_weight)*(1-logT_weight)*(logg_weight) * isochrones[feh_line2+logT_line1+logg_line+1].ha0;
	Jac1+=(feh_weight)*(1-logT_weight)*(1-logg_weight) * isochrones[feh_line2+logT_line1+logg_line].Jac
	  + (feh_weight)*(1-logT_weight)*(logg_weight) * isochrones[feh_line2+logT_line1+logg_line+1].Jac;

	if (r1<-10||i1<-10||ha1<-10||isochrones[feh_line2+logT_line1+logg_line].Jac==0.||isochrones[feh_line2+logT_line1+logg_line+1].Jac==0.){throw(7);}

	// 2,2

	Mi1+=(feh_weight)*(logT_weight)*(1-logg_weight) * isochrones[feh_line2+logT_line2+logg_line].Mi
	  + (feh_weight)*(logT_weight)*(logg_weight) * isochrones[feh_line2+logT_line2+logg_line+1].Mi;
	logAge1+=(feh_weight)*(logT_weight)*(1-logg_weight) * isochrones[feh_line2+logT_line2+logg_line].logAge
	  + (feh_weight)*(logT_weight)*(logg_weight) * isochrones[feh_line2+logT_line2+logg_line+1].logAge;
	r1+=(feh_weight)*(logT_weight)*(1-logg_weight) * isochrones[feh_line2+logT_line2+logg_line].r0
	  + (feh_weight)*(logT_weight)*(logg_weight) * isochrones[feh_line2+logT_line2+logg_line+1].r0;
	i1+=(feh_weight)*(logT_weight)*(1-logg_weight) * isochrones[feh_line2+logT_line2+logg_line].i0
	  + (feh_weight)*(logT_weight)*(logg_weight) * isochrones[feh_line2+logT_line2+logg_line+1].i0;
	ha1+=(feh_weight)*(logT_weight)*(1-logg_weight) * isochrones[feh_line2+logT_line2+logg_line].ha0
	  + (feh_weight)*(logT_weight)*(logg_weight) * isochrones[feh_line2+logT_line2+logg_line+1].ha0;
	Jac1+=(feh_weight)*(logT_weight)*(1-logg_weight) * isochrones[feh_line2+logT_line2+logg_line].Jac
	  + (feh_weight)*(logT_weight)*(logg_weight) * isochrones[feh_line2+logT_line2+logg_line+1].Jac;


	if (r1<-10||i1<-10||ha1<-10||isochrones[feh_line2+logT_line2+logg_line].Jac==0.||isochrones[feh_line2+logT_line2+logg_line+1].Jac==0.){throw(7);}
	iso_obj new_iso(targ_feh, Mi1, logAge1, targ_logT, targ_logg, r1, i1, ha1, Jac1);

//*/
	//if ((r1>15||i1>15||ha1>15)){cout <<targ_feh << " " << targ_Mi << " " << targ_logAge << " "<<feh_weight << " " <<(targ_feh+0.914)/0.136 << " " << age_weight  << " " << r1 << " " << i1 << " " << ha1 << endl;}
//	iso_obj new_iso(0., 1., 9., 3.76, 4.5, 4.5, 4.15, 4.28);	
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



vector<float> backup_A_mean_find(double l_gal, double b_gal, float s_R, float s_z, bool Sch_get)
{
	float Sch_max, density_dust;
	vector<float> backup_A_mean (150);

	// retrieve Schlegel et al limit
	if (Sch_get)
	{
		string Sch_string="./CodeC/lambert_getval CodeC/SFD_dust_4096_ngp.fits CodeC/SFD_dust_4096_sgp.fits 1 "; 	
		Sch_string.append(stringify(l_gal));
		Sch_string.append(" ");
		Sch_string.append(stringify(b_gal));
		//Sch_max=atof(getStdoutFromCommand(Sch_string).c_str())*2.944;		// 2.944 to convert E(B-V) given by Schlegel to A_6250
		Sch_max=4.0;
	}
	else {Sch_max=1.;}

	// integrate dust density to ~infinity, used to normalise the dust distribution so that at infinity it gives the Schlegel value
	density_dust=0;
	for (double d=0; d<30001.0; d+=10)
	{	// Dust scale height and lengh from Marshall et al 2006
		//dust_inf+=exp(-sqrt(pow(8080.,2)+pow(d*cos(b_gal*PI/180.),2)-2.*8080.*d*cos(b_gal*PI/180.)*cos(l_gal*PI/180.))/2500 - fabs(d*sin(b_gal*PI/180.)+17)/125)*10;
		density_dust+=exp(-sqrt(pow(8080.,2)+pow(d*cos(b_gal*PI/180.),2)-2.*8080.*d*cos(b_gal*PI/180.)*cos(l_gal*PI/180.))/s_R - fabs(d*sin(b_gal*PI/180.)+17)/s_z)*10.;
		if (d/100!=int(d/100) && d/50==int(d/50) && d<15000)				
		{
			backup_A_mean[int((d-50)/100)]=density_dust;
		}
	}

	float const_term=Sch_max/backup_A_mean[149];

	for (int it=0.0; it<150; it++)
	{
		backup_A_mean[it]*=const_term;
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
			if (fromfile.size()==4)		// check there's something in the line
			{
				totalfile.push_back(fromfile);
			}
		}
	}

	return totalfile;
}

vector <vector <double> > pack_A_mean(double input[])
{
	vector < vector <double> > output(150, vector<double> (4));
	for (int it=0; it<600; it++)
	{
		output[int(floor(it/4))][it%2]=input[it];
	}

	return output;
}

void vector_send(vector < vector <double> > input, int tag, int dest)
{
	double input1[input.size()*input[0].size()];
	for (int it1=0; it1<input.size(); it1++)
	{
		for (int it2=0; it2<input[it1].size(); it2++)
		{
			input1[it1*4+it2]=input[it1][it2];
		}
	}

	MPI_Send(&input1, 600, MPI_DOUBLE, tag, dest, MPI_COMM_WORLD);
}

vector <vector <double> > vector_recv(int tag, int source)
{
	MPI_Status Stat;
	vector <vector <double> > output;
	double output1[600];

	MPI_Recv(&output1, 600, MPI_DOUBLE, tag, source, MPI_COMM_WORLD, &Stat);

	output=pack_A_mean(output1);

	return output;
}

