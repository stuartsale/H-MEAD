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

double in_sample(double r_in, double i_in, double ha_in, double dr_in, double di_in, double dha_in)
{
	double r_s, i_s, ha_s;
	double r_ha_max, r_i_min;
	double sum=0;
	double A_int;

	for (int n=0; n<50; n++)
	{
		r_s=r_in+gsl_ran_gaussian_ziggurat(rng_handle,dr_in);
		i_s=i_in+gsl_ran_gaussian_ziggurat(rng_handle,di_in);
		ha_s=ha_in+gsl_ran_gaussian_ziggurat(rng_handle,dha_in);
		
		A_int=(+(r_s-i_s)-0.53738)/(0.233289);
		r_ha_max=-0.00248*pow(A_int,2) + 0.06848*A_int + 0.005775 + 0.32275;
		if (r_s-ha_s<r_ha_max)
		{
			sum++;
		}
	}
	if (sum>0) {return sum/50;}
	else {return 0.0002;}		
}

double in_sample2(double A_max, double mu, double sigma)
{
	double A_s;
	double sum=0;
	double A_m2=log(A_max);
	for (int n=0; n<50; n++)
	{
	//	A_s=gsl_ran_lognormal(rng_handle,mu, sigma);	
		A_s=mu+gsl_ran_gaussian_ziggurat(rng_handle,sigma);
		if (A_s<A_m2){sum++;}
	}
	if (sum>0) {return sum/50;}
	else {return 0.0002;}		
}
