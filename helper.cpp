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
	iso_obj new_iso(targ_feh, targ_Mi, targ_logAge, 0., 0., r1, i1, ha1);

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


iso_obj iso_get-Tg(double targ_feh, double targ_logT, double targ_logg, vector<iso_obj> &isochrones)
{
	int feh_line1, feh_line2;
	int logT_line1, logT_line2;
	int logg_line;
	double logg_weight, feh_weight, logT_weight;
	double r1=0, i1=0, ha1=0;

	if (targ_feh<-0.914 || targ_feh>0.417 || targ_logT< || targ_logT> || targ_logg< || targ_logg>)
	{
		throw (5);
	}


	if (targ_feh>=-0.914 && targ_feh<-0.778){feh_line1=0; feh_line2=000;feh_weight=(targ_feh+0.914)/0.136;}
	else if (targ_feh>=-0.778 && targ_feh<-0.652){feh_line1=000; feh_line2=000;feh_weight=(targ_feh+0.778)/0.126;}
	else if (targ_feh>=-0.652 && targ_feh<-0.472){feh_line1=000; feh_line2=000;feh_weight=(targ_feh+0.652)/0.180;}
	else if (targ_feh>=-0.472 && targ_feh<-0.343){feh_line1=000; feh_line2=000;feh_weight=(targ_feh+0.472)/0.129;}
	else if (targ_feh>=-0.343 && targ_feh<-0.243){feh_line1=000; feh_line2=000;feh_weight=(targ_feh+0.343)/0.100;}
	else if (targ_feh>=-0.243 && targ_feh<-0.160){feh_line1=000; feh_line2=000;feh_weight=(targ_feh+0.243)/0.083;}
	else if (targ_feh>=-0.160 && targ_feh<-0.090){feh_line1=000; feh_line2=000;feh_weight=(targ_feh+0.160)/0.070;}
	else if (targ_feh>=-0.090 && targ_feh<0.000){feh_line1=000; feh_line2=000;feh_weight=(targ_feh+0.090)/0.090;}
	else if (targ_feh>=0.000 && targ_feh<0.077){feh_line1=000; feh_line2=000;feh_weight=(targ_feh)/0.077;}
	else if (targ_feh>=0.077 && targ_feh<0.202){feh_line1=000; feh_line2=000;feh_weight=(targ_feh-0.077)/0.125;}
	else if (targ_feh>=0.202 && targ_feh<0.363){feh_line1=000; feh_line2=000;feh_weight=(targ_feh-0.202)/0.161;}
	else if (targ_feh>=0.363 && targ_feh<0.417){feh_line1=000; feh_line2=000;feh_weight=(targ_feh-0.363)/0.054;}
	else {throw(5);}


	logT_line1=floor((targ_logT-3.5)/0.01)*;
	logT_line2=logT_line1+1000;
	logT_weight=(logT_line2-targ_logT)/0.01;

	if (isochrones[feh_line1+logT_line1].logg>targ_logg || isochrones[feh_line1+logT_line1+1000].logg<targ_logg || isochrones[feh_line1+logT_line2].logg>targ_logg || isochrones[feh_line1+logT_line2+1000].logg<targ_logg ||
	   isochrones[feh_line2+logT_line1].logg>targ_logg || isochrones[feh_line2+logT_line1+1000].logg<targ_logg || isochrones[feh_line2+logT_line2].logg>targ_logg || isochrones[feh_line2+logT_line2+1000].logg<targ_logg)
	{
		throw (6);
	}

	// 1,1
	if (targ_logg>3.4)
	{
		logg_line=floor(700*(isochrones[feh_line1+logT_line1].logg - targ_logg)/(isochrones[feh_line1+logT_line1].logg -3.5));
		logg_weight=700*(isochrones[feh_line1+logT_line1].logg - targ_logg)/(isochrones[feh_line1+logT_line1].logg -3.5) - logg_line;
	}
	else
	{
		logg_line=700 + floor(300*(isochrones[feh_line1+logT_line1].logg - targ_logg)/(3.5-isochrones[feh_line1+logT_line1+999].logg));
		logg_weight=700 + 300*(isochrones[feh_line1+logT_line1].logg - targ_logg)/(3.5-isochrones[feh_line1+logT_line1+999].logg) - logg_line;
	}

	r1+=(1-feh_weight)*(1-logT_weight)*(1-logg_weight) * isochrones[feh_line1+logT_line1+logg_line].r0
	  + (1-feh_weight)*(1-logT_weight)*(logg_weight) * isochrones[feh_line1+logT_line1+logg_line+1].r0;
	i1+=(1-feh_weight)*(1-logT_weight)*(1-logg_weight) * isochrones[feh_line1+logT_line1+logg_line].i0
	  + (1-feh_weight)*(1-logT_weight)*(logg_weight) * isochrones[feh_line1+logT_line1+logg_line+1].i0;
	ha1+=(1-feh_weight)*(1-logT_weight)*(1-logg_weight) * isochrones[feh_line1+logT_line1+logg_line].ha0
	  + (1-feh_weight)*(1-logT_weight)*(logg_weight) * isochrones[feh_line1+logT_line1+logg_line+1].ha0;

	// 1,2
	if (targ_logg>3.4)
	{
		logg_line=floor(700*(isochrones[feh_line1+logT_line2].logg - targ_logg)/(isochrones[feh_line1+logT_line2].logg -3.5));
		logg_weight=700*(isochrones[feh_line1+logT_line2].logg - targ_logg)/(isochrones[feh_line1+logT_line2].logg -3.5) - logg_line;
	}
	else
	{
		logg_line=700 + floor(300*(isochrones[feh_line1+logT_line2].logg - targ_logg)/(3.5-isochrones[feh_line1+logT_line2+999].logg));
		logg_weight=700 + 300*(isochrones[feh_line1+logT_line2].logg - targ_logg)/(3.5-isochrones[feh_line1+logT_line2+999].logg) - logg_line;
	}

	r1+=(1-feh_weight)*(logT_weight)*(1-logg_weight) * isochrones[feh_line1+logT_line2+logg_line].r0
	  + (1-feh_weight)*(logT_weight)*(logg_weight) * isochrones[feh_line1+logT_line2+logg_line+1].r0;
	i1+=(1-feh_weight)*(logT_weight)*(1-logg_weight) * isochrones[feh_line1+logT_line2+logg_line].i0
	  + (1-feh_weight)*(logT_weight)*(logg_weight) * isochrones[feh_line1+logT_line2+logg_line+1].i0;
	ha1+=(1-feh_weight)*(logT_weight)*(1-logg_weight) * isochrones[feh_line1+logT_line2+logg_line].ha0
	  + (1-feh_weight)*(logT_weight)*(logg_weight) * isochrones[feh_line1+logT_line2+logg_line+1].ha0;

	// 2,1
	if (targ_logg>3.4)
	{
		logg_line=floor(700*(isochrones[feh_line2+logT_line1].logg - targ_logg)/(isochrones[feh_line2+logT_line1].logg -3.5));
		logg_weight=700*(isochrones[feh_line2+logT_line1].logg - targ_logg)/(isochrones[feh_line2+logT_line1].logg -3.5) - logg_line;
	}
	else
	{
		logg_line=700 + floor(300*(isochrones[feh_line2+logT_line1].logg - targ_logg)/(3.5-isochrones[feh_line2+logT_line1+999].logg));
		logg_weight=700 + 300*(isochrones[feh_line2+logT_line1].logg - targ_logg)/(3.5-isochrones[feh_line2+logT_line1+999].logg) - logg_line;
	}

	r1+=(feh_weight)*(1-logT_weight)*(1-logg_weight) * isochrones[feh_line2+logT_line1+logg_line].r0
	  + (feh_weight)*(1-logT_weight)*(logg_weight) * isochrones[feh_line2+logT_line1+logg_line+1].r0;
	i1+=(feh_weight)*(1-logT_weight)*(1-logg_weight) * isochrones[feh_line2+logT_line1+logg_line].i0
	  + (feh_weight)*(1-logT_weight)*(logg_weight) * isochrones[feh_line2+logT_line1+logg_line+1].i0;
	ha1+=(feh_weight)*(1-logT_weight)*(1-logg_weight) * isochrones[feh_line2+logT_line1+logg_line].ha0
	  + (feh_weight)*(1-logT_weight)*(logg_weight) * isochrones[feh_line2+logT_line1+logg_line+1].ha0;

	// 2,2
	if (targ_logg>3.4)
	{
		logg_line=floor(700*(isochrones[feh_line2+logT_line2].logg - targ_logg)/(isochrones[feh_line2+logT_line2].logg -3.5));
		logg_weight=700*(isochrones[feh_line2+logT_line2].logg - targ_logg)/(isochrones[feh_line2+logT_line2].logg -3.5) - logg_line;
	}
	else
	{
		logg_line=700 + floor(300*(isochrones[feh_line2+logT_line2].logg - targ_logg)/(3.5-isochrones[feh_line2+logT_line2+999].logg));
		logg_weight=700 + 300*(isochrones[feh_line2+logT_line2].logg - targ_logg)/(3.5-isochrones[feh_line2+logT_line2+999].logg) - logg_line;
	}

	r1+=(feh_weight)*(logT_weight)*(1-logg_weight) * isochrones[feh_line2+logT_line2+logg_line].r0
	  + (feh_weight)*(logT_weight)*(logg_weight) * isochrones[feh_line2+logT_line2+logg_line+1].r0;
	i1+=(feh_weight)*(logT_weight)*(1-logg_weight) * isochrones[feh_line2+logT_line2+logg_line].i0
	  + (feh_weight)*(logT_weight)*(logg_weight) * isochrones[feh_line2+logT_line2+logg_line+1].i0;
	ha1+=(feh_weight)*(logT_weight)*(1-logg_weight) * isochrones[feh_line2+logT_line2+logg_line].ha0
	  + (feh_weight)*(logT_weight)*(logg_weight) * isochrones[feh_line2+logT_line2+logg_line+1].ha0;


	if (r1<-10||i1<-10||ha1<-10){throw(7);}
	iso_obj new_iso(targ_feh, targ_logT, targ_logg, 0., 0., r1, i1, ha1);

//*/
	//if ((r1>15||i1>15||ha1>15)){cout <<targ_feh << " " << targ_Mi << " " << targ_logAge << " "<<feh_weight << " " <<(targ_feh+0.914)/0.136 << " " << age_weight  << " " << r1 << " " << i1 << " " << ha1 << endl;}
//	iso_obj new_iso(0., 1., 9., 3.76, 4.5, 4.5, 4.15, 4.28);	
	return new_iso;
}


