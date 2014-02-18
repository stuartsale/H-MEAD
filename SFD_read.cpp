#include "SFD_read.h"


float SFD_read(float l, float b)
{
        std::valarray<float>  contents;

	if (b>0)
	{
	// Find x & Y coords
	float x,y;
	x = sqrt(1-sin(PI*b/180.) ) * cos(PI*l/180.) * 2048 + 2048.5;
	y = -sqrt(1-sin(PI*b/180.) ) * sin(PI*l/180.) * 2048 + 2048.5;
//		cout << x << " " << y << endl;
//	x=2*sin((90-b)*PI/360) * sin((l-90)*PI/180) /(-0.0395646818624*PI/180) + 2048.5;
//	y=2*sin((90-b)*PI/360) * cos((l-90)*PI/180) /(+0.0395646818624*PI/180) + 2048.5;

//		cout << x << " " << y << endl;

	if (x<1){x=1.000001;}
	if (x>4096){x=4095.999;}
	if (y<1){y=1.000001;}
	if (y>4096){y=4095.999;}

	vector <long> start (2,0);
	start[0]=floor(x);
	start[1]=floor(y);
	vector <long> end (2,0);
	end[0]=start[0]+1;
	end[1]=start[1]+1;
	vector <long> stride (2, 1);

	// Open file
		auto_ptr<CCfits::FITS> pInfile(new CCfits::FITS(config_dir+"SFD_dust_4096_ngp.fits",CCfits::Read,true));
        	CCfits::PHDU& image = pInfile->pHDU(); 
	// Read in required section of file
		image.read(contents, start, end, stride);

	// Interpolate to get required E(B-V)
		float sum=0.;
		float weights=0;
		for (int i=0; i<contents.size(); i++)
		{
			if (contents[i]>0)
			{
				weights+=(1-abs(x-(floor(x)+(i-2*floor(i/2)) ) ) )  * (1-abs(y-(floor(y)+floor(i/2) ) ) );
				sum+=contents[i] * (1-abs(x-(floor(x)+(i-2*floor(i/2)) ) ) )  * (1-abs(y-(floor(y)+floor(i/2) ) ) );
//				cout << i << " " << contents[i] << " " << (1-abs(x-(floor(x)+(i-2*floor(i/2)) ) ) )  * (1-abs(y-(floor(y)+floor(i/2) ) ) ) << " " << x << " " << y << endl;
			}
		}
//		cout << sum << " " << weights << endl;

		return sum/weights;
	}




	else if (b<0)

	{
	// Find x & Y coords
	float x,y;
//	x=2*sin((90-b)*PI/360) * sin((l-90)*PI/180) /(-0.0395646818624*PI/180) + 2048.5;
//	y=-2*sin((90-b)*PI/360) * cos((l-90)*PI/180) /(+0.0395646818624*PI/180) + 2048.5;
	x = sqrt(1-sin(PI*(-b)/180.) ) * cos(PI*l/180.) * 2048 + 2048.5;
	y = -sqrt(1-sin(PI*(-b)/180.) ) * sin(-PI*l/180.) * 2048 + 2048.5;

//		cout << x << " " << y << endl;

	if (x<1){x=1.000001;}
	if (x>4096){x=4095.999;}
	if (y<1){y=1.000001;}
	if (y>4096){y=4095.999;}

	vector <long> start (2,0);
	start[0]=floor(x);
	start[1]=floor(y);
	vector <long> end (2,0);
	end[0]=start[0]+1;
	end[1]=start[1]+1;
	vector <long> stride (2, 1);

	// Open file
		auto_ptr<CCfits::FITS> pInfile(new CCfits::FITS(config_dir+"SFD_dust_4096_sgp.fits",CCfits::Read,true));
        	CCfits::PHDU& image = pInfile->pHDU(); 
	// Read in required section of file
		image.read(contents, start, end, stride);

	// Interpolate to get required E(B-V)
		float sum=0.;
		float weights=0;
		for (int i=0; i<contents.size(); i++)
		{
			if (contents[i]>0)
			{
				weights+=(1-abs(x-(floor(x)+(i-2*floor(i/2)) ) ) )  * (1-abs(y-(floor(y)+floor(i/2) ) ) );
				sum+=contents[i] * (1-abs(x-(floor(x)+(i-2*floor(i/2)) ) ) )  * (1-abs(y-(floor(y)+floor(i/2) ) ) );
//				cout << i << " " << contents[i] << " " << (1-abs(x-(floor(x)+(i-2*floor(i/2)) ) ) )  * (1-abs(y-(floor(y)+floor(i/2) ) ) ) << " " << x << " " << y << endl;
			}
		}
//		cout << sum << " " << weights << endl;

		return sum/weights;
	}

	else if (b==0){return SFD_read(l, 0.00000001);}
		



}
