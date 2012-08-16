#include "iphas_obj.h"

using namespace std;

#ifndef PI
#define PI 3.14159265358979323846
#endif


vector<iphas_obj> iphas_read(string filename,double &r_min1,double &i_min1,double &ha_min1,double &r_max1, double &i_max1, double &ha_max1);
vector<iphas_obj> TWOMASS_read(string filename,double &J_min1,double &H_min1,double &K_min1,double &J_max1, double &H_max1, double &K_max1);
