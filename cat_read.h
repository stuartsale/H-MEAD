#include "iphas_obj.h"

using namespace std;

#ifndef PI
#define PI 3.14159265358979323846
#endif


vector<iphas_obj> iphas_read(string filename,float &r_min1,float &i_min1,float &ha_min1,float &r_max1, float &i_max1, float &ha_max1);
vector<iphas_obj> TWOMASS_read(string filename,float &J_min1,float &H_min1,float &K_min1,float &J_max1, float &H_max1, float &K_max1);
