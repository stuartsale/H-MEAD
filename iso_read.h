#ifndef ISOGET_H_
#define ISOGET_H_

#include <vector>
#include <sstream>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

#include "iso_obj.h"

vector<iso_obj> iso_read(const string &filename);
vector<iso_obj> iso_read_Tg(const string &filename);
vector<iso_obj> iso_read_Tg_2MASS(const string &filename);

#endif
