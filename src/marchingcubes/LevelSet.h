#ifndef _LEVEL_SET_H_
#define _LEVEL_SET_H_

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include "MarchingCubes.h"

using namespace std;

int sign(double num);

int input_levelset_grid_info(std::string filename, MarchingCubes& mc);

int calculate_normal(MarchingCubes& mc);

#endif
