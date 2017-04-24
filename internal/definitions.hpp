/*
 * definitions.hpp
 *
 *  Created on: Apr 13, 2017
 *      Author: nico
 */

#ifndef INCLUDE_DEFINITIONS_HPP_
#define INCLUDE_DEFINITIONS_HPP_

#include <string>
#include <stdio.h>
#include <vector>
#include <cstring>
#include <cmath>
#include <climits>
#include <cstdlib>
#include <iostream>
#include <istream>
#include <sys/stat.h>
#include <set>
#include "stdint.h"
#include <sstream>
#include <algorithm>
#include <sys/resource.h>
#include <sys/stat.h>
#include <fstream>
#include <assert.h>
#include <tuple>

using namespace std;

namespace ri{

typedef uint64_t ulint;
typedef long int lint;
typedef unsigned int uint;
typedef unsigned short int t_char; ///< Type for char conversion
typedef unsigned short int t_errors; ///< Type for ERRORS

typedef unsigned char uchar;
typedef unsigned char symbol;
typedef unsigned char uint8;

typedef pair<ulint,ulint> range_t;

}

#endif /* INCLUDE_DEFINITIONS_HPP_ */
