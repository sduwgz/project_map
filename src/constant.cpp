#include "constant.h"

//log4cxx config file
const char* kLogConfig = "log4cxx.properties";
const char* defaultOutFile = "map_position_data";

//Length of integral interval
const int Interval = 100;

//Number of enzyme cleavage site
const int MAX_MISS_MATCH = 5;
const int MIN_MATCH_NUMBER = 5;

//Initial value of scoring matrix
const double INIT_SCORE = -100000.0;

//Parameter in deletion model
const int UNIT_LENGTH = 10000;
const int MAX_DELETION = 20;
