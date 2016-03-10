#ifndef constant_h_
#define constant_h_

//log4cxx config file
extern const char* kLogConfig;
extern const char* defaultOutFile;
extern const char* defaultParameterFile;

//Length of integral interval
extern const int Interval;

//Number of enzyme cleavage site
extern const int MAX_MISS_MATCH;
extern const int MIN_MATCH_NUMBER;

//Initial value of scoring matrix
extern const double INIT_SCORE;

//Parameter in deletion model
extern const int UNIT_LENGTH;
extern const int MAX_DELETION;

#endif
