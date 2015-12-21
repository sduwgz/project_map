// Last Update:2015-04-27 10:35:36
/**
 * @file SplitString.h
 * @brief 
 * @author Qing Xu, xuqinguestc@163.com
 * @version 0.1.00
 * @date 2015-04-20
 */

#ifndef _SPLIT_STRING_H
#define _SPLIT_STRING_H

#include "Types.h"
using namespace std;

class SplitString
{
private:
    string str;
public:
    SplitString(const string &_s){str = _s;}
    vector<int> split2Int(const char *, int rep=0);
    vector<string> split2Str(const char *, int rep=0);
    vector<double> split2Dbl(const char *, int rep=0);
};



#endif  /*_SPLIT_STRING_H*/
