/**
 * @file SplitString.cc
 * @brief 
 * @author Qing Xu, xuqinguestc@163.com
 * @version 0.1.00
 * @date 2015-04-20
 */

#include "SplitString.h"

vector<double> SplitString::split2Dbl(const char *delim, int rep)
{
    vector<double> ret;
    string buf;
    set<char> exist;
    while(*delim != '\0')
    { 
        exist.insert(*delim);
        ++delim;
    }

    for(int i = 0 ;i < str.size() ; ++i)
    {
        if (exist.find(str[i]) != exist.end())
        {
            if (!buf.empty())
            {
                ret.push_back(atof(buf.c_str()));
                buf.clear();
            }
        }
        else
        {
            buf += str[i];
        }
    }
    if (!buf.empty())
        ret.push_back(atof(buf.c_str()));
    return ret;
}


vector<int> SplitString::split2Int(const char *delim, int rep)
{
    vector<int> ret;
    string buf;
    set<char> exist;
    while(*delim != '\0')
    { 
        exist.insert(*delim);
        ++delim;
    }

    for(int i = 0 ;i < str.size() ; ++i)
    {
        if (exist.find(str[i]) != exist.end())
        {
            if (!buf.empty())
            {
                ret.push_back(atoi(buf.c_str()));
                buf.clear();
            }
        }
        else
        {
            buf += str[i];
        }
    }
    if (!buf.empty())
        ret.push_back(atoi(buf.c_str()));
    return ret;
}


vector<string> SplitString::split2Str(const char *delim, int rep) 
{
    vector<string> ret;
    string buf;
    set<char> exist;

    while(*delim != '\0')
    { 
        exist.insert(*delim);
        ++delim;
    }

    for(int i = 0 ;i < str.size() ; ++i)
    {
        if (exist.find(str[i]) != exist.end())
        {
            if (!buf.empty())
            {
                ret.push_back(buf);
                buf.clear();
            }
        }
        else
        {
            buf += str[i];
        }
    }

    if (!buf.empty())
    {
        ret.push_back(buf);
    }
    return ret;
}
