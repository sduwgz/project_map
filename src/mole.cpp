#include "mole.h"

#include <stdlib.h>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <cassert>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp> 
#include <boost/foreach.hpp>

//#include <log4cxx/logger.h>

Mole Mole::reverseMole() {
    Mole reMole;
    reMole.id = -id;
    reMole.length = length;
    reMole.enzymeNumber = enzymeNumber;
    std::vector< int > data;
    data.assign(distance.begin(), distance.end());
    reverse(data.begin(), data.end());
    reMole.distance.assign(data.begin(), data.end());
    return reMole;
}

bool Mole::getDistance() {
    distance.resize(position.size() - 2); 
    for (int i = 0; i < position.size() - 2; ++ i) {
        distance[i] = position[i + 1] - position[i];
    }
    return true;
}

bool MoleReader::read(Mole& mole) {
    if(!_stream){
        return false;
    }
    enum {
        moleId,
        molePosition,
        moleQX01,
        moleQX02,
    };

    int state = moleId;
    std::string buf;
    std::vector<std::string> data;
    mole.reset();
    while (std::getline(_stream, buf)) {
        boost::algorithm::trim(buf);
        if (buf.empty()) continue;
        if (state == moleId) {
            boost::algorithm::split(data, buf, boost::algorithm::is_any_of("\t"), boost::algorithm::token_compress_on);
            if (boost::lexical_cast<int>(data[0]) == 0) {
                mole.id = boost::lexical_cast<int> (data[1]);
                state = molePosition;
            } else {
                //LOG4CXX_WARN(logger, boost::format("bnx=>invalid line for mole id: %s") % buf);
                return false;
            }
        } else if (state == molePosition) {
            data = boost::algorithm::split(data, buf, boost::algorithm::is_any_of("\t"), boost::algorithm::token_compress_on);
            if (boost::lexical_cast<int>(data[0]) == 1) {
                for (int i = 1; i < data.size(); ++ i) {
                    mole.position.push_back(static_cast< long > (boost::lexical_cast< double >(data[i])));
                }
                if (mole.position.size() > 0) {
                    mole.getDistance();
                }
                state = moleQX01;
            } else {
                //LOG4CXX_WARN(logger, boost::format("bnx=>invalid line for mole position: %s") % buf);
                return false;
            }
        } else if (state == moleQX01) {
            state = moleQX02;
        } else if (state == moleQX02) {
            state = moleId;
            return true;
        }
    }
}
