#include "mole.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp> 
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("mole.main"));

Mole Mole::reverseMole() {
    Mole reMole;
    reMole._id = -_id;
    reMole._distance.assign(_distance.begin(), _distance.end());
    reverse(reMole._distance.begin(), reMole._distance.end());
    return reMole;
}

bool Mole::getDistance() {
    _distance.clear();
    for (int i = 0; i < _position.size() - 2; ++ i) {
        _distance.push_back(_position[i + 1] - _position[i]);
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
    while (std::getline(_stream, buf)) {
        boost::algorithm::trim(buf);
        if (buf.empty()) continue;
        if (state == moleId) {
            boost::algorithm::split(data, buf, boost::algorithm::is_any_of("\t"), boost::algorithm::token_compress_on);
            if (boost::lexical_cast<int>(data[0]) == 0) {
                mole._id = boost::lexical_cast<int> (data[1]);
                state = molePosition;
            } else {
                LOG4CXX_WARN(logger, boost::format("bnx=>invalid line for mole id: %s") % mole._id);
                return false;
            }
        } else if (state == molePosition) {
            data = boost::algorithm::split(data, buf, boost::algorithm::is_any_of("\t"), boost::algorithm::token_compress_on);
            if (boost::lexical_cast<int>(data[0]) == 1) {
                for (int i = 1; i < data.size(); ++ i) {
                    mole._position.push_back(static_cast< long > (boost::lexical_cast< double >(data[i])));
                }
                if (mole._position.size() > 0) {
                    mole.getDistance();
                }
                state = moleQX01;
            } else {
                LOG4CXX_WARN(logger, boost::format("bnx=>invalid line for mole position: %s") % buf);
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
