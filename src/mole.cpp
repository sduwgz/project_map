#include "mole.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp> 
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("mole.main"));

Mole Mole::reverseMole() {
    Mole reMole;
    reMole._id = "(-)" + _id;
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
    reset(mole);
    int state = moleId;
    std::string buf;
    std::vector<std::string> data;
    while (std::getline(_stream, buf)) {
        if(buf.size() == 0 || buf[0] == '#') continue;
        boost::algorithm::trim(buf);
        LOG4CXX_DEBUG(logger, boost::format("line: %s") % buf);

        if (buf.empty()) continue;
        if (state == moleId) {
            boost::algorithm::split(data, buf, boost::algorithm::is_any_of(" \t"), boost::algorithm::token_compress_on);
            if (boost::lexical_cast<int>(data[0]) == 0) {
                mole._id = data[1];
                state = molePosition;
            } else {
                LOG4CXX_WARN(logger, boost::format("bnx=>invalid line for mole id: %s") % mole._id);
                return false;
            }
        } else if (state == molePosition) {
            data = boost::algorithm::split(data, buf, boost::algorithm::is_any_of(" \t"), boost::algorithm::token_compress_on);
            if (boost::lexical_cast<int>(data[0]) == 1) {
                for (int i = 1; i < data.size() - 1; ++ i) {
                    mole._position.push_back(static_cast< long > (boost::lexical_cast< double >(data[i])));
                }
                state = moleQX01;
            } else {
                LOG4CXX_WARN(logger, boost::format("bnx=>invalid line for mole position: %s") % buf);
                return false;
            }
        } else if (state == moleQX01) {
            data = boost::algorithm::split(data, buf, boost::algorithm::is_any_of(" \t"), boost::algorithm::token_compress_on);
            if (data[0] == "QX11") {
                std::vector<double> qx11;
                for (int i = 1; i < data.size(); ++ i) {
                    qx11.push_back(static_cast< long > (boost::lexical_cast< double >(data[i])));
                }
                int l = 0;
                double theta = 3;
                //filter sites by QX11 score. There needs a parameter theta.
                if(qx11.size() == mole._position.size() - 1) {
                    for(int i = 0; i < qx11.size(); ++ i) {
                        if(qx11[i] >= theta) {
                            mole._position[l ++] = mole._position[i];
                        }
                    }
                    mole._position.resize(l);
                }
                if (mole._position.size() > 1) {
                    mole.getDistance();
                }
                state = moleQX02;
            } else {
                LOG4CXX_WARN(logger, boost::format("bnx=>invalid line for mole QX11: %s") % buf);
                return false;
            }

            state = moleQX02;
        } else if (state == moleQX02) {
            state = moleId;
            return true;
        }
    }
    return false;
}
