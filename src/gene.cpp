#include "gene.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp> 
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("mole.main"));

bool Gene::getDistance() {
    _distance.resize(_position.size() - 1);
    for (int i = 0; i < _position.size() - 1; ++ i) {
        _distance[i] = _position[i + 1] - _position[i];
    }
    return true;
}

bool GeneReader::read(Gene& gene) {
    if(!_stream) {
        return false;
    }
    std::string buf;
    std::vector< std::string > data;
    while (std::getline(_stream, buf)) {
        boost::algorithm::trim(buf);
        if (buf.empty()) continue;
        if(boost::algorithm::starts_with(buf, "#")) continue;
        boost::algorithm::split(data, buf, boost::algorithm::is_any_of(" \t"), boost::algorithm::token_compress_on);
        LOG4CXX_INFO(logger, boost::format("bnx=>the ref id is: %s") % data[0]);
        gene._position.push_back(static_cast< int > (boost::lexical_cast< double > (data[5])));
    }
    gene.getDistance();
    return true;
}
