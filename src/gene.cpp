#include "gene.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp> 
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("mole.main"));

bool Gene::getDistance() {
    distance.resize(position.size() - 1); 
    for (int i = 0; i < position.size() - 1; ++ i) {
        distance[i] = position[i + 1] - position[i];
    }
    return true;      
}

bool GeneReader::read(Gene& gene) {
    if(!_stream) {
        return false;
    }
    gene.position.clear();
    gene.distance.clear();
    std::string buf;
    std::vector< std::string > data;
    while (std::getline(_stream, buf)) {
        boost::algorithm::trim(buf);
        if (buf.empty()) continue;
        if(boost::algorithm::starts_with(buf, "#")) continue;
        boost::algorithm::split(data, buf, boost::algorithm::is_any_of(" "), boost::algorithm::token_compress_on);
        if (boost::lexical_cast< int > (data[0]) == 1) {
            gene.position.push_back(static_cast< int > (boost::lexical_cast< double > (data[5])));
        } else {
            LOG4CXX_WARN(logger, boost::format("bnx=>invalid line for gene file: %s") % buf);
            return false;
        }
    }
    gene.getDistance();
    return true;
}
