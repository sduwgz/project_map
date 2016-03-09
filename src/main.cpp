#include <iostream>
#include <cstdarg>
#include <unistd.h>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>

#include "mole.h"
#include "gene.h"
#include "map.h"

typedef boost::property_tree::ptree Properties;

typedef std::vector< Mole > MoleSet;

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("map.main"));

static const char *opt_string = "m:g:o:c:";


int printHelps() {
    std::cout << "USAGE : ./map -g [input gene file] -m [input mole file] -o [output file] -c [log.config file]" << std::endl;
    return 1;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        return printHelps();
    }

    Properties options;
    {
        int opt = -1;
        while ((opt = getopt(argc, argv, opt_string)) != -1) {
            std::string key = std::string(1, (char)opt);
            if (optarg == NULL) {
                options.put(key, NULL);
            } else {
                std::string val = optarg;
                options.put(key, val);
            }
        }
    }
    std::string log_config = options.get< std::string >("c", kLogConfig);
    if (boost::filesystem::exists(log_config)) {
        log4cxx::PropertyConfigurator::configure(log_config);
    } else {
        log4cxx::BasicConfigurator::configure();
    }

   
    Gene g;
    std::string gene_file = options.get< std::string >("g", "");
    if (boost::filesystem::exists(gene_file)) {
        std::ifstream geneInstream(gene_file.c_str());
        GeneReader gReader(geneInstream);
        if (!gReader.read(g)) {
            LOG4CXX_WARN(logger, boost::format("load %s failed.") % gene_file);
            return 1;
        } else {
            LOG4CXX_DEBUG(logger, boost::format("load %s successed.") % gene_file);
        }
    } else {
        LOG4CXX_WARN(logger, boost::format("%s is not existed.") % gene_file);
        return 1;
    }
    MoleSet moleSet;
    std::string mole_file = options.get< std::string >("m", "");
    if (boost::filesystem::exists(mole_file)) {
        std::ifstream moleInstream(mole_file.c_str());
        MoleReader mReader(moleInstream);
        Mole m;
        while (mReader.read(m)) {
            moleSet.push_back(m);
            Mole reMole = m.reverseMole();
            moleSet.push_back(reMole);
        }
        int mole_number = moleSet.size();
        if (mole_number == 0) {
            LOG4CXX_WARN(logger, "no mole is in moleSet");
            return 1;
        } else {
            LOG4CXX_DEBUG(logger, boost::format("%s moles have been inited.") % mole_number);
        }
    } else {
        LOG4CXX_WARN(logger, boost::format("%s is not existed.") % mole_file);
        return 1;
    }
    std::string out_file = options.get< std::string >("o", defaultOutFile);

    double mean=46.0, variance=570;
    double beta = 0.15;
    double alpha = 0.15;
 
    Map maptool(mean, variance, alpha, beta, out_file);
    maptool.whole_map_score(moleSet, g.distance);
    
    return 1;
}
