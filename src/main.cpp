#include <iostream>
#include <cstdarg>
#include <unistd.h>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/lexical_cast.hpp>

#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>

#include "mole.h"
#include "gene.h"
#include "map.h"

typedef boost::property_tree::ptree Properties;

typedef std::vector< Mole > MoleSet;

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("map.main"));

static const char *opt_string = "m:g:o:c:p:t:r:";


int printHelps() {
    std::cout << "USAGE : ./map -g [input gene file] -m [input mole file] -o [output file] -c [log.config file]" << std::endl;
    return 1;
}

int main(int argc, char* argv[]) {
    if (argc < 5) {
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
    std::string geneFile = options.get< std::string >("g", "");
    if (boost::filesystem::exists(geneFile)) {
        std::ifstream geneInstream(geneFile.c_str());
        GeneReader gReader(geneInstream);
        if (!gReader.read(g)) {
            LOG4CXX_WARN(logger, boost::format("load %s failed.") % geneFile);
            return 1;
        } else {
            LOG4CXX_INFO(logger, boost::format("load %s successed.") % geneFile);
        }
    } else {
        LOG4CXX_WARN(logger, boost::format("%s is not existed.") % geneFile);
        return 1;
    }
    MoleSet moleSet;
    std::string moleFile = options.get< std::string >("m", "");
    std::string reverseLabel = options.get< std::string >("r", "1");
    if (boost::filesystem::exists(moleFile)) {
        std::ifstream moleInstream(moleFile.c_str());
        MoleReader mReader(moleInstream);
        Mole m;
        while (mReader.read(m)) {
            moleSet.push_back(m);
            if (reverseLabel == "1") {
                Mole reMole = m.reverseMole();
                moleSet.push_back(reMole);
            }
        }
        int moleNumber = moleSet.size();
        if (moleNumber == 0) {
            LOG4CXX_WARN(logger, "no mole is in moleSet");
            return 1;
        } else {
            LOG4CXX_INFO(logger, boost::format("%s moles have been inited.") % moleNumber);
        }
    } else {
        LOG4CXX_WARN(logger, boost::format("%s is not existed.") % moleFile);
        return 1;
    }
    std::string outFile = options.get< std::string >("o", defaultOutFile);
    std::string parameterFile = options.get< std::string >("p", defaultParameterFile);
    int threadNumber = boost::lexical_cast<int> (options.get< std::string >("t", "4"));
    
    Map maptool;
    if (!maptool.initParameters(parameterFile)) {
        LOG4CXX_WARN(logger, boost::format("%s, init parameter error.") % parameterFile);
        return 1;
    }
    maptool.initPunishScore();
    //maptool.run(moleSet, g);
    maptool.multiRun(moleSet, g, threadNumber);
    maptool.printScore(moleSet);
    maptool.output(outFile, moleSet);
    
    return 1;
}
