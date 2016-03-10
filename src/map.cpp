#include "map.h"

#include <math.h>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/laplace.hpp>
#include <boost/math/distributions/exponential.hpp>

#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>

typedef std::pair< int, int > BackTrace;

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("map.main"));

bool Map::run(MoleSet& moleSet, const std::vector <int>& gene) const {
    for (int i = 0; i < moleSet.size(); i += 2) {
        if (moleSet[i]._distance.size() < MIN_MATCH_NUMBER) {
            LOG4CXX_DEBUG(logger, boost::format("mole %s is too short.") % (moleSet[i]._id));
            continue;
        };
        wholeDPscore(moleSet[i], gene);
        wholeDPscore(moleSet[i + 1], gene);
        if (!moleSet[i].mapRet.label && !moleSet[i + 1].mapRet.label) {
            LOG4CXX_DEBUG(logger, boost::format("%s and %s can not map head to tail") % (moleSet[i]._id) % moleSet[i + 1]._id);
        }
    }
    return true;
}

double Map::validScore(const Fragment& moleFragment, const Fragment& geneFragment) const {
    int moleLength = 0, geneLength = 0;
    moleLength = std::accumulate(moleFragment.begin(), moleFragment.end(), 0);
    geneLength = std::accumulate(geneFragment.begin(), geneFragment.end(), 0);
    
    int delta = moleLength - geneLength;
   
    int moleSiteNumber = moleFragment.size() - 1;
    int geneSiteNumber = geneFragment.size() - 1;
    
    if (moleSiteNumber != 0 || geneSiteNumber != 0) {
        return probLaplace(delta) + probDeletion(geneSiteNumber, moleLength) + probInsertion(moleSiteNumber) - probBackground(delta);
    } else {
        return probLaplace(delta)  - probBackground(delta);
        //return laplace(delta) + pI(0) + pD(0, moleLength) - background(delta);
    }
}

double Map::probDeletion(int siteNumber, int moleLength) const {
    int deleteNumber = static_cast< int > ((siteNumber + 0.0) / moleLength * UNIT_LENGTH + 0.5);
    if (deleteNumber < 1) {
        deleteNumber = 1;
    } else if (deleteNumber > MAX_DELETION) {
        deleteNumber = MAX_DELETION;
    }
    double lambda = _parameters.find("lambda_poisson")->second;
    boost::math::poisson_distribution<> p(lambda);
    return log(boost::math::pdf(p, deleteNumber));
}

double Map::probInsertion(int k) const {
    double lambda = _parameters.find("lambda_exponent")->second;
    boost::math::exponential_distribution<> e(lambda);
    if (k == 0) {
        return log(cdf(e, 0.5) - cdf(e, 0));
    } 
    return log(cdf(e, k + 0.5) - cdf(e, k - 0.5));
}

double Map::probBackground(int delta) const {
    double mu = _parameters.find("mu_background")->second;
    double sigma = _parameters.find("sigma_background")->second;
    boost::math::normal_distribution<> n(mu, sigma);
    int distance = delta - mu;
    int d = distance / Interval;
    int interval_left = d * Interval;
    int interval_right = (d + 1) * Interval;
    if (distance < 0) {
        interval_left = (d - 1) * Interval;
        interval_right = d * Interval;
    }
    return log(boost::math::cdf(n, interval_right) - boost::math::cdf(n, interval_left));
}

double Map::probLaplace(int delta) const {
    double mu = _parameters.find("mu_laplace")->second;
    double sigma = _parameters.find("sigma_laplace")->second;
    boost::math::laplace_distribution<> l(mu, sigma);
    int distance = delta - mu;
    int d = distance / Interval;
    int interval_left = d * Interval;
    int interval_right = (d + 1) * Interval;
    if (distance < 0) {
        interval_left = (d - 1) * Interval;
        interval_right = d * Interval;
    }
    return log(boost::math::cdf(l, interval_right) - boost::math::cdf(l, interval_left));
}

bool Map::wholeDPscore(Mole& mole, const std::vector<int>& gene) const {
    int rows = mole._distance.size() + 1, cols = gene.size() + 1;
    double scoreMatrix[rows][cols];
    BackTrace  backTrace[rows][cols];
    for (int i = 0; i < rows; ++ i) {
        for (int j = 0; j < cols; ++ j) {
            scoreMatrix[i][j] = INIT_SCORE;
            backTrace[i][j].first = -1;
            backTrace[i][j].second = -1;
        }
    }
    for (int j = 0; j < cols; ++ j) {
        scoreMatrix[0][j] = 0.0;
        scoreMatrix[1][j] = probInsertion(1);
        scoreMatrix[2][j] = probInsertion(2);
    }
    for (int i = 0; i < rows; ++ i) {
        scoreMatrix[i][0] = -3 * i;
    }
    for (int i = 1; i < rows; ++ i) {
        for (int j = 1; j < cols; ++ j) {
            Fragment moleFragment, geneFragment;
            moleFragment.push_back(mole._distance[i - 1]);
            for (int k = j - 1; k >= j - MAX_MISS_MATCH && k >= 0; -- k) {
                geneFragment.push_back(gene[k]);
                double temp = scoreMatrix[i - 1][k] + validScore(moleFragment, geneFragment);
                if (temp > scoreMatrix[i][j]) {
                    scoreMatrix[i][j] = temp;
                    backTrace[i][j].first = i - 1;
                    backTrace[i][j].second = k;
                }
            }

            moleFragment.clear();
            geneFragment.clear();
            geneFragment.push_back(gene[j - 1]);
            for (int k = i - 1; k >= i - MAX_MISS_MATCH && k >= 0; -- k) {
                moleFragment.push_back(mole._distance[k]);
                double temp = scoreMatrix[k][j - 1] + validScore(moleFragment, geneFragment);
                if (temp > scoreMatrix[i][j]) {
                    scoreMatrix[i][j] = temp;
                    backTrace[i][j].first = k;
                    backTrace[i][j].second = j - 1;
                }
            }

            if (i > 2 && j > 2){

                moleFragment.clear();
                geneFragment.clear();
                geneFragment.push_back(gene[j - 1]);
                geneFragment.push_back(gene[j - 2]);
                moleFragment.push_back(mole._distance[i - 1]);
                moleFragment.push_back(mole._distance[i - 2]);
                
                double temp = scoreMatrix[i - 2][j - 2] + validScore(moleFragment, geneFragment);
                if (temp > scoreMatrix[i][j]) {
                    scoreMatrix[i][j] = temp;
                    backTrace[i][j].first = i - 2;
                    backTrace[i][j].second = j - 2;
                }
            } 
        }
    }

    double max = INIT_SCORE;
    for (int j = 0; j < cols; ++j) {
        if (scoreMatrix[mole._distance.size()][j] > max) {
            max = scoreMatrix[mole._distance.size()][j];
            mole.mapRet.alignMolePosition.second = mole._distance.size();
            mole.mapRet.alignGenePosition.second = j;
        }
        double insertPunish1 = probInsertion(1);
        if (scoreMatrix[mole._distance.size() - 1][j] + insertPunish1 > max) {
            max = scoreMatrix[mole._distance.size() - 1][j] + insertPunish1;
            mole.mapRet.alignMolePosition.second = mole._distance.size() - 1;
            mole.mapRet.alignGenePosition.second = j;
        }
        double insertPunish2 = probInsertion(2);
        if (scoreMatrix[mole._distance.size() - 2][j] + insertPunish2 > max) {
            max = scoreMatrix[mole._distance.size() - 2][j] + insertPunish2;
            mole.mapRet.alignMolePosition.second = mole._distance.size() - 2;
            mole.mapRet.alignGenePosition.second = j;
        }
    }
    int pi = mole.mapRet.alignMolePosition.second, pj = mole.mapRet.alignGenePosition.second;
    while(true) {
        int pii = pi, pjj = pj;
        pi = backTrace[pii][pjj].first, pj = backTrace[pii][pjj].second;
        if (pi == -1 && pj == -1) {
            //have gotten the head
            mole.mapRet.alignMolePosition.first = pii;
            mole.mapRet.alignGenePosition.first = pjj;
            break;
        }

        LenNum moleLn, geneLn;

        moleLn.num = pii - pi;
        geneLn.num = pjj - pj;

        moleLn.len = accumulate(mole._distance.begin() + pi, mole._distance.begin() + pii, 0);
        geneLn.len = accumulate(gene.begin() + pj, gene.begin() + pjj, 0);

        mole.mapRet.alignLenNum.push_back(std::make_pair(moleLn, geneLn));

        mole.mapRet.moleMapPosition.push_back(std::make_pair(pi, pii - 1));
        mole.mapRet.geneMapPosition.push_back(std::make_pair(pj, pjj - 1));
    }

    //reversed result is more human-facing
    reverse(mole.mapRet.moleMapPosition.begin(), mole.mapRet.moleMapPosition.end());
    reverse(mole.mapRet.geneMapPosition.begin(), mole.mapRet.geneMapPosition.end());
    reverse(mole.mapRet.alignLenNum.begin(), mole.mapRet.alignLenNum.end());

    std::cout << mole._id << "\t" << max << "\t" << mole.mapRet.alignGenePosition.first+1 << "\t" <<  mole.mapRet.alignGenePosition.second+1 << "\t" << mole.mapRet.alignMolePosition.first << "\t" <<  mole.mapRet.alignMolePosition.second << std::endl; 
    mole.mapRet.score = max;
    mole.mapRet.label = true; 
    return true;
}

void Map::print_score(const std::string& filename, const std::vector< Mole >& moleSet) const {
    std::ofstream out;
    out.open(filename.c_str());

    for (int i=0; i<moleSet.size(); i++) {
        out << moleSet[i]._id <<"\t" << moleSet[i].mapRet.label <<"\t" << moleSet[i].mapRet.score << "\t" 
            <<  moleSet[i].mapRet.alignMolePosition.first << "\t" << moleSet[i].mapRet.alignMolePosition.second << "\t"  <<  moleSet[i].mapRet.alignGenePosition.first << "\t" << moleSet[i].mapRet.alignGenePosition.second << "\n";
        for (int j=0; j<moleSet[i].mapRet.moleMapPosition.size(); j++) {
            out << moleSet[i].mapRet.moleMapPosition[j].first << "\t" << moleSet[i].mapRet.moleMapPosition[j].second << "\t" << moleSet[i].mapRet.geneMapPosition[j].first << "\t" << moleSet[i].mapRet.geneMapPosition[j].second << "\t" << moleSet[i].mapRet.alignLenNum[j].first.len << "\t" << moleSet[i].mapRet.alignLenNum[j].second.len << "\n";
        }
    }
}

bool Map::initParameters(const std::string& parameter_file) {
    boost::property_tree::ptree parameters;
    if (boost::filesystem::exists(parameter_file)) {
        try {
            boost::property_tree::read_ini(parameter_file, parameters);
        } catch (const boost::property_tree::ini_parser_error& e) {
            LOG4CXX_ERROR(logger, boost::format("load %s failed(%s).") % parameter_file % e.what());
            return false;
        }
    } else {
        LOG4CXX_WARN(logger, boost::format("%s is not existed.") % parameter_file);
        return false;
    }
    for (boost::property_tree::ptree::const_iterator it = parameters.begin(); it != parameters.end(); it++){
       _parameters[it->first] = boost::lexical_cast<double> (it->second.data());
    }
    return true;
}
