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
#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>

typedef std::pair< int, int > BackTrace;

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("map.main"));

bool Map::multiRun(MoleSet& moleSet, const Gene& gene, int threadNumber) const {
    std::cout << moleSet[0]._id << std::endl;
    int block = moleSet.size() / threadNumber + 1;
    MoleSet *moleSetPtr = &moleSet;
    boost::thread_group group;
    for (int i = 0; i != threadNumber; ++ i){
        group.create_thread(boost::bind(&Map::start, this, moleSetPtr, gene, i, block));
    }
    group.join_all();
}

bool Map::start(MoleSet* moleSetPtr, const Gene& gene, int i, int block) const {
    MoleSet& moleSet = * moleSetPtr;
    int start = i * block;
    if (start > moleSet.size()) {
        return true;
    }
    int end = start + block;
    if (end > moleSet.size()) {
        end = moleSet.size();
    }
    for (int i = start; i != end; ++ i) {
        if (moleSet[i]._distance.size() < MIN_MATCH_NUMBER) {
            LOG4CXX_DEBUG(logger, boost::format("mole %s is too short.") % (moleSet[i]._id));
            continue;
        };
        wholeDPscore(moleSet[i], gene._distance);
    }
    return true;
}

bool Map::run(MoleSet& moleSet, const Gene& gene) const {
    for (int i = 0; i < moleSet.size(); i += 2) {
        if (moleSet[i]._distance.size() < MIN_MATCH_NUMBER) {
            LOG4CXX_DEBUG(logger, boost::format("mole %s is too short.") % (moleSet[i]._id));
            continue;
        };
        wholeDPscore(moleSet[i], gene._distance);
        wholeDPscore(moleSet[i + 1], gene._distance);
        if (moleSet[i].mapRet.score <= INIT_SCORE && moleSet[i + 1].mapRet.score <= INIT_SCORE) {
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
    int delta_laplace = abs(moleLength - geneLength - 46);
    int delta_background = abs(moleLength - geneLength - 1870);

    int moleSiteNumber = moleFragment.size() - 1;
    int geneSiteNumber = geneFragment.size() - 1;

    if (delta_laplace > 100000) {
        delta_laplace = 100000;
    }
    if (delta_background > 100000) {
        delta_background = 100000;
    }
    int deleteNumber = static_cast< int > ((geneSiteNumber + 0.0) / moleLength * UNIT_LENGTH + 0.5);
    if (deleteNumber < 1) {
        deleteNumber = 1;
    } else if (deleteNumber > MAX_DELETION) {
        deleteNumber = MAX_DELETION;
    }
    //LOG4CXX_DEBUG(logger, boost::format("_laplaceScore: %s, _backgroundScore: %s, _insertionScore: %s, _deletionScore: %s") % _laplaceScore[delta] %_backgroundScore[delta] %_insertionScore[moleSiteNumber] %_deletionScore[deleteNumber]);
    //LOG4CXX_DEBUG(logger, boost::format("probLaplace: %s, probBackground: %s, probInsertion: %s, probDeletion: %s") % probLaplace(delta) %probBackground(delta) %probInsertion(moleSiteNumber) %probDeletion(geneSiteNumber, geneLength));
    

    if (moleSiteNumber != 0 || geneSiteNumber != 0) {
        //return _laplaceScore[delta_laplace] + _deletionScore[deleteNumber] + _insertionScore[moleSiteNumber] - _backgroundScore[delta_background];
        return probLaplace(delta) + probDeletion(geneSiteNumber, moleLength) + probInsertion(moleSiteNumber) - probBackground(delta);
    } else {
        //return _laplaceScore[delta_laplace]  - _backgroundScore[delta_background];
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
    boost::math::normal_distribution<> n(0, sigma);
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
    //boost::math::laplace_distribution<> l(mu, sigma);
    boost::math::laplace_distribution<> l(0, sigma);
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
    double insertPunish = probInsertion(1);
    for (int i = 0; i < MAX_MISS_MATCH; ++ i) {
        for (int j = 0; j < cols; ++ j) {
            scoreMatrix[i][j] = insertPunish * i;
        }
    }
    for (int i = 0; i < rows; ++ i) {
        scoreMatrix[i][0] = insertPunish * i;
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
    for (int i = 0; i < MAX_MISS_MATCH; ++ i) {
        for (int j = 0; j < cols; ++ j) {
            if (scoreMatrix[mole._distance.size() - i][j] + insertPunish * i > max) {
                max = scoreMatrix[mole._distance.size() - i][j] + insertPunish * i;
                mole.mapRet.alignEndPosition.first = mole._distance.size() - i;
                mole.mapRet.alignEndPosition.second = j;
            }
        }
    }

    int pi = mole.mapRet.alignEndPosition.first, pj = mole.mapRet.alignEndPosition.second;
    while(true) {
        int pii = pi, pjj = pj;
        pi = backTrace[pii][pjj].first, pj = backTrace[pii][pjj].second;
        if (pi == -1 && pj == -1) {
            mole.mapRet.alignStartPosition.first = pii;
            mole.mapRet.alignStartPosition.second = pjj;
            break;
        }

        int moleFragmentLength = accumulate(mole._distance.begin() + pi, mole._distance.begin() + pii, 0);
        int geneFragmentLength = accumulate(gene.begin() + pj, gene.begin() + pjj, 0);

        mole.mapRet.alignFragmentLength.push_back(std::make_pair(moleFragmentLength, geneFragmentLength));

        mole.mapRet.moleMapPosition.push_back(std::make_pair(pi, pii - 1));
        mole.mapRet.geneMapPosition.push_back(std::make_pair(pj, pjj - 1));
    }

    //reversed result is more human-facing
    reverse(mole.mapRet.moleMapPosition.begin(), mole.mapRet.moleMapPosition.end());
    reverse(mole.mapRet.geneMapPosition.begin(), mole.mapRet.geneMapPosition.end());
    reverse(mole.mapRet.alignFragmentLength.begin(), mole.mapRet.alignFragmentLength.end());

    mole.mapRet.score = max;
    return true;
}

void Map::output(const std::string& filename, const std::vector< Mole >& moleSet) const {
    std::ofstream out;
    out.open(filename.c_str());

    for (int i = 0; i<moleSet.size(); ++ i) {
        out << moleSet[i]._id << "\t" << "\t" << moleSet[i].mapRet.score << "\t"
            << moleSet[i].mapRet.alignStartPosition.first << "\t" << moleSet[i].mapRet.alignEndPosition.first << "\t" << moleSet[i].mapRet.alignStartPosition.second << "\t" << moleSet[i].mapRet.alignEndPosition.second << "\n";
        for (int j = 0; j<moleSet[i].mapRet.moleMapPosition.size(); ++ j) {
            out << moleSet[i].mapRet.moleMapPosition[j].first << "\t" << moleSet[i].mapRet.moleMapPosition[j].second << "\t" << moleSet[i].mapRet.geneMapPosition[j].first << "\t" << moleSet[i].mapRet.geneMapPosition[j].second << "\t" << moleSet[i].mapRet.alignFragmentLength[j].first << "\t" << moleSet[i].mapRet.alignFragmentLength[j].second << "\n";
        }
    }
}
void Map::printScore(const MoleSet& moleSet) const {
    for (int i = 0; i < moleSet.size(); ++ i) {
        Mole mole = moleSet[i];
        std::cout << mole._id << "\t" << mole.mapRet.score << "\t" << mole.mapRet.alignStartPosition.first + 1 << "\t" <<  mole.mapRet.alignEndPosition.first + 1 << "\t" << mole.mapRet.alignStartPosition.second << "\t" <<  mole.mapRet.alignEndPosition.second << std::endl;
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
    for (boost::property_tree::ptree::const_iterator it = parameters.begin(); it != parameters.end(); it++) {
        _parameters[it->first] = boost::lexical_cast<double> (it->second.data());
    }
    return true;
}
void Map::initPunishScore() {
    for (int i = 0; i < MAX_MISS_MATCH + 1; ++ i) {
        _insertionScore.push_back(probInsertion(i));
    }
    for (int i = 0; i < MAX_DELETION + 1; ++ i) {
        _deletionScore.push_back(probDeletion(i, UNIT_LENGTH));
    }
    for (int i = 0; i <= 100000; ++ i) {
        _laplaceScore.push_back(probLaplace(i));
        _backgroundScore.push_back(probBackground(i));
    }
}
