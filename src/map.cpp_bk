/**
 * @file map.cpp
 * @author Yanbo Li, liyanbo@ict.ac.cn
 * @author Guozheng Wei, weiguozheng@ict.ac.cn
 * @date 2015-12-30
 */

#include <math.h>
#include <numeric>

#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/laplace.hpp>
#include <boost/math/distributions/exponential.hpp>

#include "map.h"
#include "Types.h"
#include "mole.h"
#include "constant.h"


bool Map::whole_map_score(std::vector<Mole>& moleSet, const std::vector <int>& gene) const {
    for (int i = 0; i < moleSet.size(); i += 2) {
        moleSet[i].resetMapRet();
        moleSet[i + 1].resetMapRet();
        //MINCNT is a threshold of the min site in mole
        if (moleSet[i].distance.size() < _MINCNT) continue;
        //dp
        whole_DP_score(moleSet[i], gene);
        whole_DP_score(moleSet[i + 1], gene);
        //between the reversed mole and the mole, choose the better one
        if (!moleSet[i].mapRet.label && !moleSet[i + 1].mapRet.label) {
            std::cerr << moleSet[i].id << " and " << moleSet[i + 1].id << " can not map head to tail" << std::endl;
        }
    }
    std::string filename = _outPrefix;
    print_score(filename, moleSet);
    return true;
}

double Map::validScore(const Fragment& moleFragment, const Fragment& geneFragment) const {
    int moleLength = 0, geneLength = 0;
    moleLength = std::accumulate(moleFragment.begin(), moleFragment.end(), 0);
    geneLength = std::accumulate(geneFragment.begin(), geneFragment.end(), 0);
    /*
         ***********************************
         *                                 *
         *                   delta         *
         *                   |   |         *
         *  MOLE:  ---*----------*---      *
         *  GENE:  ---*------*-------      *
         *                                 *
         ***********************************
     */
    
    int delta = moleLength - geneLength;
    
    /*
         ****************************************************
         *                                                  *
         *         INSERT:                                  *
         *             MOLE:  ---*----*--*---               *
         *             GENE:  ---*-------*---               *
         *         DELETE:                                  *
         *             MOLE:  ---*-------*---               *
         *             GENE:  ---*----*--*---               *
         *         MATCH:                                   *
         *             MOLE:  ---*-------*---               *
         *             GENE:  ---*-------*---               *
         *                                                  *
         ****************************************************
     */
   
    int moleSiteNumber = moleFragment.size() - 1;
    int geneSiteNumber = geneFragment.size() - 1;
    
    if (moleSiteNumber != 0 || geneSiteNumber != 0) {
        return laplace(delta) + pD(geneSiteNumber, moleLength) + probInsertion(moleSiteNumber) - background(delta);
    } else {
        //Match
        return laplace(delta)  - background(delta);
        //return laplace(delta) + pI(0) + pD(0, moleLength) - background(delta);
    }
}

double Map::pD(int siteNumber, int moleLength) const {
    int deleteNumber = (int)((siteNumber + 0.0) / moleLength * UNIT_LENGTH + 0.5);
    if (deleteNumber < 1) {
        deleteNumber = 1;
    } else if (deleteNumber > MAX_DELETION) {
        deleteNumber = MAX_DELETION;
    }
    double lambd = 0.225;
    boost::math::poisson_distribution<> p(lambd);
    return log(boost::math::pdf(p, deleteNumber));
}

double Map::probInsertion(int k) const {
    boost::math::exponential_distribution<> e(5);
    if (k == 0) {
        return log(cdf(e, 0.5) - cdf(e, 0));
    } 
    return log(cdf(e, k + 0.5) - cdf(e, k - 0.5));
}

double Map::background(int delta) const {
    double mu = 1870.0, sigma = 10840.0;
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

double Map::laplace(int delta) const {
    double mu = 46.0;
    double sigma = 403.0;
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

bool Map::whole_DP_score(Mole& mole, const std::vector<int>& gene) const {
    int rows = mole.distance.size() + 1, cols = gene.size() + 1;
    double dp[rows][cols];
    std::pair<int,int> backTrack[rows][cols];
    for (int i = 0; i <= mole.distance.size(); ++ i) {
        for (int j = 0; j <= gene.size(); ++ j) {
            dp[i][j] = INIT_SCORE;
            backTrack[i][j].first = -1;
            backTrack[i][j].second = -1;
        }
    }
    //the first row is inited to zero, and the first and second rows are inited to pI(1) and pI(2)
    for (int j = 0; j <= gene.size(); ++ j) {
        dp[0][j] = 0.0;
        dp[1][j] = probInsertion(1);
        dp[2][j] = probInsertion(2);
    }
    for (int i = 0; i <= mole.distance.size(); ++ i) {
        dp[i][0] = -3 * i;
    }
    for (int i = 1; i <= mole.distance.size(); ++ i) {
        for (int j = 1; j <= gene.size(); ++ j) {
            Fragment moleFragment, geneFragment;
            moleFragment.push_back(mole.distance[i - 1]);
            for (int tj = j - 1; tj >= j - 5 && tj >= 0; -- tj) {
                geneFragment.push_back(gene[tj]);
                double temp = dp[i - 1][tj] + validScore(moleFragment, geneFragment);
                if (temp > dp[i][j]) {
                    dp[i][j] = temp;
                    backTrack[i][j].first = i - 1;
                    backTrack[i][j].second = tj;
                }
            }

            moleFragment.clear();
            geneFragment.clear();
            geneFragment.push_back(gene[j - 1]);
            for (int ti = i - 1; ti >= i - 5 && ti >= 0; -- ti) {
                moleFragment.push_back(mole.distance[ti]);
                double temp = dp[ti][j - 1] + validScore(moleFragment, geneFragment);
                if (temp > dp[i][j]) {
                    dp[i][j] = temp;
                    backTrack[i][j].first = ti;
                    backTrack[i][j].second = j - 1;
                }
            }

            if (i > 2 && j > 2){

                moleFragment.clear();
                geneFragment.clear();
                geneFragment.push_back(gene[j - 1]);
                geneFragment.push_back(gene[j - 2]);
                moleFragment.push_back(mole.distance[i - 1]);
                moleFragment.push_back(mole.distance[i - 2]);
                
                double temp = dp[i - 2][j - 2] + validScore(moleFragment, geneFragment);
                if (temp > dp[i][j]) {
                    dp[i][j] = temp;
                    backTrack[i][j].first = i - 2;
                    backTrack[i][j].second = j - 2;
                }
            } 
        }
    }

    double max = INIT_SCORE;

    for (int j = 0; j <= gene.size(); ++j) {
        //we should find the max score in the last row
        if (dp[mole.distance.size()][j] > max) {
            max = dp[mole.distance.size()][j];
            mole.mapRet.alignMolePosition.second = mole.distance.size();
            mole.mapRet.alignGenePosition.second = j;
        }
        //if the last one is a insertion, we must find the max score in last but one row, and give a punish
        double insertPunish1 = probInsertion(1);
        if (dp[mole.distance.size() - 1][j] + insertPunish1 > max) {
            max = dp[mole.distance.size() - 1][j] + insertPunish1;
            mole.mapRet.alignMolePosition.second = mole.distance.size() - 1;
            mole.mapRet.alignGenePosition.second = j;
        }
        //if the last two is insertions
        double insertPunish2 = probInsertion(2);
        if (dp[mole.distance.size() - 2][j] + insertPunish2 > max) {
            max = dp[mole.distance.size() - 2][j] + insertPunish2;
            mole.mapRet.alignMolePosition.second = mole.distance.size() - 2;
            mole.mapRet.alignGenePosition.second = j;
        }
    }
    
    int pi = mole.mapRet.alignMolePosition.second, pj = mole.mapRet.alignGenePosition.second;
    if (max == INIT_SCORE) {
        std::cerr << "case 1 map failure" << std::endl; 
        return false;
    }
    
    //trace back, and find the path
    while(true) {
        
        int pii = pi, pjj = pj;
        pi = backTrack[pii][pjj].first, pj = backTrack[pii][pjj].second;
        
        if (pi == -1 && pj == -1) {
            //have gotten the head
            mole.mapRet.alignMolePosition.first = pii;
            mole.mapRet.alignGenePosition.first = pjj;
            break;
        }

        //moleLn and genLn are the structs of storing the result
        LenNum moleLn, geneLn;

        assert(pi == pii - 1 || pj == pjj - 1 || pi == pii - 2 || pj == pjj - 2);

        moleLn.num = pii - pi;
        geneLn.num = pjj - pj;

        moleLn.len = accumulate(mole.distance.begin() + pi, mole.distance.begin() + pii, 0);
        geneLn.len = accumulate(gene.begin() + pj, gene.begin() + pjj, 0);

        mole.mapRet.alignLenNum.push_back(std::make_pair(moleLn, geneLn));

        mole.mapRet.moleMapPosition.push_back(std::make_pair(pi, pii - 1));
        mole.mapRet.geneMapPosition.push_back(std::make_pair(pj, pjj - 1));
    }

    //reversed result is more human-facing
    reverse(mole.mapRet.moleMapPosition.begin(), mole.mapRet.moleMapPosition.end());
    reverse(mole.mapRet.geneMapPosition.begin(), mole.mapRet.geneMapPosition.end());
    reverse(mole.mapRet.alignLenNum.begin(), mole.mapRet.alignLenNum.end());

    std::cout << mole.id << "\t" << max << "\t" << mole.mapRet.alignGenePosition.first+1 << "\t" <<  mole.mapRet.alignGenePosition.second+1 << "\t" << mole.mapRet.alignMolePosition.first << "\t" <<  mole.mapRet.alignMolePosition.second << std::endl; 
    
    //max the highest score in the last row
    mole.mapRet.score = max;
    /*
    if (mole.mapRet.alignMolePosition.first < 3) {
        mole.mapRet.label = true; 
        return true;
    }

    cerr << "case 2 map failure" <<endl; 
    return false;
    */
    mole.mapRet.label = true; 
    return true;
}

void Map::print_score(const std::string& filename, const std::vector< Mole >& moleSet) const {
    std::ofstream out;
    out.open(filename.c_str());

    for (int i=0; i<moleSet.size(); i++) {
        out << moleSet[i].id <<"\t" << moleSet[i].mapRet.label <<"\t" << moleSet[i].mapRet.score << "\t" 
            <<  moleSet[i].mapRet.alignMolePosition.first << "\t" << moleSet[i].mapRet.alignMolePosition.second << "\t"  <<  moleSet[i].mapRet.alignGenePosition.first << "\t" << moleSet[i].mapRet.alignGenePosition.second << "\n";
        for (int j=0; j<moleSet[i].mapRet.moleMapPosition.size(); j++) {
            out << moleSet[i].mapRet.moleMapPosition[j].first << "\t" << moleSet[i].mapRet.moleMapPosition[j].second << "\t" << moleSet[i].mapRet.geneMapPosition[j].first << "\t" << moleSet[i].mapRet.geneMapPosition[j].second << "\t" << moleSet[i].mapRet.alignLenNum[j].first.len << "\t" << moleSet[i].mapRet.alignLenNum[j].second.len << "\n";
        }
    }
}
