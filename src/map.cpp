/**
 * @file map.cpp
 * @author Yanbo Li, liyanbo@ict.ac.cn
 * @author Guozheng Wei, weiguozheng@ict.ac.cn
 * @date 2015-12-30
 */

#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/laplace.hpp>
#include <boost/math/distributions/exponential.hpp>

#include <math.h>

#include "map.h"
#include "Types.h"
#include "mole.h"

using namespace boost::math;
using namespace std;

bool Map::whole_map_score(vector<Mole>& moleSet,vector <int>& gene) {
    double totalScore = 0.0;
    int right_count = 0;
    for (int i = 0; i < moleSet.size(); i = i + 2) {
        moleSet[i].mapRet_reset();
        moleSet[i + 1].mapRet_reset();
        //MINCNT is a threshold of the min site in mole
        if (moleSet[i].dis.size() < MINCNT) continue;
        //dp
        whole_DP_score(moleSet[i],gene);
        whole_DP_score(moleSet[i + 1],gene);
        //between the reversed mole and the mole, choose the better one
        if (moleSet[i].mapRet.label || moleSet[i + 1].mapRet.label) {
            if (moleSet[i].mapRet.label && moleSet[i + 1].mapRet.label) {
                totalScore += double(max(moleSet[i].mapRet.score,moleSet[i+1].mapRet.score));
            } else if (moleSet[i].mapRet.label) {
                totalScore += double(moleSet[i].mapRet.score);
            } else {
                totalScore += double(moleSet[i + 1].mapRet.score);
            }
            ++ right_count;
        } else {
            cerr << moleSet[i].id << " and " << moleSet[i + 1].id << " can not map head to tail" << endl;
        }
    }
    string filename = outPrefix; 
    print_score(filename,moleSet); 
    return 1;
}

double Map::validScore(int moleB, int moleE, int geneB, int geneE, const vector<int> & mole, const vector<int> & gene) {
    //B for begin and E for end, moleB and moleE should less than mole.size()
    assert(moleB <= moleE && geneB <= geneE && moleE < mole.size() && geneB < gene.size());
    
    int moleLen = 0, geneLen = 0;
    //get the length of mole
    for(int i = moleB; i <= moleE; ++ i) {
        moleLen += mole[i];
    }
    //get the length of gene
    for(int j = geneB; j <= geneE; ++ j) {
        geneLen += gene[j];
    }
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
    
    int delta = moleLen - geneLen;
    
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
    if(moleE - moleB != 0 && geneE - geneB != 0){
        int siteNumber = geneE - geneB;
        return laplace1(delta) + pD(siteNumber, moleLen) + pI(moleE - moleB) - background(delta);
    }
    if(moleE - moleB != 0) {
        //Insert
        return laplace1(delta) + pI(moleE - moleB) + pD(0, moleLen) - background(delta);
    }
    else if(geneE - geneB != 0) {
        //Delete
        int siteNumber = geneE - geneB;
        return laplace1(delta) + pD(siteNumber, moleLen) + pI(0) - background(delta);
    }
    else {
        //Match
        return laplace1(delta)  - background(delta);
        //return laplace1(delta) + pI(0) + pD(0, moleLen) - background(delta);
    }
}

double Map::guss(int delta) {
    return 0.0 - 0.5 * log(2 * 3.14) - log(sigma) - (delta - mu) * (delta - mu) / (2 * sigma * sigma);
}
/*
double Map::laplace(int delta) {
    double mu = 0.0;
    double sigma = 570.0;
    double b = sigma / sqrt(2);
    int distance = delta - mu;
    if(distance >= 0) {
        int d = distance / 100;
        return ( (1 - 0.5 * exp(0 - (d + 1) * 100 / b ) ) - (1 - 0.5 * exp(0 - d * 100 / b ) ) );
    }
    else {
        int d = (0 - distance) / 100;
        return ( (0.5 * exp(0 - d * 100 / b) ) - (0.5 * exp(0 - (d + 1) * 100 / b ) ) );
        
    }
}
*/

double probDelete(int k) {
    double lambd = 0.225;
    boost::math::poisson_distribution<> p(lambd);
    return log(boost::math::pdf(p, k));
}

double Map::pD(int siteNumber, int moleLen) {
    
    int del = (int)((siteNumber + 0.0) / moleLen * 10000 + 0.5);
    if(del < 1) {
        del = 1;
    }
    if(del > 20) {
        del = 20;
    }
    return probDelete(del);

    //double lambd = 1.82;
    //long long fact = 1;
    
    //for(int i = 1; i <= del; ++i) {
    //    fact *= i;
    //}

    //return del * log(lambd) - lambd - log(fact) - 2;
}

/*
double Map::pI(int k) {
    double lambd = 1.064;
    return log(lambd) + ((0 - k) * lambd);
}
double Map::background(int delta){
    //background distribution from fit, mean = 1870, var = 10840^2
    return 0.0 - 0.5 * log(2 * 3.14) -log(10840) - (delta - 1870) * (delta - 1870) / (2 * 10840 * 10840);
}

double Map::guss(int delta) {
    
    double mu = 46.0;
    double sigma = 570.0;
    boost::math::normal_distribution<> n(mu, sigma);
    int distance = delta - mu;
    if(distance >= 0){
        int d = distance / 100;
        int interval_left = d * 100;
        int interval_right = (d + 1) * 100;
        return log(boost::math::cdf(n, interval_right) - boost::math::cdf(n, interval_left));
    }
    else{
        int d = distance / 100;
        int interval_left = (d - 1) * 100;
        int interval_right = d * 100;
        return log(boost::math::cdf(n, interval_right) - boost::math::cdf(n, interval_left));
        
    }
}

double pD(int k) {
    double lambd = 1.82;
    poisson_distribution<> p(lambd);
    return log(boost::math::pdf(p, k));
}

*/
double Map::pI(int k) {
    double lambd = 5;
    boost::math::exponential_distribution<> e(lambd);
    if(k == 0){
        return log(cdf(e, 0.5) - cdf(e, 0));
    }
        
    return log(cdf(e, k + 0.5) - cdf(e, k - 0.5));
}
double Map::background(int delta){
    //background distribution from fit, mean = 1870, var = 10840^2
    double mu = 1870.0, sigma = 10840.0;
    boost::math::normal_distribution<> n(mu, sigma);
    int distance = delta - mu;
    if(distance >= 0){
        int d = distance / 100;
        int interval_left = d * 100;
        int interval_right = (d + 1) * 100;
        return log(boost::math::cdf(n, interval_right) - boost::math::cdf(n, interval_left));
    }
    else{
        int d = distance / 100;
        int interval_left = (d - 1) * 100;
        int interval_right = d * 100;
        return log(boost::math::cdf(n, interval_right) - boost::math::cdf(n, interval_left));
    }
}
double Map::laplace1(int delta) {
    double mu = 46.0;
    double sigma = 403.0;
    boost::math::laplace_distribution<> l(mu, sigma);
    int distance = delta - mu;
    if(distance >= 0){
        int d = distance / 100;
        int interval_left = d * 100;
        int interval_right = (d + 1) * 100;
        return log(boost::math::cdf(l, interval_right) - boost::math::cdf(l, interval_left));
    }
    else{
        int d = distance / 100;
        int interval_left = (d - 1) * 100;
        int interval_right = d * 100;
        return log(boost::math::cdf(l, interval_right) - boost::math::cdf(l, interval_left));       
    }
}

bool Map::whole_DP_score(Mole& mole, vector<int>& gene) {
    int rows = mole.dis.size() + 1, cols = gene.size() + 1;
    double dp[rows][cols];
    pair<int,int> backTrack[rows][cols];
    vector<int> tempDis;
    vector<pair<int,int> > tempNum;
    for(int i = 0; i <= mole.dis.size(); ++ i) {
        for(int j=0; j <= gene.size(); ++ j) {
            dp[i][j] = -1000000.0;
            backTrack[i][j].first = -1;
            backTrack[i][j].second = -1;
        }
    }
    //the first row is inited to zero, and the first and second rows are inited to pI(1) and pI(2)
    for(int j = 0; j <= gene.size(); ++ j) {
        dp[0][j] = 0.0;
        dp[1][j] = pI(1);
        dp[2][j] = pI(2);
    }
    for(int i = 0; i <= mole.dis.size(); ++ i) {
        dp[i][0] = -3 * i;
    }
    for(int i = 1; i <= mole.dis.size(); ++ i) {
        for(int j = 1; j <= gene.size(); ++ j) {
            for(int tj = (j - 5 > 0) ? j - 5 : 0; tj <= j - 1; ++ tj) {
                double temp = dp[i - 1][tj] + validScore(i - 1, i - 1, tj, j - 1, mole.dis, gene);
                if (temp > dp[i][j]) {
                    dp[i][j] = temp;
                    backTrack[i][j].first = i - 1;
                    backTrack[i][j].second = tj;
                }
            }
            for(int ti = ( i - 5 > 0) ? i - 5 : 0; ti <= i - 1; ++ ti) {
                double temp = dp[ti][j - 1] + validScore(ti, i - 1, j - 1, j - 1, mole.dis, gene);
                if (temp > dp[i][j]) {
                    dp[i][j] = temp;
                    backTrack[i][j].first = ti;
                    backTrack[i][j].second = j - 1;
                }
            }
            if(i > 2 && j > 2){
                double temp = dp[i - 2][j - 2] + validScore(i - 2, i - 1, j - 2, j - 1, mole.dis, gene);
                if (temp > dp[i][j]) {
                    dp[i][j] = temp;
                    backTrack[i][j].first = i - 2;
                    backTrack[i][j].second = j - 2;
                }
            } 

            /*
            for(int ki = (i - 3 > 0) ? i - 3 : 0; ki < i; ++ ki) {
                double temp = dp[ki][j] + pI(j - ki);
                if (temp > dp[i][j]) {
                    dp[i][j] = temp;
                    backTrack[i][j].first = ti;
                    backTrack[i][j].second = j;
                }
            }
            for(int kj = (j - 3 > 0) ? j - 3 : 0; kj < j; ++ kj) {
                double temp = dp[i][kj] + pD(j - ki);
                if (temp > dp[i][j]) {
                    dp[i][j] = temp;
                    backTrack[i][j].first = i;
                    backTrack[i][j].second = tj;
                }
            }
            */
        }
    }

    double max = -1000000.0;

    for(int j = 0; j <= gene.size(); ++j) {
        //we should find the max score in the last row
        if (dp[mole.dis.size()][j] > max) {
            max = dp[mole.dis.size()][j];
            mole.mapRet.alignMolePosition.second = mole.dis.size();
            mole.mapRet.alignGenePosition.second = j;
        }
        //if the last one is a insertion, we must find the max score in last but one row, and give a punish
        double insertPunish1 = pI(1);
        if (dp[mole.dis.size() - 1][j] + insertPunish1 > max) {
            max = dp[mole.dis.size() - 1][j] + insertPunish1;
            mole.mapRet.alignMolePosition.second = mole.dis.size() - 1;
            mole.mapRet.alignGenePosition.second = j;
        }
        //if the last two is insertions
        double insertPunish2 = pI(2);
        if (dp[mole.dis.size() - 2][j] + insertPunish2 > max) {
            max = dp[mole.dis.size() - 2][j] + insertPunish2;
            mole.mapRet.alignMolePosition.second = mole.dis.size() - 2;
            mole.mapRet.alignGenePosition.second = j;
        }
    }
    
    int pi = mole.mapRet.alignMolePosition.second, pj = mole.mapRet.alignGenePosition.second;
    if (max == -1000000.0) {
        cerr << "case 1 map failure" <<endl; 
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

        moleLn.len = accumulate(mole.dis.begin() + pi, mole.dis.begin() + pii, 0);
        geneLn.len = accumulate(gene.begin() + pj, gene.begin() + pjj, 0);

        mole.mapRet.alignLenNum.push_back(make_pair(moleLn, geneLn));

        mole.mapRet.moleMapPosition.push_back(make_pair(pi, pii - 1));
        mole.mapRet.geneMapPosition.push_back(make_pair(pj, pjj - 1));
    }

    //reversed result is more human-facing
    reverse(mole.mapRet.moleMapPosition.begin(), mole.mapRet.moleMapPosition.end());
    reverse(mole.mapRet.geneMapPosition.begin(), mole.mapRet.geneMapPosition.end());
    reverse(mole.mapRet.alignLenNum.begin(), mole.mapRet.alignLenNum.end());

    cout << mole.id << "\t" << max << "\t" << mole.mapRet.alignGenePosition.first+1 << "\t" <<  mole.mapRet.alignGenePosition.second+1 << "\t" << mole.mapRet.alignMolePosition.first << "\t" <<  mole.mapRet.alignMolePosition.second << endl; 
    
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

void Map::print_score(const string filename, const vector< Mole >& moleSet) {
    ofstream out;
    out.open(filename.c_str());

    for (int i=0; i<moleSet.size(); i++) {
        out << moleSet[i].get_id() <<"\t" << moleSet[i].mapRet.label <<"\t" << moleSet[i].mapRet.score << "\t" 
            <<  moleSet[i].mapRet.alignMolePosition.first << "\t" << moleSet[i].mapRet.alignMolePosition.second << "\t"  <<  moleSet[i].mapRet.alignGenePosition.first << "\t" << moleSet[i].mapRet.alignGenePosition.second <<endl;
        for (int j=0; j<moleSet[i].mapRet.moleMapPosition.size(); j++) {
            out << moleSet[i].mapRet.moleMapPosition[j].first << "\t" << moleSet[i].mapRet.moleMapPosition[j].second << "\t" << moleSet[i].mapRet.geneMapPosition[j].first << "\t" << moleSet[i].mapRet.geneMapPosition[j].second << "\t" << moleSet[i].mapRet.alignLenNum[j].first.len << "\t" << moleSet[i].mapRet.alignLenNum[j].second.len << "\n";
        }
    }
}

void Map::get_background_distribution() {
}
