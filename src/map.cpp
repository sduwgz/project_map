// Last Update:2015-11-19 22:16:48
/**
 * @file map.cpp
 * @brief 
 * @author Yanbo Li, liyanbo@ict.ac.cn
 * @version 0.1.00
 * @date 2015-10-21
 */
#include "map.h"
#include "Types.h"
#include "mole.h"
bool Map::local_map_score(vector<Mole>& moleSet,vector <int>& gene) {
    cout<<"[Info] : begin local map score:"<<endl;
    int iter = 0;
    double totalScore;
    while (true) {
        cout << iter << "th iteration " << "mu: " << mu << " sigma: " <<sigma <<" alpha: " << alpha << " beta: " << beta << endl; 
        cerr << iter << "th iteration " << "mu: " << mu << " sigma: " <<sigma <<" alpha: " << alpha << " beta: " << beta << endl; 

        int right_count = 0; 
        totalScore = 0.0;
        for (int i=0; i<moleSet.size(); i=i+2) {
            assert(moleSet[i].id == -moleSet[i+1].id);
            moleSet[i].mapRet_reset();
            moleSet[i+1].mapRet_reset();
            if (moleSet[i].dis.size()<MINCNT) continue;
            local_DP_score(moleSet[i],gene);
            local_DP_score(moleSet[i+1],gene);
            if (moleSet[i].mapRet.label || moleSet[i+1].mapRet.label) {
                if (moleSet[i].mapRet.label && moleSet[i+1].mapRet.label) {
                    totalScore += double(max(moleSet[i].mapRet.score,moleSet[i+1].mapRet.score));
                } else if (moleSet[i].mapRet.label) {
                    totalScore += double(moleSet[i].mapRet.score);
                } else {
                    totalScore += double(moleSet[i+1].mapRet.score);
                }
                right_count ++;
            } else {
                cerr << moleSet[i].id << " and " << moleSet[i+1].id << " can not map head to tail" << endl;
            }
        }
        cout<<"map number: "<<moleSet.size()<<" local map number "<<right_count<<endl;
        cout<<"total score: "<<totalScore<<endl;
        string filename = outPrefix + "mu_" + to_string(mu) + "sigma_" + to_string(sigma); 
        print_score(filename,moleSet); 
        if (change_parameter(moleSet)==false)
            break;
        iter ++; 
    }
    return true; 
}
bool Map::whole_map_score(vector<Mole>& moleSet,vector <int>& gene) {
    cout<<"[Info] : begin whole map score:"<<endl;
    int iter = 0;
    double totalScore;
    //while (true) {
        cout << iter << "th iteration " << "mu: " << mu << " sigma: " <<sigma <<" alpha: " << alpha << " beta: " << beta << endl; 
        cerr << iter << "th iteration " << "mu: " << mu << " sigma: " <<sigma <<" alpha: " << alpha << " beta: " << beta << endl; 

        int right_count = 0;
        totalScore = 0.0;
        for (int i=0; i<moleSet.size(); i=i+2) {
            assert(moleSet[i].id == -moleSet[i+1].id);
            moleSet[i].mapRet_reset();
            moleSet[i+1].mapRet_reset();
            if (moleSet[i].dis.size()<MINCNT) continue;
            whole_DP_score(moleSet[i],gene);
            whole_DP_score(moleSet[i+1],gene);
            if (moleSet[i].mapRet.label || moleSet[i+1].mapRet.label) {
                if (moleSet[i].mapRet.label && moleSet[i+1].mapRet.label) {
                    totalScore += double(max(moleSet[i].mapRet.score,moleSet[i+1].mapRet.score));
                } else if (moleSet[i].mapRet.label) {
                    totalScore += double(moleSet[i].mapRet.score);
                } else {
                    totalScore += double(moleSet[i+1].mapRet.score);
                }
                right_count ++;
            } else {
                cerr << moleSet[i].id << " and " << moleSet[i+1].id << " can not map head to tail" << endl;
            }
        }
        cout<<"map number: "<<moleSet.size()<<" whole map number "<<right_count<<endl;
        cout<<"total score: "<<totalScore<<endl;
        string filename = outPrefix + "mu_" + to_string(mu) + "sigma_" + to_string(sigma); 
        print_score(filename,moleSet); 
       // if (change_parameter(moleSet)==false) {   
       //     break;
       // }
        iter ++; 
    //}
    return 1;
}


/*bool Map::whole_map(vector<Mole>& moleSet,vector <int>& gene) {
  cout<<"[Info] : begin whole map without score:"<<endl;
  for (int i=0; i<moleSet.size(); i++) {
  if (moleSet[i].dis.size()<=MINCNT) continue;
  moleSet[i].mapRet_reset();
  for (int j=0; j<gene.size()-moleSet[i].dis.size(); j++) {
  whole_DP(moleSet[i],gene,j);
  }
  }
  string filename = outPrefix + "all possible position"; 
  print_score(filename,moleSet); 
  return true; 
  }

  bool Map::whole_DP() {

  return true;
  }*/



double Map::validScore(int a, int b) {
    if (a-b<mu+3*sigma && a-b>mu-3*sigma) { 
        return 100000*1/sqrt(2*3.14)/sigma*exp(-(a-b-mu)*(a-b-mu)/2/sigma/sigma); 
    } else {
        return -DBL_MAX;
    }
}

double Map::validScore(int moleB, int moleE, int geneB, int geneE, const vector<int> & mole, const vector<int> & gene) {
    assert(moleB<=moleE && geneB<=geneE && moleE<mole.size() && geneB<gene.size());
    int moleLen = 0, geneLen = 0;
    for(int i=moleB; i<=moleE; i++) {
        moleLen += mole[i];
    }
    for(int j=geneB; j<=geneE; j++) {
        geneLen += gene[j];
    }
    int delta = moleLen - geneLen;
    if(moleE - moleB != 0) {
        //Insert
        //cout << "DEBUG: insert:" <<endl;
        //cout << moleE - moleB + 1 << endl;
        //cout << pI(moleE - moleB + 1) << endl;

        return guss(delta) + pI(moleE - moleB + 1) - background(delta);
    }
    else if(geneE - geneB != 0) {
        //Delete
        int del = (int)((geneE - geneB + 0.0) / moleLen * 10000 + 0.5);
        if(del < 1) {
            del = 1;
        }
        if(del > 20) {
            del = 20;
        }
        return guss(delta) +  pD(del) - background(delta);
    }
    else {
        return guss(delta) - background(delta);
    }
    /*拟合得到的背景分布
        return log(1/sqrt(2*3.14)/sigma*exp(-(moleLen-geneLen-mu)*(moleLen-geneLen-mu)/2/sigma/sigma)*pow(beta,moleE-moleB)*pow(alpha,geneE-geneB)) 
            - log(1/sqrt(2*3.14)/509771*exp(-(moleLen-geneLen-1466)*(moleLen-geneLen-1466)/2/509771/509771)); 
*/
/*考虑delete
    return log(1/sqrt(2*3.14)/sigma*exp(-(moleLen-geneLen-mu)*(moleLen-geneLen-mu)/2/sigma/sigma)*(1-alpha)*pow(alpha,geneE-geneB)) 
        - log(1/sqrt(2*3.14)/10840*exp(-(moleLen-geneLen-1870)*(moleLen-geneLen-1870)/2/10840/10840)); 
*/
}
double Map::guss(int delta) {
    return 0.0 - 0.5 * log(2*3.14) - log(sigma) - (delta-mu)*(delta-mu)/2/sigma/sigma;
}

double Map::pD(int k) {
    double lambd = 1.82;
    long long fact = 1;
    for(int i = 1; i <= k; ++i) {
        fact *= i;
    }
    return k * log(lambd) - lambd - log(fact);
}

double Map::pI(int k) {
    double lambd = 1.064;
    return log(lambd) + ((0 - k) * lambd);
}

double Map::background(int delta){
    //background distribution from fit
    return log(1/sqrt(2*3.14)/10840*exp(-(delta-1870)*(delta-1870)/2/10840/10840));
}
bool Map::whole_DP_score(Mole& mole, vector<int>& gene) {

    double dp[200][700];
    pair<int,int> backTrack[200][700];
    assert(mole.dis.size()<=200 && gene.size()<=700); 
    assert(mole.mapRet.label == false && mole.mapRet.score == -DBL_MAX);
    vector<int> tempDis;
    vector<pair<int,int> > tempNum;
    for(int i=0; i<=mole.dis.size(); i++) {
        for(int j=0; j<=gene.size(); j++) {
            dp[i][j] = -DBL_MAX ;
            backTrack[i][j].first = -1;
            backTrack[i][j].second = -1;
        }
    }
    for (int j=0; j<=gene.size(); j++) {
        dp[0][j] = 0.0;
    }
    for(int i=1; i<=mole.dis.size(); i++) {
        for(int j=1; j<=gene.size(); j++) {
            for(int tj=(j-5>0)?j-5:0; tj<=j-1; tj++) {
                double temp = dp[i-1][tj] + validScore(i-1,i-1,tj,j-1,mole.dis,gene);
                if (temp > dp[i][j]) {
                    dp[i][j] = temp;
                    backTrack[i][j].first = i-1;
                    backTrack[i][j].second = tj;
                }
            }
            for(int ti=(i-5>0)?i-5:0; ti<=i-1; ti++) {
                double temp = dp[ti][j-1] + validScore(ti,i-1,j-1,j-1,mole.dis,gene);
                if (temp > dp[i][j]) {
                    dp[i][j] = temp;
                    backTrack[i][j].first = ti;
                    backTrack[i][j].second = j-1;
                }
            }
        }
    }

    double max = -DBL_MAX;

    for(int j=0; j<=gene.size(); j++) {
        if (dp[mole.dis.size()][j] > max) {
        max = dp[mole.dis.size()][j];
        mole.mapRet.alignMolePosition.second = mole.dis.size();
        mole.mapRet.alignGenePosition.second = j;
        }
    }
/*
    for (int i=0; i<=mole.dis.size(); i++) {
        for(int j=0; j<=gene.size(); j++) {
            cout << dp[i][j] << " " ;
            if (dp[i][j] > max) {
                max = dp[i][j];
                mole.mapRet.alignMolePosition.second = i;
                mole.mapRet.alignGenePosition.second = j;
            }
        }
        cout << endl;
    }*/
    int pi = mole.mapRet.alignMolePosition.second, pj = mole.mapRet.alignGenePosition.second;
    cerr << "mole" << mole.get_id() << endl; 
 //   cerr << pi << "\t" << pj << endl;
    if (max == -DBL_MAX) {
        cerr << "case 1 map failure" <<endl; 
        return false;
    }

    while(true) {
        int pii = pi;
        int pjj = pj;
        pi = backTrack[pii][pjj].first;
        pj = backTrack[pii][pjj].second;
   //     cerr << pi << "\t" << pj <<endl;
        if (pi == -1 && pj == -1) { 
            mole.mapRet.alignMolePosition.first = pii;
            mole.mapRet.alignGenePosition.first = pjj;
            break;
        }
        //应该只用正反链中分高的，而且是全局联配上的链来更新参数
        LenNum moleLn,geneLn;
        assert(pi==pii-1 || pj==pjj-1);
        moleLn.num = pii-pi;
        geneLn.num = pjj-pj;
        moleLn.len = accumulate(mole.dis.begin()+pi,mole.dis.begin()+pii,0);
        geneLn.len = accumulate(gene.begin()+pj,gene.begin()+pjj,0);
        cerr << pi << "\t" << pii-1 << "\t" << pj << "\t" <<pjj-1 << "\t" << moleLn.len <<"\t" << geneLn.len << endl; 
        mole.mapRet.alignLenNum.push_back(make_pair(moleLn,geneLn));
        mole.mapRet.moleMapPosition.push_back(make_pair(pi,pii-1));
        mole.mapRet.geneMapPosition.push_back(make_pair(pj,pjj-1));
    }
    reverse(mole.mapRet.moleMapPosition.begin(),mole.mapRet.moleMapPosition.end());
    reverse(mole.mapRet.geneMapPosition.begin(),mole.mapRet.geneMapPosition.end());
    cout << mole.id << "\t" << max << "\t" << mole.mapRet.alignGenePosition.first+1 << "\t" <<  mole.mapRet.alignGenePosition.second+1 << "\t" << mole.mapRet.alignMolePosition.first << "\t" <<  mole.mapRet.alignMolePosition.second << endl; 
    mole.mapRet.score = max;
    if (mole.mapRet.alignMolePosition.first == 0) {
        mole.mapRet.label = true; 
        //  if (mole.mapRet.score < 0) {
        /*cerr << mole.id << ": "; 
          for (int i=0; i<mole.dis.size(); i++) {
          cerr<<mole.dis[i]<<" ";
          }
          cerr << endl;
          for (int i=0; i<mole.mapRet.alignLenNum.size(); i++)
          cerr << mole.mapRet.alignLenNum[i].first.len <<" "<< mole.mapRet.alignLenNum[i].second.len << " " << mole.mapRet.alignLenNum[i].first.num <<" "<< mole.mapRet.alignLenNum[i].second.num << endl;*/
        //  }
        return true;
    }
    cerr << "case 2 map failure" <<endl; 
    return false;
}

bool Map::local_DP_score(Mole& mole, vector<int>& gene) {
    double dp[200][700];
    pair<int,int> backTrack[200][700];
    assert(mole.dis.size()<=200 && gene.size()<=700); 
    assert(mole.mapRet.label == false && mole.mapRet.score == -DBL_MAX);
    vector<int> tempDis;
    vector<pair<int,int> > tempNum;
    for(int i=0; i<=mole.dis.size(); i++) {
        for(int j=0; j<=gene.size(); j++) {
            dp[i][j] = 0;
            backTrack[i][j].first = -1;
            backTrack[i][j].second = -1;
        }
    }
    for(int i=1; i<=mole.dis.size(); i++) {
        for(int j=1; j<=gene.size(); j++) {
            for(int tj=(j-5>0)?j-5:0; tj<=j-1; tj++) {
                double temp = dp[i-1][tj] + validScore(i-1,i-1,tj,j-1,mole.dis,gene);
                if (temp > dp[i][j]) {
                    dp[i][j] = temp;
                    backTrack[i][j].first = i-1;
                    backTrack[i][j].second = tj;
                }
            }
            for(int ti=(i-5>0)?i-5:0; ti<=i-1; ti++) {
                double temp = dp[ti][j-1] + validScore(ti,i-1,j-1,j-1,mole.dis,gene);
                if (temp > dp[i][j]) {
                    dp[i][j] = temp;
                    backTrack[i][j].first = ti;
                    backTrack[i][j].second = j-1;
                }
            }
        }
    }

    double max = 0;

    for (int i=0; i<=mole.dis.size(); i++) {
        for(int j=0; j<=gene.size(); j++) {
            if (dp[i][j] > max) {
                max = dp[i][j];
                mole.mapRet.alignMolePosition.second = i;
                mole.mapRet.alignGenePosition.second = j;
            }
        }
    }
    int pi = mole.mapRet.alignMolePosition.second, pj = mole.mapRet.alignGenePosition.second;
    if (max == 0) {
        cerr << "case 1 map failure" <<endl; 
        return false;
    }

    while(true) {
        int pii = pi;
        int pjj = pj;
        pi = backTrack[pii][pjj].first;
        pj = backTrack[pii][pjj].second;
        if (pi == -1 && pj == -1) { 
            mole.mapRet.alignMolePosition.first = pii;
            mole.mapRet.alignGenePosition.first = pjj;
            break;
        }
        //应该只用正反链中分高的，而且是全局联配上的链来更新参数
        LenNum moleLn,geneLn;
        assert(pi==pii-1 || pj==pjj-1);
        moleLn.num = pii-pi;
        geneLn.num = pjj-pj;
        moleLn.len = accumulate(mole.dis.begin()+pi,mole.dis.begin()+pii,0);
        geneLn.len = accumulate(gene.begin()+pj,gene.begin()+pjj,0);
        mole.mapRet.alignLenNum.push_back(make_pair(moleLn,geneLn));
    }

    cout << mole.id << "\t" << max << "\t" << mole.mapRet.alignGenePosition.first+1 << "\t" <<  mole.mapRet.alignGenePosition.second+1 << "\t" << mole.mapRet.alignMolePosition.first << "\t" <<  mole.mapRet.alignMolePosition.second << endl; 
    mole.mapRet.score = max;
    mole.mapRet.label = true; 

    return true;

}

void Map::remove_noise(vector<Mole>& moleSet,vector<int>& gene) {
    int count = 0;
    cerr<<"remove mole: " <<endl;
    for (int i=0; i<moleSet.size(); i++) {
        moleSet[i].mapRet_reset();
        if (moleSet[i].dis.size()<=MINCNT) {
            cerr << moleSet[i].id<<endl;
            moleSet.erase(moleSet.begin()+i);
            count++;
            i--;
            continue;
        }
        if (!whole_DP_score(moleSet[i],gene)) {
            //if (moleSet[i].mapRet.score <-1000) {
            cerr << moleSet[i].id<<endl;
            count ++;
            moleSet.erase(moleSet.begin()+i);
            i--;
        }
    }
    cout<<"noise number: "<<count<<endl;
    cout<<"right mole number:" <<moleSet.size()<<endl;
}

bool Map::change_parameter(vector<Mole> & moleSet) {

    double new_mu,new_sigma;
    double new_alpha,new_beta;
    mapDis.clear();
    mapNum.clear();
    vector<int> tempDis;
    vector<pair<int,int> > tempNum;
    int betterId;
    for (int i=0; i<moleSet.size(); i=i+2) {
        assert(moleSet[i].id == -moleSet[i+1].id);
        if (moleSet[i].mapRet.label==true && (moleSet[i+1].mapRet.label==false || 
                                                  (moleSet[i+1].mapRet.label == true && moleSet[i].mapRet.score >= moleSet[i+1].mapRet.score))) {
            betterId = i;
        } else if (moleSet[i+1].mapRet.label==true && (moleSet[i].mapRet.label==false || 
                                                           (moleSet[i].mapRet.label == true && moleSet[i+1].mapRet.score >= moleSet[i].mapRet.score))) {                betterId = i+1;
        } else {
            continue;
        }
        tempDis.resize(moleSet[betterId].mapRet.alignLenNum.size());
        tempNum.resize(moleSet[betterId].mapRet.alignLenNum.size());
        int count = 0;
        for (int j=0; j<tempDis.size(); j++) {
            int dis = moleSet[betterId].mapRet.alignLenNum[j].first.len - moleSet[betterId].mapRet.alignLenNum[j].second.len; 
            tempDis[count] = dis; 
            tempNum[count] = make_pair(moleSet[betterId].mapRet.alignLenNum[j].first.num,moleSet[betterId].mapRet.alignLenNum[j].second.num);
            count ++;   
        }
        tempDis.resize(count);
        tempNum.resize(count);
        mapDis.insert(mapDis.end(),tempDis.begin(),tempDis.end());
        mapNum.insert(mapNum.end(),tempNum.begin(),tempNum.end());
    }
    if (mapDis.size() == 0) {
        cout << "error" << endl;
        return false;
    }

    double sum=0.0;
    for (int i=0; i<mapDis.size(); i++) {
        sum += mapDis[i]; 
    }
    new_mu = sum/mapDis.size();
    sum = 0;
        for (int i=0; i<mapDis.size(); i++) {
            sum += (mapDis[i]-new_mu)*(mapDis[i]-new_mu); 
        }
        new_sigma = sqrt(sum/(mapDis.size()-1));
        int second =0,first=0;
        int firstTotal=0,secondTotal=0;
        for (int i=0; i<mapNum.size(); i++) {
            firstTotal += mapNum[i].first;
            secondTotal += mapNum[i].second;
            if (mapNum[i].first != 1) {
                assert(mapNum[i].second == 1);
                first += mapNum[i].first -1;      
            }
            if (mapNum[i].second != 1) {
                assert(mapNum[i].first == 1);
                second += mapNum[i].second -1;      
            }
        }
        new_alpha = double(second)/secondTotal;
        new_beta = double(first)/firstTotal;
        assert(new_alpha<1.0 && new_beta<1.0);

        cout<<"mu: "<< new_mu <<" sigma: "<<new_sigma<<" alpha: "<< new_alpha<<" beta: "<<new_beta<<endl;
        if (new_sigma == 0 ) {
            cout<<"new_sigma ==0"<<endl; 
            return false;
        }

        if (fabs(new_mu-mu)<1&&fabs(new_sigma-sigma)<1&&fabs(alpha-new_alpha)<0.1&&fabs(beta-new_beta)<0.1) {
            return false;
        }

        mu = new_mu;
        sigma = new_sigma;
        alpha = new_alpha;
        beta = new_beta;
        return true;;
    }

void Map::print_score(const string filename, const vector< Mole >& moleSet) {
        ofstream out;
        out.open(filename.c_str());

        for (int i=0; i<moleSet.size(); i++) {
            // if (moleSet[i].mapRet.label == true) {
            out << moleSet[i].get_id() <<"\t" << moleSet[i].mapRet.label <<"\t" << moleSet[i].mapRet.score << "\t" 
                <<  moleSet[i].mapRet.alignMolePosition.first << "\t" << moleSet[i].mapRet.alignMolePosition.second << "\t"  <<  moleSet[i].mapRet.alignGenePosition.first << "\t" << moleSet[i].mapRet.alignGenePosition.second <<endl;
            for (int j=0; j<moleSet[i].mapRet.moleMapPosition.size(); j++) {
                out << moleSet[i].mapRet.moleMapPosition[j].first << "\t" << moleSet[i].mapRet.moleMapPosition[j].second << "\t" << moleSet[i].mapRet.geneMapPosition[j].first << "\t" << moleSet[i].mapRet.geneMapPosition[j].second << "\t";  
            }
            out<<endl;

            //  }
        
        }
    }

void Map::get_background_distribution() {

    }

