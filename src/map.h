// Last Update:2015-10-29 09:58:59
/**
 * @file map.h
 * @brief 
 * @author Yanbo Li, liyanbo@ict.ac.cn
 * @version 0.1.00
 * @date 2015-10-21
 */

#ifndef MAP_H
#define MAP_H
#include <vector>
#include <string>
class Mole;
//#include "mole.h" 

using namespace std;
class Map {
private:
    double mu;
    double sigma;
    double alpha;
    double beta;
    int MINCNT;
    vector < int > mapDis;
    vector < pair < int,int > > mapNum;
    string outPrefix;//输出文件路径
public:
     Map(double _mu, double _sigma, double _alpha, double _beta, int _MINCNT, string _outPrefix) : mu(_mu), sigma(_sigma), alpha(_alpha), beta(_beta), MINCNT(_MINCNT),outPrefix(_outPrefix) {}; 
     bool local_map_score(vector<Mole>& moleSet,vector<int>& gene);
     bool whole_map_score(vector<Mole>& moleSet,vector<int>& gene); 
     bool local_DP_score(Mole& mole, vector<int>& gene);
     bool whole_DP_score(Mole& mole, vector<int>& gene);
     void print_score(const string filename, const vector< Mole >& moleSet); 
     bool change_parameter(vector<Mole> & moleSet);
     double validScore(int a, int b);
     double validScore(int moleB, int moleE, int geneB, int geneE, const vector<int> & mole, const vector<int> & gene); 
     double guss(int delta);
     double pD(int k);
     double pI(int k);
     double background(int delta);
     void remove_noise(vector<Mole>& moleSet,vector<int>& gene); 
     void get_background_distribution(); 
     /*public:
     bool local_map(vector<Mole>& moleSet,vector<int>& gene);
     bool whole_map(vector<Mole>& moleSet,vector<int>& gene); 
     bool local_DP(Mole& mole, vector<int>& gene);
     bool whole_DP(Mole& mole, vector<int>& gene);
*/

};
    


#endif  /*MAP_H*/
