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

#include "mole.h" 


typedef std::vector< int > Fragment;

class Map {
public:
     Map(double mu, double sigma, double alpha, double beta, const std::string outPrefix) : _mu(mu), _sigma(sigma), _alpha(alpha), _beta(beta), _outPrefix(outPrefix) {};
     bool whole_map_score(std::vector<Mole>& moleSet, const std::vector<int>& gene) const;
     bool whole_DP_score(Mole& mole, const std::vector<int>& gene) const;
     void print_score(const std::string& filename, const std::vector< Mole >& moleSet) const;
     double validScore(const Fragment& moleFragment, const Fragment& geneFragment) const;
     
     double guss(int delta) const;
     double laplace(int delta) const;
     double pD(int siteNumber, int moleLen) const;
     double probInsertion(int k) const;
     double background(int delta) const;
private:
    double _mu;
    double _sigma;
    double _alpha;
    double _beta;
    std::string _outPrefix;
};
#endif  /*MAP_H*/
