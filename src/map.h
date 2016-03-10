#ifndef MAP_H
#define MAP_H

#include <vector>
#include <string>
#include <map>
#include "mole.h" 


typedef std::vector< int > Fragment;
typedef std::vector< Mole > MoleSet;
typedef std::map<std::string, double> ParametersList;

class Map {
public:
     Map(){};
     bool initParameters(const std::string &parameter_file);
     bool run(MoleSet& moleSet, const std::vector<int>& gene) const;
     bool wholeDPscore(Mole& mole, const std::vector<int>& gene) const;
     void print_score(const std::string& filename, const std::vector< Mole >& moleSet) const;
     double validScore(const Fragment& moleFragment, const Fragment& geneFragment) const;
     
     double probLaplace(int delta) const;
     double probDeletion(int siteNumber, int moleLength) const;
     double probInsertion(int k) const;
     double probBackground(int delta) const;
private:
    ParametersList _parameters;
};
#endif  /*MAP_H*/
