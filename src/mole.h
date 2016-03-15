#ifndef mole_h__
#define mole_h__

#include <vector>
#include <iostream>

#include "constant.h"

struct MapRet {   
    double score;
    std::vector< std::pair< int, int > > alignLenNum;
    std::pair < int, int > alignMolePosition;
    std::pair < int, int > alignGenePosition;
    std::vector< std::pair< int, int > > moleMapPosition;
    std::vector< std::pair< int, int > > geneMapPosition;
};

class Mole {
public:
    Mole() {} 
    explicit Mole(size_t id) : _id(id) {
    }
    virtual ~Mole() {}
    bool getDistance(); 
    Mole reverseMole();    
    MapRet mapRet;
public:
    int _id;
    std::vector < long > _position;
    std::vector < int > _distance;
    std::vector < int > _mapPosition;
    friend class MoleReader;
};

class MoleReader {
public:
    MoleReader(std::istream& stream) : _stream(stream) {};
    bool read(Mole& mole);

private:
    std::istream& _stream;
};
#endif //mole_h_
