#ifndef mole_h__
#define mole_h__

#include <vector>
#include <iostream>

#include "constant.h"

struct LenNum {
    int len;
    int num;
};
struct MapRet {   
    bool label = false;
    double score;
    std::vector< std::pair<LenNum,LenNum> > alignLenNum;
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
    void reset() {
        _id = -1;
        _position.clear();
        _distance.clear();
        resetMapRet();
    }
    void resetMapRet() {
        mapRet.label = false;
        mapRet.score = INIT_SCORE;
        mapRet.alignLenNum.clear();
        mapRet.alignGenePosition = std::make_pair(-1, -1);
        mapRet.alignMolePosition = std::make_pair(-1, -1);
        mapRet.moleMapPosition.clear();
        mapRet.geneMapPosition.clear();
    }
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
