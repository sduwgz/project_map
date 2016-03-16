#ifndef mole_h_
#define mole_h_

#include <vector>
#include <iostream>

#include "constant.h"

typedef std::vector< std::pair< int , int > > FragmentLength;
typedef std::vector< std::pair< int , int > > MapPosition;
typedef std::pair < int , int > AlignPosition;

struct MapRet {   
    double score;
    FragmentLength alignFragmentLength;
    AlignPosition alignStartPosition;
    AlignPosition alignEndPosition;
    MapPosition moleMapPosition;
    MapPosition geneMapPosition;
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
