#ifndef mole_h_
#define mole_h_

#include <vector>
#include <string>
#include <iostream>

#include "constant.h"

typedef std::vector< std::pair< int , int > > FragmentLength;
typedef std::vector< std::pair< std::string , std::string > > FragmentCond;
typedef std::vector< std::pair< int , int > > MapPosition;
typedef std::pair < int , int > AlignPosition;

struct MapRet {   
    double score;
    FragmentLength alignFragmentLength;
    FragmentCond alignFragmentCond;
    AlignPosition alignStartPosition;
    AlignPosition alignEndPosition;
    MapPosition moleMapPosition;
    MapPosition geneMapPosition;
};

class Mole {
public:
    Mole() {} 
    explicit Mole(std::string id) : _id(id) {
    }
    virtual ~Mole() {}
    bool getDistance(); 
    Mole reverseMole();    
public:
    std::string _id;
    std::vector < long > _position;
    std::vector < int > _distance;
    std::vector < int > _mapPosition;
    MapRet mapRet;
    friend class MoleReader;
};

class MoleReader {
public:
    MoleReader(std::istream& stream) : _stream(stream) {};
    bool read(Mole& mole);
    void reset(Mole& mole) {
        mole._id = "";
        mole._distance.clear();
        mole._position.clear();
        mole._mapPosition;
    }

private:
    std::istream& _stream;
};
#endif //mole_h_
