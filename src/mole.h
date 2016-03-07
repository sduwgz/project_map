#ifndef mole_h__
#define mole_h__

#include <vector>
#include <cfloat>
#include <iostream>

struct LenNum {
    int len;
    int num;
};
struct MapRet {   
    bool label = false;
    double score;
    std::vector< std::pair<LenNum,LenNum> > alignLenNum;
    std::pair <int,int> alignMolePosition;
    std::pair <int,int> alignGenePosition;
    std::vector< std::pair<int,int> > moleMapPosition;
    std::vector< std::pair<int,int> > geneMapPosition;
};

class Mole {
public:
    Mole() {} 
    explicit Mole(size_t _id) : id(_id) {
    }
    virtual ~Mole() {}
    bool getDistance(); 
    //bool moleToFLES();
    Mole reverseMole();    
    void resetMole() {
        id = -1;
        position.clear();
        distance.clear();
        resetMapRet();
    }
    void resetMapRet() {
        mapRet.label = false;
        mapRet.score = -DBL_MAX;
        mapRet.alignLenNum.clear();
        mapRet.alignGenePosition = std::make_pair(-1,-1);
        mapRet.alignMolePosition = std::make_pair(-1,-1);
        mapRet.moleMapPosition.clear();
        mapRet.geneMapPosition.clear();
    }
    MapRet mapRet;
public:
    int id;
    std::vector < long > position;
    
    std::vector < int > distance;
    std::vector<int> mapPosition;
    int length;
    int enzymeNumber;
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
