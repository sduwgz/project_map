#ifndef mole_h__
#define mole_h__

#include <vector>
#include <cfloat>
#include <iostream>
using namespace std;
using namespace std;
struct LenNum {
    int len;
    int num;
};
struct MapRet {   
    bool label = false;//全局map上没有
    double score;
    vector< pair<LenNum,LenNum> > alignLenNum;//分子上某一段长度和数量对应基因上的某段长度和数量
    pair <int,int> alignMolePosition;//局部联配分子联配的起点和终点
    pair <int,int> alignGenePosition;//分子联配到基因上的起点和终点
    vector< pair<int,int> > moleMapPosition;
    vector< pair<int,int> > geneMapPosition;
};

class Mole {
public:
    Mole() {} 
    explicit Mole(size_t _id) : id(_id) {
    }
    virtual ~Mole() {}
    bool getDis(); 
    //bool moleToFLES();
    int get_enzyme_dis_num() const{ return dis.size(); }   
    int get_id() const{ return id; }
    Mole& operator = (const Mole& o);
    Mole reverseMole();    
    void print();
    void mole_reset() {
        id = -1;
        pos.clear();
        dis.clear();
        mapRet_reset();
    }
    void mapRet_reset() {
        mapRet.label = false;
        mapRet.score = -DBL_MAX;
        mapRet.alignLenNum.clear();
        mapRet.alignGenePosition = make_pair(-1,-1);
        mapRet.alignMolePosition = make_pair(-1,-1);
        mapRet.moleMapPosition.clear();
        mapRet.geneMapPosition.clear();
    }
    friend std::ostream& operator<<(std::ostream& os, const Mole& component) ;
    MapRet mapRet;
public:
    int id;
    std::vector < long > pos;
    
    std::vector < int > dis;  //distance of Enzyme
    vector<int> mapping_pos;  //position on ref 
    int len;
    int enzymeDisNum;
    friend class MoleReader;
};

class MoleReader {
public:
    MoleReader(std::istream& _stream) : stream(_stream) {};
    bool read(Mole& mole);

private:
    std::istream& stream;
};

// component Format Specification : file name :component_iter
// Syntax
//    <component_id> :=    [0-9]+
//    <contig_id>    :=    [0-9]+ ... [0-9]+
//    <gap>          :=    [0-9]+ ... [0-9]+
// Requirements
//    1. The component id appears right after '>component' 
//    2. The number of gap must equal the number of cotigs_id minus 1


#endif //mole_h_
