// Last Update:2015-10-23 15:37:26
/**
 * @file gene.h
 * @brief 
 * @author Yanbo Li, liyanbo@ict.ac.cn
 * @version 0.1.00
 * @date 2015-10-21
 */

#ifndef GENE_H
#define GENE_H

#include <vector>
#include <iostream>
using namespace std;
class Gene {
public:
    Gene() {}
    virtual ~Gene() {}
    friend std::ostream& operator<<(std::ostream& os, const Gene& component) ;
    bool getDis();
    void print();
    //private:
public:  
    std::vector < long > pos;
    std::vector < int > dis;
    int len;
    friend class GeneReader;
};

class GeneReader {
public:
    GeneReader(std::istream& _stream) : stream(_stream) {};
    bool read(Gene& gene);

private:
    std::istream& stream;
};
#endif  /*GENE_H*/
