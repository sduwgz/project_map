#ifndef GENE_H
#define GENE_H

#include <vector>
#include <iostream>

class Gene {
public:
    typedef std::vector < long > Position;
    typedef std::vector < int > Distance;
    Gene() {}
    virtual ~Gene() {}
    bool getDistance();
public:  
    Position position;
    Distance distance;
    friend class GeneReader;
};

class GeneReader {
public:
    GeneReader(std::istream& stream) : _stream(stream) {};
    bool read(Gene& gene);

private:
    std::istream& _stream;
};
#endif  /*GENE_H*/
