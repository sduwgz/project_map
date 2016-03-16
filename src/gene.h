#ifndef gene_h_
#define gene_h_

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
    Position _position;
    Distance _distance;
    friend class GeneReader;
};

class GeneReader {
public:
    GeneReader(std::istream& stream) : _stream(stream) {};
    bool read(Gene& gene);

private:
    std::istream& _stream;
};
#endif  //gene_h_
