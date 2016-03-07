// Last Update:2015-11-10 17:00:32
/**
 * @file gene.cpp
 * @brief 
 * @author Yanbo Li, liyanbo@ict.ac.cn
 * @version 0.1.00
 * @date 2015-10-21
 */

#include "gene.h"

#include <stdlib.h>
#include <numeric>
#include <iostream>
#include <cassert>

#include "SplitString.h"

bool Gene::getDistance() {
    distance.resize(position.size() - 1); 
    for (int i = 0; i < position.size() - 1; ++ i) {
        distance[i] = position[i + 1] - position[i];
    }
    return true;      
}

bool GeneReader::read(Gene& gene) {
    if(!_stream) {
        return false;
    }
    gene.position.clear();
    gene.distance.clear();
    std::string line;
    std::vector< double > data;
    while (std::getline(_stream, line)) {
        if(line[0] == '#') continue;
        data = SplitString(line).split2Dbl("\t ,");
        if (static_cast< int > (data[0]) == 1) {
            gene.position.push_back(static_cast< int > (data[5]));
        } else {
            std::cerr << "Ref=>invalid line: " << line << std::endl;
            return false;
        }
    }
    gene.getDistance();
    return true;
}
