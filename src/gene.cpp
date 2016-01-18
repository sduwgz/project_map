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
using namespace std;
bool Gene::getDis() {
    assert(pos.size()>0);
    dis.resize(pos.size()-1); 
    int j = 0;
    for (int i = 0; i < pos.size() - 1; ++i) {
        int disOne = pos[i+1] - pos[i]; 
        if (disOne != 0) 
            dis[j++] = disOne;    
    }
    dis.resize(j);
    return true;      
}

bool GeneReader::read(Gene& gene) {
    if(stream) {
        gene.pos.clear();
        gene.dis.clear();
        string line;
        vector< double > tmp;
        while (std::getline(stream, line)) {
            if(line[0] == '#') continue;
            tmp = SplitString(line).split2Dbl("\t ,");
            if (static_cast< int >(tmp[0]) == 1) {
                gene.pos.push_back(static_cast< int > (tmp[5]));
            } else {
                std::cerr << "Ref=>invalid line: " << line << endl;
                return false;
            }
        }
        gene.getDis();
        return true;
    }
    return false;
}

void Gene::print() {
    cout<<"[gene pos]"<<endl;
    for (int i=0; i<pos.size(); i++) 
        cout<< pos[i] << "\t";
    cout<<endl;
    cout<<"[gene dis]"<<endl;
    for (int i=0; i<dis.size(); i++) 
        cout<< dis[i] << "\t";
    cout<<endl;
}
