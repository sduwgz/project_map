#include "mole.h"

#include <stdlib.h>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <cassert>

#include "SplitString.h"
using namespace std;

Mole& Mole::operator = (const Mole & o) {
    if(&o != this) {
        id = o.id;
        len = o.len;
        enzymeDisNum = o.enzymeDisNum;
        dis.assign(o.dis.begin(), o.dis.end());
        pos.assign(o.pos.begin(), o.pos.end());
    }
    return *this;
}

Mole Mole::reverseMole() {
    Mole revM;
    revM.id = -id;
    revM.len = len;
    revM.enzymeDisNum = enzymeDisNum;
    vector < int > tmp;
    tmp.assign(dis.begin(), dis.end());
    reverse(tmp.begin(), tmp.end());
    revM.dis.assign(tmp.begin(), tmp.end());
    return revM;
}
/*
   bool Mole::operator != (const Mole& o) const {
   return id != o.id;
   }*/


bool Mole::getDis() {
    assert(pos.size()>0);
    dis.resize(pos.size()-1); 
    int j =0;
    for (int i=0; i+1 < pos.size(); ++i) {
        int disOne = pos[i+1]-pos[i]; 
        if (disOne != 0) 
            dis[j++] = disOne;    
    }
    dis.resize(j);
    return true;      
}

void Mole::print() {
     cout<<"[mole Index] "<<id<< " pos size"<<pos.size()<<" dis size" << dis.size() << endl;
     cout<<"[mole position] "<<endl;
     for (int i=0; i<pos.size(); i++)
         cout<<pos[i]<<"\t";
     cout<<endl;
     cout<<"[mole dis]" <<endl;
     for (int i=0; i<dis.size(); i++)
         cout<<dis[i]<<"\t";
     cout<<endl;
}

bool MoleReader::read(Mole& mole) {
    enum {
        eId,
        ePos,
    };

    if(stream) {
        int state = eId;
        mole.mole_reset();
        string line;
        vector<double> tmp;
        while (std::getline(stream,line)) {
            if (state == eId) {
                tmp = SplitString(line).split2Dbl("\t ,");
                if (static_cast<int>(tmp[0])==0) {
                    mole.id = static_cast<int>(tmp[1]); 
                    state = ePos;
                } else {
                    std::cerr<<"input.bin=>invalid line for not start with 0: "<<line<<std::endl;
                    return false;
                }
            } else if (state = ePos) {
                tmp = SplitString(line).split2Dbl("\t ,");
                if (static_cast<int>(tmp[0])==1) {
                    mole.pos.resize(tmp.size()-1,0);    
                    for (int i = 1; i < tmp.size(); i++) { 
                         mole.pos[i-1] = static_cast<long>(tmp[i]); 
                    }
                    if (mole.pos.size() > 0) {
                        mole.getDis();
                    }
                    std::getline(stream,line);
                    std::getline(stream,line);
                    state = eId;
                    return true;  
                } else {
                    std::cerr<<"input.bin=>invalid line for not start with 1: "<<line<<endl;
                    return false;
                }
            }
        }
        return false;
    }
    return false;

}
