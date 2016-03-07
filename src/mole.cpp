#include "mole.h"

#include <stdlib.h>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <cassert>

#include "SplitString.h"

Mole Mole::reverseMole() {
    Mole reMole;
    reMole.id = -id;
    reMole.length = length;
    reMole.enzymeNumber = enzymeNumber;
    std::vector< int > data;
    data.assign(distance.begin(), distance.end());
    reverse(data.begin(), data.end());
    reMole.distance.assign(data.begin(), data.end());
    return reMole;
}
/*
   bool Mole::operator != (const Mole& o) const {
   return id != o.id;
   }*/


bool Mole::getDistance() {
    distance.resize(position.size() - 2); 
    for (int i = 0; i < position.size() - 2; ++ i) {
        distance[i] = position[i + 1] - position[i];
    }
    return true;
}

bool MoleReader::read(Mole& mole) {
    if(!_stream){
        return false;
    }
    enum {
        eId,
        ePosition,
    };

    int state = eId;
    mole.resetMole();
    std::string line;
    std::vector<double> data;
    while (std::getline(_stream, line)) {
        if (state == eId) {
            data = SplitString(line).split2Dbl("\t ,");
            if (static_cast<int> (data[0]) == 0) {
                mole.id = static_cast<int> (data[1]);
                state = ePosition;
            } else {
                return false;
            }
        } else if (state == ePosition) {
            data = SplitString(line).split2Dbl("\t ,");
            if (static_cast<int>(data[0]) == 1) {
                mole.position.resize(data.size() - 1,0);
                for (int i = 1; i < data.size(); ++ i) {
                     mole.position[i - 1] = static_cast<long> (data[i]);
                }
                if (mole.position.size() > 0) {
                    mole.getDistance();
                }
                std::getline(_stream,line);
                std::getline(_stream,line);
                state = eId;
                return true;
            } else {
                return false;
            }
        }
    }
}
