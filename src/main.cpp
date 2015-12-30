// Last Update:2015-11-13 15:28:03
/**
 * @file main.cpp
 * @brief 
 * @author Yanbo Li, liyanbo@ict.ac.cn
 * @version 0.1.00
 * @date 2015-10-22
 */


#include "mole.h"
#include "gene.h"
#include "map.h"

#include "SplitString.h"

#include <cstdarg>
#include <climits>

using namespace std;



static const char *optStr = "o:c:m:g:";
void printHelp() {
    cout << "USAGE : ./Nano-ARCS -g [input gene file name] -m [input mole file name][other options]" << endl;
    cout << "-g\tgene file name" << endl;
    cout << "-m\tmole file name" << endl;
    cout << "-c\tmin length to map" << endl;
    cout << "-o\toutput directory[default .]" << endl;
    cout << endl;
}

int main(int argc, char ** argv)
{
    string geneFileName;
    string moleFileName;

    int opt;
    opt = getopt(argc, argv, optStr);
    if(argc < 4) {
        printHelp();
        exit(EXIT_FAILURE);
    }
    int K;
    int MINCNT = 3;
    string outPrefix = "../whole_map/";

    while(opt != -1){
        switch(opt){
            case 'g': geneFileName = optarg; break;
            case 'm': moleFileName = optarg; break;  
            case 'c': MINCNT = atoi(optarg); break;
            case 'o': outPrefix = optarg; break;
            default: printHelp(); exit(EXIT_FAILURE);
        }
        opt = getopt(argc, argv, optStr);
    }
    
    cout<< "[INFO]" << "Read whole gene. File name: " << geneFileName << endl;
    Gene g;
    ifstream geneIn(geneFileName.c_str());
    GeneReader gReader(geneIn);
    if(gReader.read(g))
        cout << "[REPORT]" << "The gene has been inited." << endl;
    else{ 
        cout << "[REPORT]" << "The gene is bad." << endl;
        exit(EXIT_FAILURE);
    }
    cout<< "[INFO]" << "Read mole. File name: " << moleFileName << endl;

    Mole m;
    vector<Mole> moleSet;
    ifstream moleIn(moleFileName.c_str());
    MoleReader mReader(moleIn);
    while (mReader.read(m)) {
        moleSet.push_back(m);
        Mole revM = m.reverseMole();
        moleSet.push_back(revM);
    }
    int mole_number = moleSet.size();
    if(mole_number == 0){
        cout << "[REPORT]" << "No mole is in moleSet." << endl; 
        exit(EXIT_FAILURE);
    }
    cout << "[REPORT]" << mole_number << " moles have been inited." << endl;

    double mean=46.0, variance=570;
    double beta = 0.15;
    double alpha = 0.15;
 
    Map maptool(mean, variance, alpha, beta, MINCNT,outPrefix);
    maptool.whole_map_score(moleSet, g.dis);
    return 1;
}

