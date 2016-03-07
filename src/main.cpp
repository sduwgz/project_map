// Last Update:2015-11-13 15:28:03
/**
 * @file main.cpp
 * @brief 
 * @author Yanbo Li, liyanbo@ict.ac.cn
 * @version 0.1.00
 * @date 2015-10-22
 */


#include <cstdarg>
#include <climits>

#include "mole.h"
#include "gene.h"
#include "map.h"

#include "SplitString.h"

static const char *optStr = "o:c:m:g:";

int printHelp() {
    std::cout << "USAGE : ./Nano-ARCS -g [input gene file name] -m [input mole file name][other options]" << std::endl;
    std::cout << "-g\tgene file name" << std::endl;
    std::cout << "-m\tmole file name" << std::endl;
    std::cout << "-c\tmin length to map" << std::endl;
    std::cout << "-o\toutput directory[default .]" << std::endl;
    std::cout << std::endl;
    return 1;
}

int main(int argc, char ** argv) {
    if (argc < 4) {
        return printHelp();
    }

    std::string geneFileName;
    std::string moleFileName;
    int MINCNT = 3;
    std::string outPrefix = "../whole_map/";

    int opt = -1;
    while ((opt = getopt(argc, argv, optStr)) != -1) {
        switch(opt){
            case 'g': geneFileName = optarg; break;
            case 'm': moleFileName = optarg; break;  
            case 'c': MINCNT = atoi(optarg); break;
            case 'o': outPrefix = optarg; break;
            default: printHelp(); 
        }
    }
    
    //std::cout<< "[INFO]" << "Read whole gene. File name: " << geneFileName << std::endl;
    //load ref_file
    Gene g;
    {
        ifstream geneIn(geneFileName.c_str());
        GeneReader gReader(geneIn);
        if (!gReader.read(g)) {
            //std::cout << "[REPORT]" << "The gene has been inited." << std::endl;
        //else{ 
            std::cout << "[REPORT]" << "The gene is bad." << std::endl;
            exit(EXIT_FAILURE);
        }
        //std::cout<< "[INFO]" << "Read mole. File name: " << moleFileName << std::endl;
    }
    //load mole_file
    vector<Mole> moleSet;
    {
        Mole m;
        ifstream moleIn(moleFileName.c_str());
        MoleReader mReader(moleIn);
       
       while (mReader.read(m)) {
            moleSet.push_back(m);
            Mole revM = m.reverseMole();
            moleSet.push_back(revM);
        }
        int mole_number = moleSet.size();
        if(mole_number == 0){
            std::cout << "[REPORT]" << "No mole is in moleSet." << std::endl; 
            exit(EXIT_FAILURE);
        }
    //std::cout << "[REPORT]" << mole_number << " moles have been inited." << std::endl;
    }

    double mean=46.0, variance=570;
    double beta = 0.15;
    double alpha = 0.15;
 
    Map maptool(mean, variance, alpha, beta, MINCNT,outPrefix);
    maptool.whole_map_score(moleSet, g.dis);
    
    return 1;
}

