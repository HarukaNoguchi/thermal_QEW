#include <iostream>
#include <assert.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <random>
#include <omp.h> 
#include "particle.hpp"
#include <string>
#include <sys/stat.h>
using namespace std;

double line::Gauss(int seed){
    thread_local static std::mt19937_64 mt64(seed);
    std::normal_distribution<> dist(0,1);
    return dist(mt64);
}


void line::add_site(double v,double w,int L,int seed){
    Site a;
    a.u=v;
    a.w=w;
    for (int i=0;i<L;i++){
        a.noise.push_back(Gauss(seed));
    }
    site.push_back(a);//push_back(a)はvectorのメンバ関数
};
