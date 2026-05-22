#pragma once
#include <iostream>
#include <assert.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <chrono>
#include <omp.h> 
#include <mutex>


struct Site{
    double u=0;
    double w=0;
    std::vector<double> noise;
};

struct noise{
    double F=0;
};

class line{
    public:
    std::vector<Site> site;
    double Gauss(int);
    void add_site(double x,double w,int,int);
    //int number_of_atoms(void){return static_cast<int>(atoms.size());}
};