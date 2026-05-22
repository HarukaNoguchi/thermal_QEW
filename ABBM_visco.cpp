#include <fstream>
#include <iostream>
#include <sstream>
#include <assert.h>
#include <math.h>
#include <unordered_set>
#include <chrono>
#include <random>
#include <functional>
#include <omp.h>
#include "ABBM_visco.hpp"
#include <filesystem>
#include <sys/stat.h>
#include <vector>
#include <string>
#include <iomanip>
#include <mutex>
#include <algorithm>
#include <limits>
#include <set>

using namespace std;

double monte_range(int seed,double min_val,double max_val)///熱活性化の判定に用いる
{
    // 乱数生成器
    thread_local static std::mt19937_64 mt64(seed);
    // [min_val, max_val] の一様分布整数 (real) の分布生成器
    std::uniform_real_distribution<> get_rand_uni_int(min_val, max_val);
    return get_rand_uni_int(mt64);
}


MD::MD(void)
{  lat = new line();
}
//------------------------------------------------------------------------
MD::~MD(void)
{delete lat;
}

double Gauss_dist(int seed){
    thread_local static std::mt19937_64 mt64(seed);
    std::normal_distribution<> dist(0,0.3);
    return dist(mt64);
}



void MD::makeini(int seed) {
    //std::ifstream ifs1("/home/noguchi/creep_modified/QEW_G_10_V_0.0000010000_N_1024_dt_0.1000000000_c_0.300000_G_10.000000_D_1.000000/T_3.000000/front_seed1.dat");
    std::ostringstream oss;
    // 初期状態を平坦な状態から始める場合
    
    int h=int(L/du);
    for (int i=0;i<N;++i){
        (lat->add_site)(0.0,0,h,seed); ///  u=0,w=0から始める
    }


}

double MD::random_force(int I){
    double u_I = (lat->site)[I].u;
    double shifted_u = u_I + 10.0;

    const auto& noise = (lat->site)[I].noise;
    double umax = (noise.size() - 1) * du;

    if (shifted_u < 0.0) shifted_u = 0.0;
    if (shifted_u > umax) shifted_u = umax;

    int i = int(floor(shifted_u / du));

    if (i < 0) i = 0;
    if (i + 1 >= noise.size()) i = noise.size() - 2;

    double x = shifted_u - i * du;
    double y = (i + 1) * du - shifted_u;

    double f0 = noise[0];
    return ( (noise[i+1] - noise[i])*x ) / du + noise[i];
}

double MD::thermal_noise(int seed){
    thread_local static std::mt19937_64 mt64(seed);
    std::normal_distribution<> dist(0.0,1.0);
    return sqrt(2.0 * gamma * T / dt) * dist(mt64);
}




void MD::Eular(double u_av,double t,int i,int seed){
    double elast_term=(c)*((lat->site)[((i+1)%N)].u+ (lat->site)[((i-1)%N+N)%N].u-2*((lat->site)[i].u));
    double u_new=(lat->site)[i].u+(dt)*((k*(V*t-(lat->site)[i].u)) + elast_term  
    +  std::sqrt(2*D)*random_force(i) +std::sqrt(2*T*dt)*thermal_noise(seed));
    (lat->site)[i].u=u_new;
}



/*void MD::one_step(int n,std::ofstream& ofs,int seed){
    double u_av=0;
    double t=dt*n;

    for (int i=0;i<N;i++){
        u_av+=(lat->site)[i].u/N;
    }

    const double eps = du;
    if (u_av>5 && (u_av - u_old) > eps && u_av<1000 ){
        ofs<<n<<" "<<u_av<<" "<<k*(V*t-u_av)<<" "<<(u_av - u_old)<<endl;
        u_old=u_av;

    }
    
    for (int i=0;i<N;i++){
        Eular(u_av,t*dt,i,seed);
    }
}
*/



std::string to_string_with_precision(double value, int precision) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << value;
    return oss.str();
}


void MD::run(int seed)
{   
    ///make file///
    int steps=1e8;
    std::string filename1="./V_"+to_string_with_precision(V, 10) +"_N_"+std::to_string(N)+"_gamma_"+to_string_with_precision(gamma,6)+
    "_c_"+to_string_with_precision(c,6)+"_D_"+to_string_with_precision(D,6)+"_ck"+to_string_with_precision(k,6)+"_T_"+to_string_with_precision(T,6);
    cout<<filename1<<endl;
    
    mkdir((filename1).c_str(),0777);
    
    std::string filename2=filename1+"/T_"+std::to_string(T);
    mkdir((filename2).c_str(),0777);

    string SD_all2 = "./" + filename2 + "/SD_all2_"+std::to_string(seed)+".dat";
    //ofstream ofs_SD_all2_damy(SD_all2);
    string velocity = "./" + filename2 + "/_velocity_"+std::to_string(seed)+".dat";
    ofstream ofs_velocity(velocity);


    makeini(seed);



    for (int t=0;t<steps;t++){
        one_step(t,ofs_velocity,seed);
    }

    string front = "./" + filename2 + "/front_seed"+std::to_string(seed)+".dat";
    ofstream ofs_front(front);
    //ofs_front<<"index"<<" "<<"z"<<" "<<"pin_position"<<" "<<"th"<<" "<<"n_redistributerd"<<" "<<"distance_for_pin_position"<<" "<<"pinned_or_not"<<endl;
    for(int i=0;i<N;i++){
        ofs_front<<std::setprecision(20)<<i<<" "<<(lat->site)[i].u<<endl;
    }

};