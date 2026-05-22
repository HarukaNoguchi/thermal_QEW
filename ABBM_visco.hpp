#include <unordered_set>
#include <vector>
#include <omp.h> 
#include "particle.hpp"
class MD {
private:
  line *lat;
  void makeini(int);
  void Eular(double,double,int,int);
  void one_step(int,std::ofstream&,int);
  double random_force(int);
  double thermal_noise(int);


  public:
  MD(void);
  ~MD(void);
  void run(int);

  int seed;
  double gamma=1;
  double dt=0.01;
  double k=1;
  double D=100;
  double V=1e-4;
  double du=0.1;
  double c=1;
  double u_old=0;
  double T=0;
  int N=256;
  double L=10000;

  std::vector<double> break_h;
  
  std::string stabile_num="stability";
};
