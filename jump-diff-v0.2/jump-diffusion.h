//jump-diffusion.h

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <math.h>
#include <ctype.h> 
#include <string>
#include <map>
#include <vector>

// GSL headers...
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_multimin.h"

using namespace std;

class JumpDiffusion{
public:
  JumpDiffusion(Emission * emit);
  ~JumpDiffusion();
  Emission * myEmit;
  int nSamples;
  double sigma, jump;
  int FwdBwd_done, wTotal;
  int gridSize;
  double dx;
  int * nSites;
  double ** dist;
  unsigned int ** loci;
  double ** pstay;
  double ** pjump;
  double ** pnojump;
  double * xgrid;
  void get_EmitProb(int read, int depth, double * xgrid, gsl_vector * eprob);// emission probability
  gsl_vector * proposal;
  gsl_matrix ** alpha;
  gsl_matrix ** gamma;
  gsl_matrix ** total;
  void do_Fwd();
  void do_Bwd();
  void set_DiffProp(gsl_matrix * propagator, double variance);
  double total_llh;
};

struct fpar{
  JumpDiffusion * myJD;
  vector<int> to_opt;
};

double find_parameter(double * param, JumpDiffusion * myJD, double low_stop, double high_stop);
double Qall( const gsl_vector * x, void * p);
