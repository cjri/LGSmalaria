//emission.h

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
#include "gsl/gsl_sf_gamma.h"

using namespace std;


class Emission{
public:
  Emission();
  void set(int nsamples, vector<int>& nsites, int grid);
  ~Emission();
  void set_dist();
  int dist_set;
  gsl_matrix ** EmitProb;
  int EmitProb_set;
  double rnd_emit, shrink;
  int mode;
  void set_EmitProb();
  int nSamples;
  int gridSize;
  double dx;
  unsigned int ** reads;
  unsigned int ** depths;
  unsigned int ** loci;
  double ** dist;
  double * xgrid;
  int * nSites;
  int total_loci;
};
