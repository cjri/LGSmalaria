#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <map>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>


using namespace std;

struct rec {
	int chr;
	int index;
	int pos;
	int obs;
	double q;
	int N;
};

struct opt_data {
	int i;
	double x;
	double e;
	double rho1;
	double rho2;
	double rho3;
	double xf;
	double rx1;
	double rx2;
};


void ImportData (string input, vector<rec>& traj);

void SetXXVector(opt_data od, gsl_vector* xx);
void SetParams (int betaparm, int driver, double *params, vector<rec> traj);
void GetParams (double *p, int& betaparm, int& driver, vector<rec>& traj);
int CheckParams(const gsl_vector *xx, double *p);
double get_best_fit(const gsl_vector *xx, void *params);
void FindFrequencies (int driver, vector<rec> traj, vector<rec>& traj_inf, const gsl_vector *xx);
void FindFrequenciesFinal (int driver, vector<rec> traj, vector<rec>& traj_inf, const gsl_vector *xx, ofstream& out_file);
double GetLogLikelihood(int betaparm, vector<rec> traj, vector<rec> traj_inf, vector<double> fact_store);
void FindLogFact(vector<double>& fact_store,int N);
double BetaBinomCalc(int N, int r, float p, int betaparm, vector<double> fact_store);
