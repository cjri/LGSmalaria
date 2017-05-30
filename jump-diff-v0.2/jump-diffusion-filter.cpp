/*
  jump-diffusion-filter.cpp
 
  Input: a file with four columns:
  chr locus reads depth

  Each chromosome is treated as an independent sample from the same process.
  The hidden state sequence is modeled as a jump-difusion trajectory in [0,1].
  Its value x at any point in time is the rate of a (Beta-) Binomial emission process that 
  generates the observed data.

  Output:
  The path of (local) maximum likelihood frequency.
  On request (--dist), the posterior distribution at each point (may be large files!).

  Options:
  --data        The input observed data file
  --pre         The prefix to put before all output files
  --grid [int]  The grid size for the distributions (partitions [0,1] into [grid] bins).
  --dist        Prints all the posterior distributions as well.
  --total       A fixed distribution for the total frequency
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <math.h>
#include <ctype.h> 
#include <string>
#include <map>
#include <vector>
#include <list>


// GSL headers...
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_blas.h"

//own headers...
#include "emission.h"
#include "jump-diffusion.h"

using namespace std;


struct cmdl_opts{
  const char * data_fn;
  const char * pre;
  const char * rec_fn;
  const char * total_fn;
  int grid, dist, nojump, mode,seed;
  double sigma, jump, rnd_emit, shrink;
};


//*** OWN FUNCTIONS ***
void get_opts( int argc, const char ** argv, cmdl_opts& opts);

void get_dims( const char * data_fn,
	       vector<int>& nSites,
	       vector<int>& chrs, 
	       vector<int>& idx_of
	       );

void get_data( const char * data_fn, Emission * myEmit,vector<int>& idx_of);
void include_rec(const char * rec_fn, JumpDiffusion * myJD, vector<int>& idx_of);
void get_total(gsl_matrix ** total, const char * total_fn, Emission * myEmit, vector<int>& idx_of);


double get_mean(gsl_vector * dist);
double get_var(gsl_vector * dist, double mean);

// *** MAIN START***
int main (int argc, const char * argv[]){
  cmdl_opts opts;
  get_opts( argc, argv, opts);
  srand(opts.seed);
  vector<int> chrs;
  vector<int> idx_of;
  vector<int> nSites;
  get_dims( opts.data_fn, nSites, chrs, idx_of);
  int nSamples = (int) chrs.size();
  int total_nLoci=0;
  for (int s=0; s<nSamples; s++){
    total_nLoci += nSites[s];
  }
  printf("Fitting jump-diffusion model to data at %i sites with a ", total_nLoci);
  if (opts.mode == 1){
    printf("binomial ");
  }
  else if (opts.mode == 2){
    printf("beta-binomial ");
  }
  printf("emission model:\n");
  // read in data: time points, read depths, reads
  Emission myEmit;
  myEmit.set(nSamples,nSites,opts.grid);
  myEmit.mode = opts.mode;  // mode 1: binomial, 2: beta-binomial
  get_data( opts.data_fn, &myEmit, idx_of);
  JumpDiffusion myJD(&myEmit);
  //
  // *** EXPLICIT TOTAL FREQ DISTRIBUTION: SCORE DATA ***
  // if explicit total frequency distribution is given, simply score the data against that
  if (opts.total_fn != NULL){
    // *** READ IN TOTAL DIST ***
    //alocate...
    gsl_matrix ** total = new gsl_matrix * [nSamples];
    for (int s=0; s<nSamples; s++){
      total[s] = gsl_matrix_calloc( nSites[s], myJD.gridSize+1);
    }
    //get total distribution...
    get_total( total, opts.total_fn, &myEmit, idx_of);
    myJD.wTotal   = 1;
    myJD.total    = total;
    myEmit.rnd_emit = opts.rnd_emit;
    myEmit.shrink   = opts.shrink;
    myEmit.set_EmitProb();
    myJD.do_Fwd();
    printf("LLH = %.5e, rnd = %.2e", myJD.total_llh, myEmit.rnd_emit);
    if (myEmit.mode == 2){
      printf(" shrink = %.2e\n",  myEmit.shrink);
    }
    else{
      cout<<endl;
    }
    return(0);//stop here!
  }
  //
  // *** LEARN PARAMETERS OF THE JUMP-DIFFUSION MODEL ***
  double * param = new double [4];
  int nvar=0;
  if (opts.nojump == 1 || opts.jump > 0.0){
    myJD.jump = (opts.nojump==1) ? 0.0 : opts.jump;
    param[0] = 0.0;
  }
  else{//jump probability
    param[0] = 1.0e-5;
    nvar++;
  }
  if (opts.sigma > 0.0){//diffusion constant
    myJD.sigma = opts.sigma;
    param[1] = 0.0;
  }
  else{
    param[1] = 1.0e-4;
    nvar++;
  }
  if (opts.rnd_emit > 0.0){//random error rate
    myEmit.rnd_emit = opts.rnd_emit;
    param[2] = 0.0;
  }
  else{
    param[2] = 0.01;
    nvar++;
  }
  if (opts.shrink > 0.0 || myEmit.mode == 1){//shrinkage
    myEmit.shrink = opts.shrink;
    param[3] = 0.0;
  }
  else{
    param[3] = 10.0;
    nvar++;
  }
  //
  if (opts.rec_fn != NULL) include_rec( opts.rec_fn, &myJD, idx_of);
  //
  double low_stop  = 1.0e-10;
  double high_stop = 1.0e4;
  while (nvar>0){
    find_parameter( param, &myJD, low_stop, high_stop);
    if (nvar==1) break;
    // check if needs to be redone...
    int redo=0;
    for (int i=0;i<4;i++){
      if (param[i]>0.0 && (param[i] < low_stop || param[i] > high_stop)){
	double fix = (param[i] < low_stop) ? low_stop : high_stop;
	if (i==0) myJD.jump       = fix;
	if (i==1) myJD.sigma      = fix;
	if (i==2) myEmit.rnd_emit = fix;
	if (i==3){
	  myEmit.shrink   = fix;
	  myEmit.mode     = 1;//switch to binomial
	}
	param[i] = 0.0;
	nvar--;
	redo=1;
      }
    }
    if (redo==0) break;
  }
  //
  if ( param[0] > 0.0) myJD.jump       = param[0];
  if ( param[1] > 0.0) myJD.sigma      = param[1];
  if ( param[2] > 0.0) myEmit.rnd_emit = param[2];
  if ( param[3] > 0.0) myEmit.shrink   = param[3];
  // do forward-backward for jump-diffusion...
  myJD.do_Fwd();
  myJD.do_Bwd();
  printf("Jump-Diffusion process: LLH = %.5e, --jump %.2e --sigma %.2e --rnd %.2e",
	 myJD.total_llh, myJD.jump,  myJD.sigma, myEmit.rnd_emit);
  if (myEmit.mode == 2){
    printf(" --shrink = %.2e\n",  myEmit.shrink);
  }
  else{
    cout<<endl;
  }
  // *** DONE LEARNING JUMP-DIFFUSION PARAMETERS ***
  //
  // print total mean frequencies and standard deviations...
  char buff[128];
  sprintf(buff,"%s.total.txt", opts.pre);
  FILE * total_fp = fopen(buff,"w");
  fprintf(total_fp, "#sample site mean std-dev jump-prob posterior\n");
  for (int s=0; s < myJD.nSamples; s++){
    for (int l=0; l < myJD.nSites[s]; l++){
      gsl_vector_view post = gsl_matrix_row(myJD.gamma[s],l);
      double mean = get_mean(&post.vector);
      double var  = get_var(&post.vector,mean);
      fprintf(total_fp, "%i %6i %.6e %.6e %.6e", chrs[s], myJD.loci[s][l], mean, sqrt(var), exp(myJD.pjump[s][l]));
      if(opts.dist==1){// full posterior distribution? LARGE!
	for (int i=0; i <= myJD.gridSize; i++){
	  double p = gsl_matrix_get( myJD.gamma[s], l, i);
	  //p = (p < 1.0e-10) ? 0 : p;
	  fprintf(total_fp, " %.2e", p); 
	}
      }
      fprintf(total_fp,"\n");
    }
  }
  fclose(total_fp);
  return (0);
}
// *** MAIN END ***


// get command line arguments...
void get_opts( int argc, const char ** argv, cmdl_opts& opts){
  int opt_idx = 1;
  string opt_switch;
  opts.data_fn  = NULL;
  opts.rec_fn   = NULL;
  opts.total_fn   = NULL;
  opts.pre    = "./out";
  opts.grid   = 100;
  opts.dist   = 1;
  opts.nojump = 0;
  opts.sigma    = 0.0;
  opts.jump     = 0.0;
  opts.rnd_emit = 0.0;
  opts.shrink   = 0.0;
  opts.mode = 1;
  opts.seed = (int) time(NULL);
  while ( opt_idx < argc && (argv[opt_idx][0] == '-')){
    opt_switch = argv[opt_idx];
    if ( opt_switch.compare("--data") == 0){
      opt_idx++;
      opts.data_fn = argv[opt_idx];
    }
    else if ( opt_switch.compare("--pre") == 0){
      opt_idx++;
      opts.pre = argv[opt_idx];
    }
    else if ( opt_switch.compare("--rec") == 0){
      opt_idx++;
      opts.rec_fn = argv[opt_idx];
    }
    else if ( opt_switch.compare("--total") == 0){
      opt_idx++;
      opts.total_fn = argv[opt_idx];
    }
    else if ( opt_switch.compare("--grid") == 0){
      opt_idx++;
      opts.grid = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--seed") == 0){
      opt_idx++;
      opts.seed = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--mode") == 0){
      opt_idx++;
      opts.mode = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--sigma") == 0){
      opt_idx++;
      opts.sigma = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--rnd") == 0){
      opt_idx++;
      opts.rnd_emit = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--jump") == 0){
      opt_idx++;
      opts.jump = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--shrink") == 0){
      opt_idx++;
      opts.shrink = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--dist") == 0){
      opts.dist = 1;
    }
    else if ( opt_switch.compare("--nojump") == 0){
      opts.nojump = 1;
    }
    else {
      cout << "Usage: jump-diffusion-filter --data input.txt --pre ./out-dir/out (--grid [int] --sigma [double] ";
      cout<<"--jump [double] --rnd [double] --nojump --rec [file] --dist [0,1]  --total [file] "<<endl;
      exit(1);
    }
    opt_switch.clear();
    opt_idx++;
  }
}


void get_dims( const char * data_fn, 
	       vector<int>& nSites,
	       vector<int>& chrs,
	       vector<int>& idx_of
	       ){
  ifstream data_ifs;
  string line;
  stringstream line_ss;
  data_ifs.open( data_fn, ios::in);
  if (data_ifs.fail()){
    printf("ERROR in get_data(): file %s cannot be opened.\n", data_fn);
    exit(1);
  }
  nSites.clear();
  int ct=0;
  int chr,old=-1;
  int maxChr=0;
  while( data_ifs.good()){
    line.clear();
    getline( data_ifs, line);
    if (line.empty()) break;
    line_ss.clear();
    line_ss.str(line);
    line_ss >> chr;
    if (chr != old ){//new chromosome encounter
      if (ct>0) nSites.push_back(ct);
      chrs.push_back(chr);
      if (chr > maxChr) maxChr = chr;
      ct=0;
    }
    old=chr;
    ct++;
  }
  nSites.push_back(ct);
  idx_of.clear();
  for (int c=0; c <= maxChr; c++){
    idx_of.push_back(-1);
  }
  for (int c=0; c < (int) chrs.size(); c++){
    idx_of[chrs[c]] = c;
  }
  data_ifs.close();
}


// read in data: expects columns to be "chr location depth reads"
void get_data( const char * data_fn, 
	       Emission * myEmit,
	       vector<int>& idx_of
	      ){
  //
  ifstream data_ifs;
  string line;
  stringstream line_ss;
  data_ifs.open( data_fn, ios::in);
  if (data_ifs.fail()){
    printf("ERROR in get_data(): file %s cannot be opened.\n", data_fn);
    exit(1);
  }
  int nchr = myEmit->nSamples;
  int ct=0,l;
  int chr,old=-1;
  int d,r;
  //now collect all data...
  while( data_ifs.good()){
    line.clear();
    getline( data_ifs, line);
    if (line.empty()) break;
    line_ss.clear();
    line_ss.str(line);
    line_ss >> chr >> l >> r >> d;
    if (chr != old){
      ct  = 0;
      old = chr;
    }
    if( idx_of[chr] < 0 || idx_of[chr] >= nchr){
      cout<<endl<<line<<endl;
      printf("%i %i %i\n", chr, idx_of[chr], nchr);
      exit(1);
    }

    if (d == 0 && r > 0){
      printf("ERROR: depth = 0 but read = %i in chr %i locus %i\n", r, chr, l);
      cout<<flush;
      exit(1);
    }
    myEmit->loci[ idx_of[chr] ][ct] = l;
    myEmit->reads[idx_of[chr]][ct]  = r;
    myEmit->depths[idx_of[chr]][ct] = d;
    ct++;
  }  
  data_ifs.close();
  myEmit->set_dist();
}


void include_rec(const char * rec_fn, JumpDiffusion * myJD, vector<int>& idx_of){
  ifstream rec_ifs;
  string line;
  stringstream line_ss;
  rec_ifs.open( rec_fn, ios::in);
  if (rec_ifs.fail()){
    printf("ERROR in include_rec(): file %s cannot be opened.\n", rec_fn);
    exit(1);
  }
  int chr, ochr=0, l=0, locus;//, ol=-1;
  double rec;
  while( rec_ifs.good()){
    line.clear();
    getline( rec_ifs, line);
    if (line.empty()) break;
    line_ss.clear();
    line_ss.str(line);
    line_ss >> chr >> locus >> rec;
    //if (locus == ol) break;
    if (chr != ochr){
      l=0;
      ochr=chr;
    }
    else{
      l++;
      //ol = locus;
    }
    if (l > myJD->nSites[ idx_of[chr] ] - 1){
      printf("ERROR: chr %i only has %i sites, but I am here at %i!\n", 
	     chr, myJD->nSites[ idx_of[chr] ], l+1);
      exit(1);
    }
    if ( fabs( myJD->loci[ idx_of[chr] ][l] - locus) > 1){
      printf("ERROR: locus %i at coordinate %i in chr %i is not the same as %i\n", 
	     l+1, myJD->loci[ idx_of[chr] ][l], chr, locus);
      exit(1);
    }
    // SCALE HERE
    if (rec > 0.0){
      myJD->dist[ idx_of[chr] ][l+1] *= rec;
      //printf("%.2e\n", myHBJD->dist[chr-1][l+1]);
    }
  }
  rec_ifs.close();
}



void get_total(gsl_matrix ** total, const char * total_fn, Emission * myEmit, vector<int>& idx_of){
  printf("Using data in %s to score the data (site-by-site)...\n", total_fn);
  ifstream ifs;
  string line;
  stringstream line_ss;
  ifs.open( total_fn, ios::in);
  if (ifs.fail()){
    printf("ERROR in get_bulk(): file %s cannot be opened.\n", total_fn);
    exit(1);
  }
  int chr = 0, old_chr = -1, l=0, locus;
  double mn, sd, jp;
  while( ifs.good() ){
    line.clear();
    getline( ifs, line);
    if (line.empty()) break;
    if (line[0] == '#') continue;
    line_ss.clear();
    line_ss.str(line);
    line_ss >> chr >> locus >> mn >> sd >> jp;
    l = (chr == old_chr) ? l+1 : 0;
    if (chr != old_chr && chr > (int) idx_of.size() - 1){
      printf("ERROR 1a in get_bulk()\n");
      cout<<line<<endl;
      exit(1);
    }
    if (chr != old_chr && idx_of[chr] == -1){
      printf("ERROR 1 in get_bulk()\n");
      printf("chr=%i, idx=%i\n", chr, idx_of[chr]);
      exit(1);
    }
    // exit(0);
    old_chr = chr;
    if( (int) myEmit->loci[ idx_of[chr] ][l] != locus){
      printf("ERROR 2 in get_bulk()\n");
      cout<<line<<endl;
      printf( "%i, %i, %i vs %i\n", idx_of[chr], l, myEmit->loci[ idx_of[chr] ][l], locus);
      exit(1);
    }
    double p;
    for (int i=0; i < myEmit->gridSize; i++){
      if (line_ss.good() != true){
	printf("ERROR 3 in get_bulk()\n");
	exit(1);
      }
      line_ss >> p;
      gsl_matrix_set( total[ idx_of[chr] ], l, i, p);
    }

  }
  ifs.close();
}





double get_mean(gsl_vector * dist){
  double mean=0.0,P1,P2;
  int n = (int) dist->size;
  double dx = 1.0 / double(n-1);
  for (int i=0; i <= n-2; i++){
    P1 = gsl_vector_get(dist,i);
    P2 = gsl_vector_get(dist,i+1);
    mean += ( 0.5*(P1+P2)*double(i) + (P1+2.0*P2) / 6.0 ) * pow(dx,2);
  }
  return(mean);
}


double get_var(gsl_vector * dist, double mean){
  double var=0.0, P1, P2,dev;
  int n = (int) dist->size;
  double dx = 1.0 / double(n-1);
  for (int i=0; i<n-1; i++){
    P1 = gsl_vector_get(dist,i);
    P2 = gsl_vector_get(dist,i+1);
    dev = double(i)*dx - mean;
    var += (1.0/12.0)*dx*((P1+3.0*P2)*dx*dx + 4.0*dx*(P1+2.0*P2)*dev + 6.0*(P1+P2)*pow(dev,2));
  }
  return(var);
}
