//emission.cpp

//own headers...
#include "emission.h"

// Constructor
Emission::Emission(){
  rnd_emit = 1.0e-10;
  shrink   = 1.0;
  mode     = 0;// 1: binomial, 2: beta-binomial
  EmitProb_set = 0;
  dist_set = 0;
}


// real constructor
void Emission::set(int nsamples, vector<int>& nsites, int grid){
  nSamples  = nsamples;
  if (nSamples != (int) nsites.size()){
    cout<<"ERROR-1 in Emission::Emission()\n";
    exit(1);
  }
  nSites = new int [nSamples];
  for (int s=0; s<nSamples; s++){
    nSites[s] = nsites[s];
  }
  // set xgrid...
  gridSize = grid;
  dx = 1.0 / double(gridSize);
  xgrid    = new double[gridSize+1];
  for (int i=0; i<= gridSize; i++){
    xgrid[i] = double(i) * dx;
  }
  EmitProb  = new gsl_matrix * [nSamples];
  reads     = new unsigned int * [nSamples];
  depths    = new unsigned int * [nSamples];
  loci      = new unsigned int * [nSamples];
  for (int s=0; s<nSamples; s++){
    EmitProb[s] = gsl_matrix_alloc( nSites[s], gridSize + 1);
    reads[s]    = new unsigned int [nSites[s]];
    depths[s]   = new unsigned int [nSites[s]];
    loci[s]     = new unsigned int [nSites[s]];
  }
}

void Emission::set_dist(){
  if (loci == NULL){
    cout<<"ERROR-1 in Emission::set_dist()\n";
    exit(1);
  }
  dist = new double * [nSamples];
  total_loci=0;
  for (int s=0; s<nSamples; s++){
    total_loci += nSites[s];
    dist[s] = new double [nSites[s]];
    for (int l=1; l <nSites[s]; l++){
      dist[s][l] = fabs(double(loci[s][l] - loci[s][l-1]));
      if (dist[s][l] == 0.0){
	printf("ERROR: dist=0 in chr %i at locus %i\n", s+1, loci[s][l]);
	exit(1);
      }
    }
    dist[s][0] = 0.0; 
  }
  dist_set = 1;
  //printf("Data in %i sample(s) and with %i sites.\n", nSamples, total_loci);
}


Emission::~Emission(){
  for (int s=0; s<nSamples; s++){
    gsl_matrix_free(EmitProb[s]);
  }
  delete [] EmitProb;
  for (int s=0; s<nSamples; s++){
    delete [] dist[s];
  }
  delete [] dist;
  delete [] xgrid;
}

//emission probability as a function of total freq x
void Emission::set_EmitProb(){
  if (mode == 0){
    printf("ERROR-1 in Emission::set_EmitProb(): mode not set.\n");
    exit(1);
  }
  int s;
#pragma omp parallel for schedule( dynamic, 1) default(shared)
  for ( s=0; s<nSamples; s++){
    double x,f, p0,p;
    int n,N;
    for (int i=0; i < nSites[s]; i++){
      n = reads[s][i];
      N = depths[s][i];
      if (mode == 1){//binomial emission model
	for (int j=0; j<=gridSize; j++){
	  x = double(j)*dx;
	  f = (1.0-rnd_emit) * gsl_ran_binomial_pdf(n, x, N) + rnd_emit / double(N+1);
	  if (f<0.0 || f!= f){
	    printf("ERROR-2 in Emission::set_EmitProb(): %e\n", f);
	  }
	  gsl_matrix_set(EmitProb[s], i, j, f);
	}
      }
      else if (mode == 2){//beta-binomial emission model
	p0 = gsl_sf_lngamma(double(N+1)) - gsl_sf_lngamma(double(n+1)) - gsl_sf_lngamma(double(N-n+1)); 
	for (int j=0; j<=gridSize; j++){
	  x = double(j)*dx;
	  if (x==0.0){
	    p = (n==0) ? 1.0 : 0.0;
	  }
	  else if (x==1.0){
	    p = (n==N) ? 1.0 : 0.0;
	  }
	  else{
	    p  = p0 + gsl_sf_lnbeta(double(n) + shrink*x, double(N-n) + shrink*(1.0-x));
	    p -= gsl_sf_lnbeta( shrink*x,  shrink*(1.0-x));
	    if (p>0.0){
	      printf("ERROR in BetaBinomial::set_EmitProb(): p = %e\n", p);
	      exit(1);
	    }
	    p = exp(p);
	  }
	  f  = (1.0-rnd_emit) * p + rnd_emit / double(N+1);
	  gsl_matrix_set(EmitProb[s], i, j, f);
	}
      }
      else{
	printf("ERROR in Emission::set_EmitProb(): mode %i does not exist.\n", mode);
	exit(1);
      }
    }
  }
  EmitProb_set=1;
}
