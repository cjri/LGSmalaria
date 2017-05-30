//jump-diffusion.cpp

//own headers...
#include "emission.h"
#include "jump-diffusion.h"

// Constructor
JumpDiffusion::JumpDiffusion( Emission * emit){
  myEmit   = emit;
  nSamples = myEmit->nSamples;
  nSites   = myEmit->nSites;
  dist = myEmit->dist;
  loci = myEmit->loci;
  jump      = 1.0e-4;
  sigma     = 1.0e-3;
  wTotal    = 0;
  pstay     = new double * [nSamples];
  pjump     = new double * [nSamples];
  pnojump   = new double * [nSamples];
  // matrices...
  alpha  = new gsl_matrix * [nSamples];
  gamma  = new gsl_matrix * [nSamples];  
  gridSize = myEmit->gridSize;
  dx = 1.0/double(gridSize);
  xgrid    = myEmit->xgrid;
  proposal = gsl_vector_calloc(gridSize+1);
  gsl_vector_set_all( proposal, 1.0);// uniform proposal distribution
  for (int s=0; s<nSamples; s++){
    pstay[s]   = new double [nSites[s]];
    pjump[s]   = new double [nSites[s]];
    pnojump[s] = new double [nSites[s]];
    alpha[s] = gsl_matrix_alloc(nSites[s],gridSize+1);
    gamma[s] = gsl_matrix_alloc(nSites[s],gridSize+1);
  }
  FwdBwd_done = 0;
}

JumpDiffusion::~JumpDiffusion(){
  for (int s=0; s<nSamples; s++){
    delete [] pstay[s];
    delete [] pjump[s];
    delete [] pnojump[s];
    gsl_matrix_free(alpha[s]);
    gsl_matrix_free(gamma[s]);
  }
  delete [] alpha;
  delete [] pstay;
  delete [] pjump;
  delete [] pnojump;
  delete [] gamma;
}


// Diffusion propagator
void JumpDiffusion::set_DiffProp(gsl_matrix * propagator, double sd){
  int range = 3 * ceil(sd / dx);// this means, only fluctuations up to 3 sigma are possible
  if (2*range > gridSize+1){
    sd = dx * double(gridSize) / 6.0;
    range = 3 * ceil(sd / dx);
  }
  gsl_vector * gauss = gsl_vector_alloc(2*range+1);
  for (int i=0; i<2*range+1; i++){
    gsl_vector_set(gauss,i, gsl_ran_gaussian_pdf( double(i-range)*dx, sd));
  }
  gsl_matrix_set_zero(propagator);
  double val = 0, norm=0;
  for (int i=0; i<= gridSize; i++){
    norm = 0.0;
    for (int j=0; j<= gridSize; j++){
      if ( abs(i-j) <= range){
	val = gsl_vector_get( gauss, j-i+range);
	if( i+1 <= range && j+1 <= range - i){
	  val += gsl_vector_get( gauss, range - i - j - 1);
	}
	if( i+1 >= gridSize+1 - range && i+j >= 2*(gridSize+1) - range - 1){
	  val += gsl_vector_get( gauss, 2*(gridSize+1) + range - i - j - 1);
	}
	gsl_matrix_set(propagator,i,j,val);
	norm += (j==0 || j==gridSize) ? 0.5*val : val;
      }
    }
    gsl_vector_view row = gsl_matrix_row(propagator,i);
    //double norm = gsl_blas_dasum(&row.vector);
    //norm -= 0.5*(gsl_matrix_get(propagator,i,0) + gsl_matrix_get(propagator,i,grid_size));
    norm *= dx;
    gsl_vector_scale(&row.vector,1.0/norm);
    if (i==0 || i==gridSize) gsl_vector_scale( &row.vector, 0.5);//trapezoidal integration rule
  }
  gsl_vector_free(gauss);
}



void JumpDiffusion::do_Fwd(){
  total_llh   = 0.0;
  if (myEmit->EmitProb_set == 0){
    myEmit->set_EmitProb();
  }
  int s;
  // SAMPLES:
#pragma omp parallel for schedule( dynamic, 1) default(shared)
  for (s=0; s<nSamples; s++){
    // first compute the staying probabilities..
    if (jump > 0.0){
      for (int l=1; l< nSites[s]; l++){
	//pstay[s][l] = pow( 1.0-stiffness, int(dist[s][l]));
	pstay[s][l] = exp( - jump * dist[s][l]);
      }
    }
    else{
      for (int l=1; l<nSites[s]; l++){
	pstay[s][l] = 1.0;
      }
    }
    pstay[s][0]   = 1.0;
    // prepare fwd...
    gsl_vector * eprob = gsl_vector_alloc(gridSize+1);
    gsl_vector * prior = gsl_vector_alloc(gridSize+1);
    gsl_vector * post  = gsl_vector_alloc(gridSize+1);
    gsl_vector * mem    = gsl_vector_alloc(gridSize+1);
    gsl_matrix * DiffProp = gsl_matrix_calloc(gridSize+1,gridSize+1);
    double norm;
    // Forward Pass
    for (int l=0; l<nSites[s]; l++){
      gsl_vector_memcpy( prior, proposal);
      if (wTotal == 0 && l>0){
	// get diffusion propagator...
	JumpDiffusion::set_DiffProp( DiffProp, sigma * sqrt(dist[s][l]));
	//diffusion scaled by time difference
	gsl_blas_dgemv(CblasTrans, dx*pstay[s][l], DiffProp, post, 0.0, mem);
	gsl_vector_memcpy( post, mem);
	gsl_vector_scale( prior, 1.0 - pstay[s][l]);
	gsl_vector_add( prior, post);
      }
      else if (wTotal == 1){// if explicit distribution is given
	gsl_matrix_get_row( prior, total[s], l);
      }
      //emission probability
      gsl_matrix_get_row( eprob, myEmit->EmitProb[s], l);
      gsl_vector_mul( prior, eprob);// at this time it is the posterior!
      norm  = gsl_blas_dasum(prior);
      norm -= 0.5*(prior->data[0] + prior->data[gridSize]);
      norm *= dx;
      if (norm <=0 || norm != norm){
	printf ("ERROR-1 in  JumpDiffusion::do_Fwd()\n");
	exit(1);
      }
      gsl_vector_scale(prior, 1.0 / norm);
      gsl_vector_memcpy(post,prior);
#pragma omp critical
      {
	total_llh += log(norm);// get part of the total log-likelihood
	gsl_matrix_set_row(alpha[s], l, post);// set forward variable
      }      
    }
    // clean up...
    gsl_vector_free(eprob);
    gsl_vector_free(prior);
    gsl_vector_free(post);
    gsl_vector_free(mem);
    gsl_matrix_free(DiffProp);
  }
  // FWD DONE
}


// Backward Pass (don't really need to save beta, maybe not even gamma)
// but we need to collect all possible observables here!
void JumpDiffusion::do_Bwd(){
  int s;
  // SAMPLES:
#pragma omp parallel for schedule( dynamic,1) default(shared)
  for (s=0; s<nSamples; s++){
    // prepare bwd...
    gsl_vector * eprob = gsl_vector_alloc(gridSize+1);
    gsl_vector * prior = gsl_vector_alloc(gridSize+1);
    gsl_vector * post  = gsl_vector_alloc(gridSize+1);
    gsl_vector * beta  = gsl_vector_alloc(gridSize+1);
    gsl_vector * last_beta  = gsl_vector_alloc(gridSize+1);
    gsl_vector * last_eprob = gsl_vector_alloc(gridSize+1);
    gsl_vector * mem    = gsl_vector_alloc(gridSize+1);
    gsl_vector * mem2    = gsl_vector_alloc(gridSize+1);
    gsl_matrix * DiffProp = gsl_matrix_calloc(gridSize+1,gridSize+1);
    double x,y;
    for (int l = nSites[s]-1; l>=0; l--){
      gsl_vector_memcpy(prior,proposal);
      if (l<nSites[s]-1){
	JumpDiffusion::set_DiffProp( DiffProp, sigma*sqrt(dist[s][l+1]));
	//diffusion scaled by time difference
	gsl_blas_dgemv(CblasTrans, dx*pstay[s][l+1], DiffProp, post, 0.0, mem);
	gsl_vector_memcpy( post, mem);
	gsl_vector_scale( prior, 1.0 - pstay[s][l+1]);
	gsl_vector_add( prior, post);
      }
      gsl_vector_memcpy( beta, prior);
      // get gamma, i.e. the total posterior probability vector
      gsl_vector_view alph = gsl_matrix_row(alpha[s],l);
      gsl_vector_memcpy( post, beta);
      gsl_vector_mul( post, &alph.vector);
      double norm = gsl_blas_dasum(post);
      norm -= 0.5*(post->data[0] + post->data[gridSize]);
      norm *= dx;
      //
      gsl_vector_scale( post, 1.0 / norm);
#pragma omp critical
      {
	gsl_matrix_set_row( gamma[s], l, post);// the posterior on-site sojourn probability
      }
      // get the emission probabilty for this site...
      gsl_matrix_get_row( eprob, myEmit->EmitProb[s], l);
      // get posterior update for the next step...
      gsl_vector_mul( prior, eprob);// now it is the posterior!
      norm = gsl_blas_dasum(prior);
      norm -= 0.5*(prior->data[0] + prior->data[gridSize]);
      norm *= dx;
      //
      gsl_vector_scale(prior, 1.0 / norm);
      gsl_vector_memcpy( post, prior);
      // posterior jump-probability...
      if (l<nSites[s]-1){
        gsl_vector_mul(last_beta, last_eprob);
        gsl_vector_memcpy(mem2,last_beta);
        gsl_blas_ddot( proposal, last_beta, &x);
        x *= (1.0 - pstay[s][l+1]);
        gsl_blas_dgemv( CblasNoTrans, pstay[s][l+1], DiffProp, mem2, 0.0, mem);
        gsl_blas_ddot( mem, &alph.vector, &y);
        // log(pjump) and log(1.0-pjump) for the transition l->l+1
        if (jump > 0.0){
          pjump[s][l+1]   = log(x) - log(x+y);
          pnojump[s][l+1] = log(y) - log(x+y);
        }
        else{
          pjump[s][l+1]   = -1.0e6;
          pnojump[s][l+1] = 0.0;
        }
      }
      gsl_vector_memcpy(last_beta,beta);
      gsl_vector_memcpy(last_eprob,eprob);
    }
    //clean up...
    gsl_vector_free(eprob);
    gsl_vector_free(prior);
    gsl_vector_free(post);
    gsl_vector_free(beta);
    gsl_vector_free(mem);
    gsl_vector_free(mem2);
    gsl_vector_free(last_beta);
    gsl_vector_free(last_eprob);
    gsl_matrix_free(DiffProp);
  }
  FwdBwd_done = 1;
}


// The  numerical minimization routine to get the clone frequencies...
double find_parameter(double * param, JumpDiffusion * myJD, double low_stop, double high_stop){
  // Here starts the minimizing step...
  int iter = 0, max_iter = 1.0e3; // max no. iterations
  int status;
  double val;
  gsl_multimin_function my_func;
  vector<int> to_opt;
  int npar = 4;
  string * par_names = new string [npar];
  par_names[0] = "jump";
  par_names[1] = "sigma";
  par_names[2] = "rnd";
  par_names[3] = "shrink";
  for (int i=0; i<npar; i++){
    if(param[i] > 0.0) to_opt.push_back(i);
  }
  int nvar = (int) to_opt.size();
  if (nvar==0){
    printf("ERROR: nvar = 0.\n");
    exit(1);
  }
  else{
    //printf("nvar = %i.\n", nvar);
  }
  gsl_vector * x  = gsl_vector_alloc(nvar);;
  gsl_vector * dx = gsl_vector_alloc(nvar);;
  my_func.f = &Qall;
  my_func.n  = nvar;
  //
  for (int i=0; i<nvar; i++){
    x->data[i]  = log(param[to_opt[i]]);
    dx->data[i] = max(0.1,0.1*fabs(x->data[i]));
  }
  fpar myfpar;
  myfpar.myJD    = myJD;
  myfpar.to_opt  = to_opt;
  my_func.params = static_cast<void*>(&myfpar);
  // Define type of minimization procedure...
  const gsl_multimin_fminimizer_type * T = gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer * s = gsl_multimin_fminimizer_alloc( T, nvar);
  // Set and initialize the minimizer with LOG-TRANSFORMED DATA...
  gsl_multimin_fminimizer_set( s, &my_func, x, dx);
  // Now iterate to find the minimum...
  do{
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);     
    if (status) break;
    status = gsl_multimin_test_size( gsl_multimin_fminimizer_size(s), 1.0e-2);
    printf("\rStep %3i: llh = %+.6e ", iter, - gsl_multimin_fminimizer_minimum(s));
    for (int i=0; i<nvar; i++){
      printf("%s = %.5e, ", par_names[to_opt[i]].c_str(), exp((gsl_multimin_fminimizer_x(s))->data[i]));
      val = exp((gsl_multimin_fminimizer_x(s))->data[i]);
      if ( val < low_stop || val > high_stop) status = 2;
    }
    cout<<flush;
  } 
  while (status == GSL_CONTINUE && iter < max_iter);
  cout<<endl;
  //
  double f = gsl_multimin_fminimizer_minimum(s);
  // copy the final result back into the argument, after EXPONENTIATION
  for (int i=0; i<nvar; i++) param[to_opt[i]] = exp((gsl_multimin_fminimizer_x(s))->data[i]);
  // cleanup...
  gsl_multimin_fminimizer_free(s);
  gsl_vector_free(x);
  gsl_vector_free(dx);
  return(-f);
}


double Qall( const gsl_vector * x, void * p){
  fpar * myfpar = static_cast<fpar*> (p);
  for (int i=0; i<(int) x->size; i++){
    // jump probability
    if ((myfpar->to_opt)[i] == 0){
      if (x->data[i] > 0.0){
	return(1.0e10);
      }
      else{
	(myfpar->myJD)->jump = exp(x->data[i]);
      }
    }
    // diffusion constant
    else if ((myfpar->to_opt)[i] == 1){
      (myfpar->myJD)->sigma = exp(x->data[i]);
    }
    // random emissions
    else if ((myfpar->to_opt)[i] == 2){
      if (x->data[i] > 0.0){
	return(1.0e10);
      }
      else{
	double rnd = exp(x->data[i]);
	 (myfpar->myJD)->myEmit->rnd_emit = rnd;
	 (myfpar->myJD)->myEmit->set_EmitProb();
      }
    }
    // shrinkage
    else if ((myfpar->to_opt)[i] == 3){
      double shrink = exp(x->data[i]);
      (myfpar->myJD)->myEmit->shrink = shrink;
      (myfpar->myJD)->myEmit->set_EmitProb();
    }
  }
  // DO FWD TO GET LLH
  (myfpar->myJD)->do_Fwd();
  return(-(myfpar->myJD)->total_llh);
}
