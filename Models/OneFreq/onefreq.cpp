#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <deque>
#include <map>
using namespace std;
#include "shared.h"

int main(int argc, char *argv[]){

	vector<rec> traj;
	vector<rec> traj_inf;
	double x_final=0;
	
	string input(argv[1]);
	int seed=atoi(argv[2]);

	//Set up random number generator
	srand((unsigned int)(seed)); 
	gsl_rng *rgen = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (rgen, seed);
		
	ImportData(input,traj);
	
	//return 0;
	
	//Number of data points
	int p=traj.size();
	//cout << p << "\n";
	
	//Parameters are l (locus), x (cross frequency), e (vertical movement), rho (recombination), xf (driver frequency)

	opt_data od;
	double best_fit=-1e10;  
	double last_fit=-1e10;

	int reps=atoi(argv[3]);
	int bparm=atoi(argv[4]);

	for (int r=1; r<=reps; r++) {

		//Set up initial parameters
		od.x = 0.5+(0.5*gsl_rng_uniform(rgen));
		
		//General setup for optimisation routine
		
		//Optimisation processs
		for (int it=0;it<1000000;it++){
			
			//cout << "it " << it << "\n";
			size_t iter=0; //GSL iteration number
			
			//CheckLastLikelihood
			
			//Vector x contains parameters to be optimised.  Parameters are locus and two frequencies
			int xx_size=1;
			gsl_vector *xx = gsl_vector_calloc(xx_size);
			SetXXVector(od,xx);
			//cout << "Done set xx\n";

			//Optimisation bit goes here
			
			int p_size = 4*(traj.size())+25;
			
			double *params;
			params=(double *)calloc(p_size,sizeof(double));
			
			SetParams(bparm,params,traj);
			//cout << "Done set params\n";

			//Define the optimisation function
			gsl_multimin_function my_func;
			my_func.n=xx_size;
			my_func.f=&get_best_fit;
			my_func.params=params;
			
			gsl_multimin_fminimizer *s;
			const gsl_multimin_fminimizer_type *T=gsl_multimin_fminimizer_nmsimplex;
			s=gsl_multimin_fminimizer_alloc(T,xx_size);
			//cout << "Done define\n";

			//Starting change magnitudes
			gsl_vector* ss=gsl_vector_alloc(xx_size);
			gsl_vector_set(ss,0,0.1);
			//cout << "Done set ss\n";

			gsl_multimin_fminimizer_set (s,&my_func,xx,ss);
			//cout << "Done set mmin\n";
			int status;
			double size;
			do {
				iter++;
				status=gsl_multimin_fminimizer_iterate(s);
				if (status) {break;}
				
				size=gsl_multimin_fminimizer_size(s);
				status = gsl_multimin_test_size(size, 1e-4);
				
				if (status == GSL_SUCCESS) {
				}
				
			} while (status==GSL_CONTINUE && iter<1000);
			
			last_fit = s->fval;	
			last_fit=-last_fit;
			cout << "Last fit score = " << last_fit <<"\n";

			if (last_fit>best_fit) {
				best_fit=last_fit;
				x_final=gsl_vector_get(s->x,0);
				
			} else {
				break;
			}
			
			gsl_vector_free(xx);
			gsl_vector_free(ss);
			gsl_multimin_fminimizer_free(s);
			
		
		}
	}
	
	
	//Print values and frequencies
	cout << "Optimised values\n";
	cout << " x " << x_final << " Log L " << best_fit << "\n";
	od.x = x_final;
	int xx_size=1;
	gsl_vector *xx = gsl_vector_calloc(xx_size);
	SetXXVector(od,xx);
	//cout << "Allele frequencies\n";
	ofstream out_file;		
	out_file.open("both");
	FindFrequenciesFinal(traj,traj_inf,xx,out_file);

	
	return 0;
}