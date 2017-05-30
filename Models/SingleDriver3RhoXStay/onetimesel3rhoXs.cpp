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
	double i_final=0;
	double x_final=0;
	double e_final=0;
	double rho1_final=0;
	double rho2_final=0;
	double rho3_final=0;
	double rho4_final=0;
	double rx1_final=0;  //Point of change in recombination rate
	double rx2_final=0;  //Point of change in recombination rate
	double rx3_final=0;  //Point of change in recombination rate
	double xf_final=0;
	
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
	int driver=atoi(argv[5]);

	for (int r=1; r<=reps; r++) {

		//Set up initial parameters
		od.i = floor(gsl_rng_uniform(rgen)*p);
		od.x = 0.5+(0.5*gsl_rng_uniform(rgen));
		double m=min(pow(od.x,2),pow((1-od.x),2));
		od.e = (m*gsl_rng_uniform(rgen));
		od.rho1 = 1e-6;
		od.rho2 = 1e-6;
		od.rho3 = 1e-6;
		od.rho4 = 1e-6;
		od.rx1 = floor(0.33*gsl_rng_uniform(rgen)*p);
		od.rx2 = floor((0.33*p)+(0.33*gsl_rng_uniform(rgen)*p));
		od.rx3 = floor((0.66*p)+(0.33*gsl_rng_uniform(rgen)*p));
		od.xf =0.5+(0.5*gsl_rng_uniform(rgen));
		
		//General setup for optimisation routine
		
		//Optimisation processs
		for (int it=0;it<1000000;it++){
			
			//cout << "it " << it << "\n";
			size_t iter=0; //GSL iteration number
			
			//CheckLastLikelihood
			
			//Vector x contains parameters to be optimised
			int xx_size=10;
			gsl_vector *xx = gsl_vector_calloc(xx_size);
			SetXXVector(od,xx);
			//cout << "Done set xx\n";

			//Optimisation bit goes here
			
			int p_size = 4*(traj.size())+25;
			
			double *params;
			params=(double *)calloc(p_size,sizeof(double));
			
			SetParams(bparm,driver,params,traj);
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
			gsl_vector_set(ss,1,0.1);
			gsl_vector_set(ss,2,od.rho1/2.);
			gsl_vector_set(ss,3,od.rho2/2.);
			gsl_vector_set(ss,4,od.rho3/2.);
			gsl_vector_set(ss,5,od.rho4/2.);
			gsl_vector_set(ss,6,100);
			gsl_vector_set(ss,7,100);
			gsl_vector_set(ss,8,100);
			gsl_vector_set(ss,9,0.1);
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
				e_final=gsl_vector_get(s->x,1);
				rho1_final=gsl_vector_get(s->x,2);
				rho2_final=gsl_vector_get(s->x,3);
				rho3_final=gsl_vector_get(s->x,4);
				rho4_final=gsl_vector_get(s->x,5);
				rx1_final=gsl_vector_get(s->x,6);
				rx2_final=gsl_vector_get(s->x,7);
				rx3_final=gsl_vector_get(s->x,8);
				xf_final=gsl_vector_get(s->x,9);
				
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
	cout << "index " << floor(driver) << " locus " << traj[driver].pos << " x " << x_final << " rx1 " << rx1_final << " rx2 " << rx2_final << " rx3 " << rx3_final << " e " << e_final << " rho1 " << rho1_final << " rho2 " << rho2_final << " rho3 " << rho3_final << " rho4 " << rho4_final << " xf " << xf_final << " Log L " << best_fit << "\n";
	od.i = i_final;
	od.x = x_final;
	od.e = e_final;
	od.rho1 = rho1_final;
	od.rho2 = rho2_final;
	od.rho3 = rho3_final;
	od.rho4 = rho4_final;
	od.rx1 = rx1_final;
	od.rx2 = rx2_final;
	od.rx3 = rx3_final;
	od.xf = xf_final;
	cout << od.xf << "\n";
	int xx_size=10;
	gsl_vector *xx = gsl_vector_calloc(xx_size);
	SetXXVector(od,xx);
	//cout << "Allele frequencies\n";
	ofstream out_file;		
	out_file.open("both");
	FindFrequenciesFinal(driver,traj,traj_inf,xx,out_file);

	
	return 0;
}