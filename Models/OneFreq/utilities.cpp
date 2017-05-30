using namespace std;
#include "shared.h"

void ImportData (string input, vector<rec>& traj) {
	ifstream pass_file;		
	pass_file.open(input.c_str());

	int i=0;
	do {
		rec new_traj;
		pass_file >> new_traj.chr;
		pass_file >> new_traj.pos;
		pass_file >> new_traj.obs;
		pass_file >> new_traj.N;
		new_traj.index=i;
		i++;
		new_traj.q=(new_traj.obs+0.)/(new_traj.N+0.);
		//cout << i << " " << new_traj.q << "\n";
		traj.push_back(new_traj);
		
	} while (!pass_file.eof());
}

void SetXXVector(opt_data od, gsl_vector* xx) {
	gsl_vector_set(xx,0,od.x);	//Frequency
}

void SetParams (int betaparm, double *params, vector<rec> traj) {
	params[0]=traj.size();
	params[1]=betaparm;
	params[2]=traj.size();
	params[3]=traj.size();
	for (unsigned int i=1;i<traj.size();i++) {
		params[(4*i)]=traj[i].chr;
		params[(4*i)+1]=traj[i].pos;
		params[(4*i)+2]=traj[i].obs;
		params[(4*i)+3]=traj[i].N;
	}
}

void GetParams (double *p, int& betaparm, vector<rec>& traj) {
	unsigned int ts=p[0];
	betaparm=p[1];
	for (unsigned int i=1;i<ts;i++) {
		rec new_rec;
		new_rec.chr=p[(4*i)];
		new_rec.pos=p[(4*i)+1];
		new_rec.obs=p[(4*i)+2];
		new_rec.N=p[(4*i)+3];
		traj.push_back(new_rec);
	}
}

double get_best_fit(const gsl_vector *xx, void *params) {
	double L=0;
	double *p = (double *)params;
	int betaparm=0;
	int chk=CheckParams(xx,p);
	if (chk==0) {
		return 1e10;
	}
	
	vector<rec> traj;
	vector<rec> traj_inf;
	vector<double> fact_store;

	GetParams(p,betaparm,traj);
	//cout << "i "  << gsl_vector_get(xx,0) << " x " << gsl_vector_get(xx,1) << " e " << gsl_vector_get(xx,2) << " rho " << gsl_vector_get(xx,3) << " xf " << gsl_vector_get(xx,4) << "\n";
	FindFrequencies(traj,traj_inf,xx);
	//cout << "Done FindFrequencies\n";
	L=GetLogLikelihood(betaparm,traj,traj_inf,fact_store);
	//cout << "L= " << L << "\n";
	L=-L;
	return L;
}

int CheckParams(const gsl_vector *xx, double *p) {
	int chk=1;
	double x=gsl_vector_get(xx,0);
	if (x<=0.0 || x>=1.0) {
		chk=0;
	}
	return chk;
}

void FindFrequencies (vector<rec> traj, vector<rec>& traj_inf, const gsl_vector *xx) {
	//Set up traj_inf
	traj_inf=traj;

	//Optimisation parameters
	double x=gsl_vector_get(xx,0);

	//Frequency at each of the loci
	for (unsigned int j=0;j<traj.size();j++) {
		traj_inf[j].q=x;
	}
}

void FindFrequenciesFinal (vector<rec> traj, vector<rec>& traj_inf, const gsl_vector *xx, ofstream& out_file) {
	//Set up traj_inf
	traj_inf=traj;
	//cout << "i "  << gsl_vector_get(xx,0) << " x " << gsl_vector_get(xx,1) << " e " << gsl_vector_get(xx,2) << " rho " << gsl_vector_get(xx,3) << " xf " << gsl_vector_get(xx,4) << "\n";

	//Optimisation parameters
	double x=gsl_vector_get(xx,0);
	
	//Frequency at each of the loci
	for (unsigned int j=0;j<traj.size();j++) {
		traj_inf[j].q=x;
		out_file << j << " " << traj[j].pos << " " << traj[j].q << " " << traj_inf[j].q << "\n";
	}
}

double GetLogLikelihood(int betaparm, vector<rec> traj, vector<rec> traj_inf, vector<double> fact_store) {
	double L=0;
	for (unsigned int i=0;i<traj.size();i++) {
		//cout << "i= " << i << " L= " << L << "\n";
		double q=traj_inf[i].q;
		double n=traj[i].N;
		double r=traj[i].obs;
		double p_next=(1.0)*BetaBinomCalc(n,r,q,betaparm,fact_store);
		L=L+log(p_next);
	}
	return(L);
}

double BetaBinomCalc(int N, int r, float p, int betaparm, vector<double> fact_store) {
	double alpha=betaparm*p;
	double beta=betaparm*(1-p);
	double bin=gsl_sf_lngamma(N+1)+gsl_sf_lngamma(r+alpha)+gsl_sf_lngamma(N-r+beta)+gsl_sf_lngamma(alpha+beta)-gsl_sf_lngamma(r+1)-gsl_sf_lngamma(N-r+1)-gsl_sf_lngamma(N+alpha+beta)-gsl_sf_lngamma(alpha)-gsl_sf_lngamma(beta);
	bin=exp(bin);
	if (bin<1e-300) {
		bin=1e-300;
	}
	//cout<<"Binomial N r p0 = "<<N<<" "<<r<<" "<<p<<" "<<bin<<"\n";
	return(bin);
}


void FindLogFact(vector<double>& fact_store,int N){
	double logN=0;
	for (int i=1;i<=N;i++) {
		logN=logN+log(i);
		fact_store.push_back(logN);
		//cout << "fact_store "<<i<<" "<<gsl_vector_get(fact_store,i)<<"\n";
	}
}
		
		

