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
	gsl_vector_set(xx,0,od.i);	//Driver index
	gsl_vector_set(xx,1,od.x);	//Frequency at cross
	gsl_vector_set(xx,2,od.e);	//Vertical movement
	gsl_vector_set(xx,3,od.rho);	//Recombination rate
	gsl_vector_set(xx,4,od.xf);	//Driver frequency after selection
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
		//cout << "i " << gsl_vector_get(xx,0) << " x " << gsl_vector_get(xx,1) << " e " << gsl_vector_get(xx,2) << " rho " << gsl_vector_get(xx,3) << " xf " << gsl_vector_get(xx,4) << "\n";
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
	double i=gsl_vector_get(xx,0);
	double x=gsl_vector_get(xx,1);
	double e=gsl_vector_get(xx,2);
	double rho=gsl_vector_get(xx,3);
	double xf=gsl_vector_get(xx,4);
	if (rho<0) {
		//cout << "Rho " << rho << "\n";
		chk=0;
	}
	if (i<0 || i>p[0]) {
		//cout << "i " << i << "\n";
		chk=0;
	}
	if (x<=0.5 || x>=1) {
		//cout << "x " << x << "\n";
		chk=0;
	}
	if (e< -pow(x,2) || e> pow((1-x),2)) {
		//cout << "e " << e << " " << -pow(x,2) << " " << pow((1-x),2) << "\n";
		chk=0;
	}
	if (xf<pow(x,2) || xf > 1-pow(1-x,2)) {
		//cout << "xf " << xf << " " << pow(x,2) << " " << 1-pow(1-x,2) << "\n";
		chk=0;
	}
	if (xf+e>=1 || xf+e <=0) {
		//cout << "xf+e " << xf+e << "\n";
		chk=0;
	}
	return chk;
}

void FindFrequencies (vector<rec> traj, vector<rec>& traj_inf, const gsl_vector *xx) {
	//Set up traj_inf
	traj_inf=traj;
	//Frequency at driver locus
	unsigned int i=floor(gsl_vector_get(xx,0)); //Index of driver
	//cout << i << " " << traj_inf.size() << " " << traj.size() << "\n";
	double x=gsl_vector_get(xx,1);
	double e=gsl_vector_get(xx,2);
	double rho=gsl_vector_get(xx,3);
	double xf=gsl_vector_get(xx,4);   
	//	double xf=gsl_vector_get(xx,4);   //?? Sort out something around here?
	
	//cout << "FindFrequencies " << i << " " << x << " " << e << " " << rho << " " << xf << "\n";

	
	traj_inf[i].q=xf+e;
	
	//Frequency at passenger loci
	for (unsigned int j=0;j<traj.size();j++) {
		//cout << j << "\n";
		if (j!=i) {
			double Delta_ij=abs(traj_inf[i].pos-traj[j].pos);
			if (Delta_ij<0) {
				cout << "Error\n";
			}
			double ExpRhoD=exp((-1)*rho*Delta_ij);
			//traj_inf[j].q= ((1+ExpRhoD)*(2*x*(1-x)*xf))+((1-ExpRhoD)*(1-(2*x*(1-x)*xf)))+e;
			traj_inf[j].q= ((x+(0.5*(1-x)*(1+ExpRhoD)))*xf) + ((0.5*x*(1-ExpRhoD))*(1-xf)) + e;
			if (traj_inf[j].q>1) {
				cout << "Error2 " << traj_inf[i].pos << " " << traj_inf[j].pos << " " << x << " " << Delta_ij << " " << rho << " " << ExpRhoD << " " << xf << " " << e << " " << traj_inf[j].q << "\n"; ;
			}
		}
	}
}

void FindFrequenciesFinal (vector<rec> traj, vector<rec>& traj_inf, const gsl_vector *xx, ofstream& out_file) {
	//Set up traj_inf
	traj_inf=traj;
	//cout << "i "  << gsl_vector_get(xx,0) << " x " << gsl_vector_get(xx,1) << " e " << gsl_vector_get(xx,2) << " rho " << gsl_vector_get(xx,3) << " xf " << gsl_vector_get(xx,4) << "\n";

	//Frequency at driver locus
	unsigned int i=floor(gsl_vector_get(xx,0)); //Index of driver
	//cout << i << " " << traj_inf.size() << " " << traj.size() << "\n";
	double x=gsl_vector_get(xx,1);
	double e=gsl_vector_get(xx,2);
	double rho=gsl_vector_get(xx,3);
	double xf=gsl_vector_get(xx,4);   
	
	traj_inf[i].q=xf+e;
	//cout << traj_inf[i].q << "\n";
	out_file << i << " " << traj[i].pos << " " << traj[i].q << " " << traj_inf[i].q << "\n";
	
	//Frequency at passenger loci
	for (unsigned int j=0;j<traj.size();j++) {
		//cout << j << "\n";
		if (j!=i) {
			double Delta_ij=abs(traj_inf[i].pos-traj[j].pos);
			if (Delta_ij<0) {
				cout << "Error\n";
			}
			double ExpRhoD=exp((-1)*rho*Delta_ij);
			//traj_inf[j].q= ((1+ExpRhoD)*(2*x*(1-x)*xf))+((1-ExpRhoD)*(1-(2*x*(1-x)*xf)))+e;
			traj_inf[j].q= ((x+(0.5*(1-x)*(1+ExpRhoD)))*xf) + ((0.5*x*(1-ExpRhoD))*(1-xf)) + e;
			out_file << j << " " << traj[j].pos << " " << traj[j].q << " " << traj_inf[j].q << "\n";
			if (traj_inf[j].q>1) {
				cout << "Error3 " << x << " " << Delta_ij << " " << ExpRhoD << " " << xf << " " << e << "\n"; ;
			}
		}
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
	double bin=1e10;
	if ((alpha==0)||(beta==0)) {
		bin=1e10;
	} else {
		bin=gsl_sf_lngamma(N+1)+gsl_sf_lngamma(r+alpha)+gsl_sf_lngamma(N-r+beta)+gsl_sf_lngamma(alpha+beta)-gsl_sf_lngamma(r+1)-gsl_sf_lngamma(N-r+1)-gsl_sf_lngamma(N+alpha+beta)-gsl_sf_lngamma(alpha)-gsl_sf_lngamma(beta);
		bin=exp(bin);
	}
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
		
		

