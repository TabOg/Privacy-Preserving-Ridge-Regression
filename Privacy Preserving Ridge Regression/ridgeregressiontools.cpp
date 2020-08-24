#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "math.h"
#include "seal/seal.h"
#include "time.h"
#include "databasetools.h"
#include <valarray>

using namespace std;
using namespace seal;


double Rsquared(dMat data, dVec results, dVec weights) {
	//data should be fed to this function with each record in the form (1,x1,...,xd) and the
	//results vector separate		
	//first calculate mean
	double mean = 0;
	for (int i = 0; i < results.size(); i++)mean += results[i] / data.size();
	/*cout << "mean is:" << mean << "\n";*/
	//now initialise SS_tot and SS_res & loop over dataset
	double SStot = 0;
	double SSres = 0;
	for (int i = 0; i < data.size(); i++) {
		SStot += (results[i] - mean) * (results[i] - mean);
		
		SSres += (results[i] - inner_prod(weights, data[i])) * (results[i] - inner_prod(weights, data[i]));
	}

	/*cout << "total sum of squares: " << SStot << "\n";
	cout << "residual sum of squares: " << SSres << "\n";*/
	return  100*(1 - (SSres / SStot));

}

void RR_GD_iteration(dVec& weights, dMat data, dVec results, double alpha, double lambda) {
	//processes data with entries of the form (1,x1,x2,..,xd)
	dVec grad(data[0].size(),0.0);
	double temp;	
	for (int i = 0; i < data.size(); i++) {
		temp = results[i] - inner_prod(data[i], weights);
		for (int j = 0; j < data[0].size(); j++) {
			grad[j] += temp * data[i][j];
		}
	}

	weights[0] = weights[0] + (grad[0]) * (alpha);
	for (int j = 1; j < data[0].size(); j++) {
		weights[j] =((1-(alpha)*lambda)*weights[j] + (alpha) * (grad[j]));
	}
}

void RR_NAG_iteration(dVec& beta, dVec& v, dMat data, dVec results, double alpha, double gamma, double lambda) {
	dVec grad(data[0].size(), 0.0);
	double temp;
	dVec weights = beta;

	for (int i = 0; i < data.size(); i++) {
		temp = results[i] - inner_prod(data[i], v);
		for (int j = 0; j < data[0].size(); j++) {
			grad[j] += temp * data[i][j];
		}
	}
	beta[0] = v[0] + (grad[0]) * (alpha);
	v[0] = (1 - gamma) * beta[0] + gamma * weights[0];

	for (int j = 1; j < data[0].size(); j++) {
		beta[j] = ((1 - (alpha)*lambda) * v[j] + (alpha) * (grad[j]));
		v[j] = (1 - gamma) * beta[j] + gamma * weights[j];
	}	
}

void RR_fixed_Hessian_iteration(dVec& beta, dMat data, dVec results, dVec H, double lambda) {
	dVec grad(data[0].size(), 0.0);
	double temp;

	for (int i = 0; i < data.size(); i++) {
		temp = results[i] - inner_prod(data[i], beta);
		for (int j = 0; j < data[0].size(); j++) {
			grad[j] += temp * data[i][j];
		}
	}
	beta[0] = beta[0] + (grad[0]) * (H[0]);
	for (int j = 1; j < data[0].size(); j++) {

		beta[j] = ((1 - H[j] * lambda) * beta[j] + H[j] * (grad[j]));
	}
}