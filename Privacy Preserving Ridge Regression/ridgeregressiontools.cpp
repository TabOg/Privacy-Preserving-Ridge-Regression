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
	int n = data.size();
	int nfeatures = data[0].size();
	//first calculate mean
	double mean = 0;
	for (int i = 0; i < results.size(); i++)mean += results[i] / n;

	//now initialise SS_tot and SS_res & loop over dataset
	double SStot = 0;
	double SSres = 0;
	for (int i = 0; i < n; i++) {
		SStot += (results[i] - mean) * (results[i] - mean);
		SSres += (results[i] - inner_prod(weights, data[i])) * (results[i] - inner_prod(weights, data[i]));
	}
	return  1 - SSres / SStot;

}

void RR_GD_iteration(dVec& weights, dMat data, dVec results, double alpha, double lambda) {
	//need to data with entries of the form (1,x1,x2,..,xd)
	dVec grad(data[0].size());
	double temp;
	for (int i = 0; i < data.size(); i++) {
		temp = results[i] - inner_prod(data[i], weights);
		for (int j = 0; j < data[0].size(); j++) {
			grad[j] += temp * data[i][j];
		}
	}
	weights[0] += 2 * alpha * grad[0];
	for (int j = 1; j < data[0].size(); j++) {
		weights[j] += 2 * alpha * (grad[j] - lambda * weights[j]);
	}
}