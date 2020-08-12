#include <iostream>
#include "databasetools.h"
#include "ridgeregressiontools.h"

using namespace std;

int main(){
	dMat data;
	dVec results;
	ImportDataRR(data, results, "boston_housing.txt");

	dMatMat cvtrain, cvtest;
	dMat cvtrainresults, cvtestresults;
	cout << "1";
	CVrandomSampling(cvtrain, cvtrainresults, cvtest, cvtestresults, data, results);
	cout << "2";
	dVec weights(data[0].size());
	cout << "3";
	for (int k = 1; k < 100; k++) {
		cout << k << ",";
		RR_GD_iteration(weights, cvtrain[0], cvtrainresults[0], 10 / (k + 1), 0.1);
	}
	cout << "r^2 is: " << Rsquared(cvtest[0], cvtestresults[0], weights);

	return 0;
}