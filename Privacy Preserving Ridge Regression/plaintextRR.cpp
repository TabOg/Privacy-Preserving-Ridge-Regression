#include <iostream>
#include "databasetools.h"
#include "ridgeregressiontools.h"


int RR_gradient_descent(string file, char split_char, bool y_last,double lambda,double alpha, int iternum) {
	dMat data;
	dVec results;
	ImportDataRR(data, results, file, split_char, y_last);

	dMatMat cvtrain, cvtest;
	dMat cvtrainresults, cvtestresults;

	//CVrandomSampling(cvtrain, cvtrainresults, cvtest, cvtestresults, data, results);
	//cout << "fold sizes are:\n1 -- " << cvtrain[0].size() << "--" << cvtrainresults[0].size() << "\n";
	//cout << "2 -- " << cvtrain[1].size() << "--" << cvtrainresults[1].size() << "\n";
	//cout << "3 -- " << cvtrain[2].size() << "--" << cvtrainresults[2].size() << "\n";
	//cout << "4 -- " << cvtrain[3].size() << "--" << cvtrainresults[3].size() << "\n";
	//cout << "5 -- " << cvtrain[4].size() << "--" << cvtrainresults[4].size() << "\n";

	//cout << "test sizes are:\n1 -- " << cvtest[0].size() << "--" << cvtestresults[0].size() << "\n";
	//cout << "2 -- " << cvtest[1].size() << "--" << cvtestresults[1].size() << "\n";
	//cout << "3 -- " << cvtest[2].size() << "--" << cvtestresults[2].size() << "\n";
	//cout << "4 -- " << cvtest[3].size() << "--" << cvtestresults[3].size() << "\n";
	//cout << "5 -- " << cvtest[4].size() << "--" << cvtestresults[4].size() << "\n";

	dVec weights(data[0].size(), 0.0);
	dVec a, b;
	double mean;
	double avgr2 = 0;

	for (int j = 0; j < 5; j++) {
		fill(weights.begin(), weights.end(), 0.0);
		mean = 0;
		center(cvtrainresults[j], mean);
		scale_fit(cvtrain[j], a, b);
		scale_columns(cvtest[j], a, b);
		shift_results(cvtestresults[j], mean);
		/*cout << "iteration " << 0 << " r^2 : " << Rsquared(cvtrain[j], cvtrainresults[j], weights) << "%\n";*/
		for (int k = 0; k < iternum; k++) {
			RR_GD_iteration(weights, cvtrain[j], cvtrainresults[j], alpha, lambda);
			cout << "fold" << j + 1 << ", iteration " << k + 1 << " CV r^2: " << Rsquared(cvtest[j], cvtestresults[j], weights) << "%\n";
			//cout << "(";
			//for (int i = 0; i < weights.size() - 1; i++)cout << weights[i] << ",";
			//cout << weights[weights.size() - 1] << ")\n";
			/*cout << "iteration" << k + 1<< "r^2:" << Rsquared(cvtrain[j], cvtrainresults[j], weights) << "%\n";*/
		}
		cout << "fold " << j << " CV r^2: " << Rsquared(cvtest[j], cvtestresults[j], weights) << "%\n";
		avgr2 += Rsquared(cvtest[j], cvtestresults[j], weights);
		a.clear();
		b.clear();		
	}
	avgr2 /= 5;
	cout << "average cross validation accuracy: " << avgr2 << "%\n";
	weights.clear();
	return 0;
}
int RR_NAG_descent(string file, char split_char, bool y_last, double lambda, double alpha, int iternum) {
	dMat data;
	dVec results;
	ImportDataRR(data, results, file, split_char, y_last);

	dMatMat cvtrain, cvtest;
	dMat cvtrainresults, cvtestresults;

	CVrandomSampling(cvtrain, cvtrainresults, cvtest, cvtestresults, data, results);
	//cout << "fold sizes are:\n1 -- " << cvtrain[0].size() << "--" << cvtrainresults[0].size() << "\n";
	//cout << "2 -- " << cvtrain[1].size() << "--" << cvtrainresults[1].size() << "\n";
	//cout << "3 -- " << cvtrain[2].size() << "--" << cvtrainresults[2].size() << "\n";
	//cout << "4 -- " << cvtrain[3].size() << "--" << cvtrainresults[3].size() << "\n";
	//cout << "5 -- " << cvtrain[4].size() << "--" << cvtrainresults[4].size() << "\n";

	//cout << "test sizes are:\n1 -- " << cvtest[0].size() << "--" << cvtestresults[0].size() << "\n";
	//cout << "2 -- " << cvtest[1].size() << "--" << cvtestresults[1].size() << "\n";
	//cout << "3 -- " << cvtest[2].size() << "--" << cvtestresults[2].size() << "\n";
	//cout << "4 -- " << cvtest[3].size() << "--" << cvtestresults[3].size() << "\n";
	//cout << "5 -- " << cvtest[4].size() << "--" << cvtestresults[4].size() << "\n";

	dVec beta(data[0].size(), 0.0);
	dVec v = beta;
	dVec a, b;
	double mean = 0;
	double gamma, T;
	double t = 1.;
	double avgr2 = 0;
	dVec sum;
	for (int j = 0; j < 5; j++) {
		fill(beta.begin(), beta.end(), 0.0);
		fill(v.begin(), v.end(), 0.0);
		mean = 0;
		t = 1;
		center(cvtrainresults[j], mean);
		scale_fit(cvtrain[j], a, b);
		cout << "iteration " << 0 << " r^2 : " << Rsquared(cvtrain[j], cvtrainresults[j], v) << "%\n";
		for (int k = 0; k < iternum; k++) {
			T = (1. + sqrt(1. + 4 * t * t)) / 2.;
			gamma = (1 - t) / T;
			t = T;
			RR_NAG_iteration(beta, v, cvtrain[j], cvtrainresults[j], alpha, gamma, lambda);
		}
	/*	cout << "iteration "<<iternum<<" r^2: " << Rsquared(cvtrain[j], cvtrainresults[j], v) << "%\n";*/
		scale_columns(cvtest[j], a, b);
		shift_results(cvtestresults[j], mean);
		cout <<"fold "<< j << " r^2: " << Rsquared(cvtest[j], cvtestresults[j], v) << "%\n";
		avgr2 += Rsquared(cvtest[j], cvtestresults[j], v);
		a.clear();
		b.clear();
	}
	avgr2 /= 5;
	cout << "average cross validation r^2: " << avgr2 << "%";
}

int RR_fixed_Hessian(string file, char split_char, bool y_last, double lambda, int iternum) {
	dMat data;
	dVec results;
	ImportDataRR(data, results, file, split_char, y_last);

	dMatMat cvtrain, cvtest;
	dMat cvtrainresults, cvtestresults;

	CVrandomSampling(cvtrain, cvtrainresults, cvtest, cvtestresults, data, results);
	cout << "fold sizes are:\n1 -- " << cvtrain[0].size() << "--" << cvtrainresults[0].size() << "\n";
	cout << "2 -- " << cvtrain[1].size() << "--" << cvtrainresults[1].size() << "\n";
	cout << "3 -- " << cvtrain[2].size() << "--" << cvtrainresults[2].size() << "\n";
	cout << "4 -- " << cvtrain[3].size() << "--" << cvtrainresults[3].size() << "\n";
	cout << "5 -- " << cvtrain[4].size() << "--" << cvtrainresults[4].size() << "\n";

	cout << "test sizes are:\n1 -- " << cvtest[0].size() << "--" << cvtestresults[0].size() << "\n";
	cout << "2 -- " << cvtest[1].size() << "--" << cvtestresults[1].size() << "\n";
	cout << "3 -- " << cvtest[2].size() << "--" << cvtestresults[2].size() << "\n";
	cout << "4 -- " << cvtest[3].size() << "--" << cvtestresults[3].size() << "\n";
	cout << "5 -- " << cvtest[4].size() << "--" << cvtestresults[4].size() << "\n";

	dVec beta(data[0].size(), 0.0);
	dVec a, b, H;
	double mean = 0;
	double avgr2 = 0;
	double alpha;
	dVec sum;
	for (int j = 0; j < 5; j++) {
		fill(beta.begin(), beta.end(), 0.0);
		mean = 0;
		center(cvtrainresults[j], mean);
		scale_fit(cvtrain[j], a, b);
		H.clear();
		
		for (int i = 0; i < cvtrain[j].size(); i++) {
			alpha = 0;
			for (int k = 0; k < cvtrain[j][0].size(); k++) {
				alpha += cvtrain[j][i][k];
			}
			sum.push_back(alpha);
		}

		alpha = 0;
		for (int i = 0; i < cvtrain[j].size(); i++) {
			alpha += sum[i];
		}
		cout << alpha << ",";
		H.push_back(1. / alpha);
		for (int i = 1; i < cvtrain[j][0].size(); i++) {
			alpha = 0;
			for (int k = 0; k < cvtrain[j].size(); k++) {
				alpha += cvtrain[j][k][i] * sum[k];
			}
			alpha += lambda;
			cout << alpha << ",";
			H.push_back(1. / alpha);
		}
		cout << "\n";
		cout << "(";
		for (int i = 0; i < H.size() - 1; i++)cout << H[i] << ",";
		cout << H[H.size() - 1] << ")\n";

		for (int k = 0; k < iternum; k++) {
			RR_fixed_Hessian_iteration(beta, cvtrain[j], cvtrainresults[j], H, lambda);
			cout << "iteration " << k + 1 << " r^2: " << Rsquared(cvtrain[j], cvtrainresults[j], beta) << "%\n";
		}
		cout << "iteration " << iternum << " r^2: " << Rsquared(cvtrain[j], cvtrainresults[j], beta) << "%\n";
		scale_columns(cvtest[j], a, b);
		shift_results(cvtestresults[j], mean);
		cout << j << "th fold r^2: " << Rsquared(cvtest[j], cvtestresults[j], beta) << "%\n";
		avgr2 += Rsquared(cvtest[j], cvtestresults[j], beta);
		a.clear();
		b.clear();
	}
	avgr2 /= 5;
	cout << "average cross validation accuracy: " << avgr2 << "%";
	return 0;
}
