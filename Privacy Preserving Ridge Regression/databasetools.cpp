#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <typeinfo>

#include "math.h"
#include "seal/seal.h"
#include "time.h"
#include "databasetools.h"

using namespace std;
using namespace seal;

int ImportData(dMat& Matrix, string filename) {
	//open file
	ifstream inFile;
	inFile.open(filename);
	//check file is open
	if (!inFile) {
		cout << "unable to open file";
		return 0;
	}

	string line;
	char split_char = '\t';
	int long ncolumns{};
	//process first row, moving class to the front and extracting number of columns
	if (getline(inFile, line)) {
		istringstream split(line);
		vector<string> record;
		for (string entry; getline(split, entry, split_char); record.push_back(entry));
		ncolumns = record.size();
		vector<double> entry1;
		entry1.push_back(stod(record[1.*(ncolumns - 1)]) * 2 - 1);
		for (int i = 0; i < ncolumns - 1; i++) entry1.push_back(stod(record[i]));

		//add to matrix
		Matrix.push_back(entry1);
	}
	else {
		cout << "could not read file" << exit;
	}
	//process rest of the data
	while (getline(inFile, line)) {
		istringstream split(line);
		vector<string> record;
		for (string entry; getline(split, entry, split_char); record.push_back(entry));
		//record should have the same number of features
		if (record.size() != ncolumns) {
			cout << "database dimension error" << exit;
		}
		//define a new entry
		vector<double> entryi;
		entryi.push_back(stod(record[ncolumns - 1]) * 2 - 1);
		for (int i = 0; i < ncolumns - 1; i++) entryi.push_back(stod(record[i]));
		//add it to the matrix
		Matrix.push_back(entryi);
	}

	return 0;
}
int ImportDataRR(dMat& Matrix, dVec& results, string filename,char split_char,bool ylast) {
	//this file imports data and prepares it for Ridge Regression, putting the covariates into a matrix
	//with rows of the form(1,x1,x2,...,xd), and creates a results vector (y1,...,yn). Split_char is the 
	//character than separates entries in the .txt file, and ylast indicates whether the dependent variable
	//is listed first or last
	//open file
	ifstream inFile;
	inFile.open(filename);
	//check file is open
	if (!inFile) {
		cout << "unable to open file";
		return 0;
	}
	string line;
	
	int ncolumns{};
	
	vector<string> record;
	vector<double> entry1;
	int y = (ylast == true)? 1 : 0;
	//process first row, 
	if (getline(inFile, line)) {
		istringstream split(line);
		for (string entry; getline(split, entry, split_char); record.push_back(entry));
		ncolumns = record.size();
		//push back (1,x1,x2,..,xd)
		entry1.push_back(1.0);
		for (int i = 1-y; i < ncolumns - y; i++) entry1.push_back(stod(record[i]));

		//add to matrix, and put y in Result
		Matrix.push_back(entry1);
		results.push_back(stod(record[y * (ncolumns - 1)]));
		
		record.clear();
		entry1.clear();
	}
	while (getline(inFile, line)) {
		istringstream split(line);
		for (string entry; getline(split, entry, split_char); record.push_back(entry));
		//record should have the same number of features
		if (record.size() != ncolumns) {
			cout << "database dimension error" << exit;
		}
		//define a new entry		
		entry1.push_back(1.0);
		for (int i = 1-y; i < ncolumns - y; i++) entry1.push_back(stod(record[i]));
		
		//add entry to matrix, and result to vector
		Matrix.push_back(entry1);
		results.push_back(stod(record[y*(ncolumns - 1)]));
		
		record.clear();
		entry1.clear();
	}
}

double inner_prod(dVec v, dVec u, int start) {
	if (v.size() != u.size()) {
		cout << "error - vectors are different sizes";
		return 0;
	}
	else {
		double prod = 0.0;
		for (int i = start; i < u.size(); i++) prod += v[i] * u[i];
		return prod;
	}

}

void CVrandomSampling(dMatMat& CVtrain, dMat& CVtrainresults, dMatMat& CVtest, dMat& CVtestresults, dMat data, dVec results) {
	/*srand(time(NULL));*/
	dMat train, test;
	dVec resultstemp;
	int n = data.size();
	int m = floor(n / 5);

	int n_test[5];
	n_test[0] = m;
	n_test[1] = m;
	n_test[2] = m;
	n_test[3] = m;
	n_test[4] = n - 4 * m;

	
	//label all pieces of vector as "unchosen"
	dVec sort(n, -1);
	
	//decide where each record will go
	for (int i = 0; i < 5; i++) {
		
		//start a counter
		int counter = 0;
		while (counter < n_test[i]) {
			//sample a random number from [data.size()]
			int j = rand() % data.size();
			//if it's unchosen, add it to the fold
			
			if (sort[j] == -1) {
				sort[j] += 1 * (i + 1);
				
				//now add record to testing fold
				test.push_back(data[j]);
				
				//and add result to resultstemp
				resultstemp.push_back(results[j]);
				
				counter++;

			}
		}
		
		CVtest.push_back(test);
		CVtestresults.push_back(resultstemp);
		test.clear();
		resultstemp.clear();

	}
	//form the training sets. 
	for (int m = 0; m < 5; m++) {
		for (int l = m + 1; l < 5; l++) {
			for (int i = 0; i < n_test[l]; i++) {
				train.push_back(CVtest[l][i]);
				resultstemp.push_back(CVtestresults[l][i]);
			}
		}
		for (int l = 0; l < m; l++) {
			for (int i = 0; i < n_test[l]; i++) {
				train.push_back(CVtest[l][i]);
				resultstemp.push_back(CVtestresults[l][i]);
			}
		}
		CVtrain.push_back(train);
		CVtrainresults.push_back(resultstemp);
		train.clear();
		resultstemp.clear();
	}
}

void scale_fit(dMat& data, dVec& a, dVec& b,double k) {
	//this function takes a matrix with rows in the form (1,x1,x2,...,xd) and scales each column so
	//that every entry is in [0,k]. It keeps these scaling factors aX+b in the vectors a and b so that
	//the same scaling can be applied to the testing data
	double maxtemp, mintemp,temp;
	for (int i = 1; i < data[0].size(); i++) {
		maxtemp = data[0][i];
		mintemp = data[0][i];
		for (int j = 0; j < data.size(); j++) {
			if (data[j][i] < mintemp)mintemp = data[j][i];
			if (data[j][i] > maxtemp)maxtemp = data[j][i];
		}
		
		temp = (maxtemp - mintemp);
		
		a.push_back(k/temp);
		/*cout << "pushing back: " << a[i-1] << "\n";*/
		b.push_back((k*mintemp) / (mintemp-maxtemp));
		/*cout << "pushing back: " << b[i-1] << "\n";*/
		for (int j = 0; j < data.size(); j++) {
			data[j][i] = a[i-1] * data[j][i] + b[i-1];
		}
	}

}
void scale_columns(dMat& data, dVec a, dVec b) {
	//this function takes scaling vectors a and b and maps column i X -> a[i]X+b[i], leaving the first column
	for(int j = 0; j<data.size();j++){
		for (int i = 1; i < data[0].size(); i++) {
			data[j][i] = a[i-1] * data[j][i] + b[i-1];
		}
	}
}

void center(dVec& v,double& mean) {
	//this function centers a results vector: explicitly, it shifts each y by ybar
	mean = v[0];
	for (int i = 1; i < v.size(); i++)mean += v[i];
	mean /= v.size();
	for (int i = 0; i < v.size(); i++)v[i] -= mean;
}

void shift_results(dVec& v,double mean) {
	//this function shifts a results vector by a constant. Used to test performance of a shifted model
	for (int j = 0; j < v.size(); j++)v[j] -= mean;
}
//
//for (int i = 0; i < D.size(); i++) {
//	cout << "(";
//	for (int j = 0; j < data[0].size() - 1; j++)cout << D[i][j] << ",";
//	cout << D[i][data[0].size() - 1] << ")\n";
//}

bool is_numeric(std::string const& str)
{
	auto result = double();
	auto i = std::istringstream(str);

	i >> result;

	return !i.fail() && i.eof();
}

bool is_number(const std::string& s) {
	return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
}




