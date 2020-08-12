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
	int long ncolumns;
	//process first row, moving class to the front and extracting number of columns
	if (getline(inFile, line)) {
		istringstream split(line);
		vector<string> record;
		for (string entry; getline(split, entry, split_char); record.push_back(entry));
		ncolumns = record.size();
		vector<double> entry1;
		entry1.push_back(stod(record[ncolumns - 1]) * 2 - 1);
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
int ImportDataRR(dMat& Matrix, dVec results, string filename) {
	//open file
	ifstream inFile;
	inFile.open(filename);
	//check file is open
	if (!inFile) {
		cout << "unable to open file";
		return 0;
	}
	string line;
	char split_char = ',';
	int long ncolumns{};
	vector<string> record;
	vector<double> entry1;
	//process first row, 
	if (getline(inFile, line)) {
		istringstream split(line);
		for (string entry; getline(split, entry, split_char); record.push_back(entry));
		ncolumns = record.size();
		//push back (1,x1,x2,..,xd)
		entry1.push_back(1);
		for (int i = 0; i < ncolumns - 1; i++) entry1.push_back(stod(record[i]));

		//add to matrix, and put y in Result
		Matrix.push_back(entry1);
		results.push_back(stod(record[1 * (ncolumns - 1)]));
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
		entry1.push_back(1);
		for (int i = 0; i < ncolumns - 1; i++) entry1.push_back(stod(record[i]));

		//add entry to matrix, and result to vector
		Matrix.push_back(entry1);
		results.push_back(stod(record[1 * (ncolumns - 1)]));
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
	cout << "1";
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

	cout << "2";
	//label all pieces of vector as "unchosen"
	dVec sort(n, -1);
	cout << "3";
	//decide where each record will go
	for (int i = 0; i < 5; i++) {
		cout << "4";
		//start a counter
		int counter = 0;
		while (counter < n_test[i]) {
			//sample a random number from [data.size()]
			int j = rand() % data.size();
			//if it's unchosen, add it to the fold
			cout << "5";
			if (sort[j] == -1) {
				sort[j] += 1 * (i + 1);
				cout << "6";
				//now add record to testing fold
				test.push_back(data[j]);
				cout << "7";
				//and add result to resultstemp
				resultstemp.push_back(results[j]);
				cout << "8";
				counter++;

			}
		}
		cout << "6";
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


