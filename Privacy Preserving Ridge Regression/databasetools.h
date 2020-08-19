#ifndef DATABASE_TOOLS
#define DATABASE_TOOLS

#include <vector>
#include <iostream>
#include <string>
#pragma once
#include <fstream>
#include <sstream>

#include "math.h"
#include "seal/seal.h"
using namespace std;
using namespace seal;

typedef vector<double> dVec;
typedef vector<vector<double>> dMat;
typedef vector<Ciphertext> cVec;
typedef vector<Plaintext> PVec;
typedef vector<dMat> dMatMat;

int ImportData(dMat& Z, string filename);
int ImportDataRR(dMat& Matrix, dVec& results, string filename, char split_char = '\t', bool ylast = true);
double inner_prod(dVec v, dVec u, int start = 0);
void CVrandomSampling(dMatMat& CVtrain, dMat& CVtrainresults, dMatMat& CVtest, dMat& CVtestresults, dMat data, dVec results);
void scale_fit(dMat& data, dVec& a, dVec& b,double k =1);
void scale(dMat& data, dVec a, dVec b);
void center(dVec& v, double& mean);
void shift_results(dVec& v, double mean);
#endif