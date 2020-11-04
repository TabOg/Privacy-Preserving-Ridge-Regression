#pragma once

#include "seal/seal.h"
#include <iostream>
#include "databasetools.h"
#include "ridgeregressiontools.h"

using namespace std;
using namespace seal;

int PP_Gradient_Descent_RR(double alpha,double lambda, bool bit_size);
int PP_Nesterov_Gradient_Descent_RR(double alpha, double lambda, bool bit_size);
int PP_Fixed_Hessian_RR(double lambda, bool bit_size);
