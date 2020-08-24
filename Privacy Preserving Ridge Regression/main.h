#pragma once

#include "seal/seal.h"
#include <iostream>
#include "databasetools.h"
#include "ridgeregressiontools.h"

using namespace std;
using namespace seal;

int Gradient_Descent_RR(double alpha,double lambda);
