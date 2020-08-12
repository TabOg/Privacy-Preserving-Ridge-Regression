#pragma once
#include "databasetools.h"


double Rsquared(dMat data, dVec results, dVec weights);
void RR_GD_iteration(dVec& weights, dMat data, dVec results, double alpha, double lambda);