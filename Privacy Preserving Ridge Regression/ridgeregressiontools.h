#pragma once
#include "databasetools.h"


double Rsquared(dMat data, dVec results, dVec weights);
void RR_GD_iteration(dVec& weights, dMat data, dVec results, double alpha, double lambda);
void RR_NAG_iteration(dVec& beta, dVec& v, dMat data, dVec results, double alpha, double gamma, double lambda);
void RR_fixed_Hessian_iteration(dVec& beta, dMat data, dVec results, dVec H, double lambda);