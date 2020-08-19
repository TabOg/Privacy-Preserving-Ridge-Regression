#pragma once
#include "databasetools.h"
#include "ridgeregressiontools.h"

int RR_gradient_descent(string file, char split_char, bool y_last, double lambda,double alpha, int iternum);
int RR_NAG_descent(string file, char split_char, bool y_last, double lambda,double alpha, int iternum);
int RR_fixed_Hessian(string file, char split_char, bool y_last, double lambda, int iternum);
