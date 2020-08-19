#include <iostream>
#include "databasetools.h"
#include "ridgeregressiontools.h"
#include "plaintextRR.h"

using namespace std;

int main() {
	RR_gradient_descent("boston_housing.txt", ',', true, 1, 0.0009, 15);
	RR_NAG_descent("boston_housing.txt", ',', true, 1, 0.0009, 15);
	RR_fixed_Hessian("boston_housing.txt", ',', true, 1, 15);
	return 0;
}
		
	
