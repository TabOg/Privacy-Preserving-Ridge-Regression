# Privacy-Preserving-Ridge-Regression

This repo implements privacy preserving ridge regression training using homomorphic encryption, and compares three different minimisation algorithms, Gradient Descent, Nesterov Accelerated Gradient Descent, and a fixed Hessian version of the Newton Raphson method. All methods are implemented with 5-fold cross validation on the Boston Housing dataset. The code is built on top of SEAL (https://github.com/Microsoft/SEAL) which will need to be installed separately and then linked.
