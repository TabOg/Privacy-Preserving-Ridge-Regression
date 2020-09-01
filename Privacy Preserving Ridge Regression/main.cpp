#include <iostream>
#include "databasetools.h"
#include "ridgeregressiontools.h"
#include "plaintextRR.h"
#include "seal/seal.h"
#include "main.h"

using namespace std;
using namespace seal;

int main() {

    cout << "+---------------------------------------------------------+" << endl;
    cout << "| Select a procedure to run:                              |" << endl;
    cout << "+---------------------------------------------------------+" << endl;
    cout << "| Methods                    | Source Files               |" << endl;
    cout << "+----------------------------+----------------------------+" << endl;
    cout << "| 1. Plaintext Gradient      | plaintextRR.cpp            |" << endl;
    cout << "|    Descent                 |                            |" << endl;
    cout << "| 2. Plaintext Nesterov      | plaintextRR.cpp            |" << endl;
    cout << "|    Accelerated GD          |                            |" << endl;
    cout << "| 3. Plaintext Fixed Hessian | plaintextRR.cpp            |" << endl;
    cout << "|    Taylor approximation    |                            |" << endl;
    cout << "| 4. Encrypted Gradient      | GD.cpp                     |" << endl;
    cout << "|    Descent                 |                            |" << endl;
    cout << "| 5. Encrypted Nesterov      | NAG.cpp                    |" << endl;
    cout << "|    Accelerated GD          |                            |" << endl;
    cout << "| 6. Encrypted Fixed Hessian | FH.cpp                     |" << endl;
    cout << "+----------------------------+----------------------------+" << endl;
    int selection = 0;
    bool invalid = true;
    string alpha, lambda,iternum;
    do
    {
        cout << endl << "> Run example (1 ~ 6) or exit (0): ";
        if (!(cin >> selection))
        {
            invalid = false;
        }
        else if (selection < 0 || selection > 7)
        {
            invalid = false;
        }
        else
        {
            invalid = true;
        }
        if (!invalid)
        {
            cout << "  [Beep~~] Invalid option: type 0 ~ 6" << endl;
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
        }
    } while (!invalid);

    switch (selection)
    {
    case 1:
        cout << "alpha = ";
        cin >> alpha;
        while (!is_numeric(alpha)) {
                cout << "learning rate must be a real number!\n";
                cout << "alpha = ";
                cin >> alpha;
        }
        cout << "lambda = ";
        cin >> lambda;
        while (!is_numeric(lambda)) {
            cout << "regularisation parameter must be a real number!\n";
            cout << "lambda = ";
            cin >> lambda;
        }
        cout << "number of iterations = ";
        cin >> iternum;
        while (!is_number(iternum)) {
            cout << "iteration number must be an integer!\n";
            cout << "number of iterations = ";
            cin >> iternum;
        }
        RR_gradient_descent("boston_housing.txt", ',', true, stod(lambda), stod(alpha), stoi(iternum));
        
        break;

    case 2:
        cout << "alpha = ";
        cin >> alpha;
        while (!is_numeric(alpha)) {
            cout << "learning rate must be a real number!\n";
            cout << "alpha = ";
            cin >> alpha;
        }
        cout << "lambda = ";
        cin >> lambda;
        while (!is_numeric(lambda)) {
            cout << "regularisation parameter must be a real number!\n";
            cout << "lambda = ";
            cin >> lambda;
        }
        cout << "number of iterations = ";
        cin >> iternum;
        while (!is_number(iternum)) {
            cout << "iteration number must be an integer!\n";
            cout << "number of iterations = ";
            cin >> iternum;
        }
        RR_NAG_descent("boston_housing.txt", ',', true, stod(lambda),stod(alpha), stoi(iternum));
        break;

    case 3:
        cout << "lambda = ";
        cin >> lambda;
        while (!is_numeric(lambda)) {
            cout << "regularisation parameter must be a real number!\n";
            cout << "lambda = ";
            cin >> lambda;
        }
        cout << "number of iterations = ";
        cin >> iternum;
        while (!is_number(iternum)) {
            cout << "iteration number must be an integer!\n";
            cout << "number of iterations = ";
            cin >> iternum;
        }
        RR_fixed_Hessian("boston_housing.txt", ',', true, stod(lambda), stoi(iternum));
        break;
    
    case 4:
        cout << "alpha = ";
        cin >> alpha;
        while (!is_numeric(alpha)) {
            cout << "learning rate must be a real number!\n";
            cout << "alpha = ";
            cin >> alpha;
        }
        cout << "lambda = ";
        cin >> lambda;
        while (!is_numeric(lambda)) {
            cout << "regularisation parameter must be a real number!\n";
            cout << "lambda = ";
            cin >> lambda;
        }
        break;
        PP_Gradient_Descent_RR(stod(alpha), stod(lambda));
    case 5:
        cout << "alpha = ";
        cin >> alpha;
        while (!is_numeric(alpha)) {
            cout << "learning rate must be a real number!\n";
            cout << "alpha = ";
            cin >> alpha;
        }
        cout << "lambda = ";
        cin >> lambda;
        while (!is_numeric(lambda)) {
            cout << "regularisation parameter must be a real number!\n";
            cout << "lambda = ";
            cin >> lambda;
        }
        PP_Nesterov_Gradient_Descent_RR(stod(alpha), stod(lambda));
        break;

    case 6:
        cout << "lambda = ";
        cin >> lambda;
        while (!is_numeric(lambda)) {
            cout << "regularisation parameter must be a real number!\n";
            cout << "lambda = ";
            cin >> lambda;
        }
        PP_Fixed_Hessian_RR(stod(lambda));
        break;

    case 0:
        return 0;
    }
}
		
	
