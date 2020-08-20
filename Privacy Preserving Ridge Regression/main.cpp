#include <iostream>
#include "databasetools.h"
#include "ridgeregressiontools.h"
#include "plaintextRR.h"
#include "seal/seal.h"

using namespace std;
using namespace seal;

int main() {
    EncryptionParameters parms(scheme_type::CKKS);
    size_t poly_modulus_degree = 16384;

    parms.set_poly_modulus_degree(poly_modulus_degree);

    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 50,40,40,50 }));
    cout << "Generating context...";
    auto start = chrono::steady_clock::now();
    auto context = SEALContext::Create(parms);
    KeyGenerator keygen(context);
    PublicKey public_key = keygen.public_key();
    SecretKey secret_key = keygen.secret_key();
    RelinKeys relin_keys = keygen.relin_keys_local();
    GaloisKeys gal_keys = keygen.galois_keys_local();

    Encryptor encryptor(context, public_key);
    Decryptor decryptor(context, secret_key);
    Evaluator evaluator(context);

    auto end = chrono::steady_clock::now();
    auto diff = end - start;
    cout << "KeyGen time = " << chrono::duration <double, milli>(diff).count() / 1000.0 << " s \n";
    CKKSEncoder encoder(context);
    size_t slot_count = encoder.slot_count();
    cout << "Number of slots: " << slot_count << "\n";

    dMat data;
    dVec results;
    ImportDataRR(data, results, "boston_housing.txt", ',', true);
    
    double alpha = 0.0009;

    //define the matrix C, which is zero apart from alphas on the diagonal
    dMat C;
    for (int i = 0; i < data[0].size(); i++) {
        dVec ctemp(data[0].size(), 0.0);
        ctemp[i] += alpha;
        C.push_back(ctemp);
    }

    //set 40 bits of precision
    double scale = pow(2.0, 40);

    //encode the rows of the matrix C
    pVec Cplain;
    Plaintext ptemp;
    for (int i = 0; i < C.size(); i++) {
        encoder.encode(C[i], scale, ptemp);
        Cplain.push_back(ptemp);
    }
    dMat D;
    dVec Dtemp;
    for (int i = 0; i < C.size(); i++) {
        encoder.decode(Cplain[i], Dtemp);
        D.push_back(Dtemp);
    }

    dMatMat cvtrain, cvtest;
    dMat cvtrainresults, cvtestresults;

    CVrandomSampling(cvtrain, cvtrainresults, cvtest, cvtestresults, data, results);
    cout << "fold sizes are:\n1 -- " << cvtrain[0].size() << "--" << cvtrainresults[0].size() << "\n";
    cout << "2 -- " << cvtrain[1].size() << "--" << cvtrainresults[1].size() << "\n";
    cout << "3 -- " << cvtrain[2].size() << "--" << cvtrainresults[2].size() << "\n";
    cout << "4 -- " << cvtrain[3].size() << "--" << cvtrainresults[3].size() << "\n";
    cout << "5 -- " << cvtrain[4].size() << "--" << cvtrainresults[4].size() << "\n";

    cout << "test sizes are:\n1 -- " << cvtest[0].size() << "--" << cvtestresults[0].size() << "\n";
    cout << "2 -- " << cvtest[1].size() << "--" << cvtestresults[1].size() << "\n";
    cout << "3 -- " << cvtest[2].size() << "--" << cvtestresults[2].size() << "\n";
    cout << "4 -- " << cvtest[3].size() << "--" << cvtestresults[3].size() << "\n";
    cout << "5 -- " << cvtest[4].size() << "--" << cvtestresults[4].size() << "\n";
    dVec a, b;
    double mean;

    for (int j = 0; j < 1; j++) {
        a.clear();
        b.clear();
        mean = 0;
        cout << "starting fold" << j + 1 << "\n";
        //start by centering the results, and fitting the scale to the data
        center(cvtrainresults[j], mean);
        scale_fit(cvtrain[j], a, b);
        //store a copy of the matrix C for this fold
        pVec Cj = Cplain;
        //extract dimensions
        int nfeatures = cvtrain[j][0].size();
        int n = cvtrain.size();
        //split the data into columns, leaving the first as it's a column of 1s. This would happen client side
        dVec input;
        pVec dataplain;

        cout << "Encoding...";
        start = chrono::steady_clock::now();
        for (int i = 1; i < nfeatures; i++) {
            input.clear();
            for (int k = 0; k < n; k++)input.push_back(cvtrain[j][k][i]);
            encoder.encode(input, scale, ptemp);
            dataplain.push_back(ptemp);
        }
        Plaintext resultsplain;
        Ciphertext results;

        encoder.encode(cvtrainresults[j], scale, resultsplain);
        end = chrono::steady_clock::now();
        diff = end - start;
        cout << "Encoding time = " << chrono::duration <double, milli>(diff).count() / 1000.0 << " s \n";

        cout << "Encrypting...";
        start = chrono::steady_clock::now();
        cVec dataenc;
        Ciphertext ctemp;
        for (int i = 0; i < dataplain.size(); i++) {
            encryptor.encrypt(dataplain[i], ctemp);
            dataenc.push_back(ctemp);
        }

        Ciphertext resultsenc;
        encryptor.encrypt(resultsplain, resultsenc);
        end = chrono::steady_clock::now();
        diff = end - start;
        cout << "Encrypting time = " << chrono::duration <double, milli>(diff).count() / 1000.0 << " s \n";

        cout << "Forming Y...";
        Ciphertext Y;
        //first we create a zero vector. As the results vector is centered, Y[0]=0.
        encoder.encode(0,scale, ptemp);
        encryptor.encrypt(ptemp, Y);
        evaluator.mod_switch_to_next_inplace(Y);
        evaluator.mod_switch_to_next_inplace(Y);
        //now we add each of the entries
        Ciphertext allsumtemp;
        for (int l = 1; l < nfeatures; l++) {
            //multiply Xl and the results vector
            evaluator.multiply(dataenc[l - 1], resultsenc, ctemp);
            evaluator.relinearize_inplace(ctemp, relin_keys);
            evaluator.rescale_to_next_inplace(ctemp);
            allsumtemp = ctemp;
            for (int k = 0; k < log2(slot_count); k++) {
                ctemp = allsumtemp;
                evaluator.rotate_vector_inplace(ctemp, pow(2, k), gal_keys);
                evaluator.add_inplace(allsumtemp, ctemp);
            }
            cout << "1";
            ptemp = Cplain[l];
            cout << "2";
            evaluator.mod_switch_to_inplace(ptemp, allsumtemp.parms_id());
            cout << "3";
            evaluator.multiply_plain_inplace(allsumtemp, ptemp);
            cout << "4";
            evaluator.rescale_to_next_inplace(allsumtemp);
            cout << "5";
            allsumtemp.scale() = Y.scale();
            cout << "6";
            evaluator.add_inplace(Y, allsumtemp);
            cout << "7";
        }


            decryptor.decrypt(Y, ptemp);
            encoder.decode(ptemp, input);
            cout << "(";
            for (int i = 0; i < 20; i++)cout << input[i] << ",";
        
    }
    




    ////encode the data, a column at a time. The ith element of data corresponds to the ith feature
    //pVec dataplain;
    //dataplain.reserve(nfeatures);
    //Plaintext plain;
    //dVec input;
    //for (int i = 0; i < nfeatures; i++) {
    //    input.clear();
    //    for (int j = 0; j < Matrix.size(); j++)input.push_back(Matrix[j][i]);
    //    encoder.encode(input, scale, plain);
    //    dataplain.push_back(plain);
    //}

    //end = chrono::steady_clock::now();
    //diff = end - start;
    //cout << "Encoding time = " << chrono::duration <double, milli>(diff).count() / 1000.0 << " s \n";

    //cout << "Encoding...";
    //start = chrono::steady_clock::now();
    ////encrypt these plaintexts
    //cVec dataenc;
    //dataenc.reserve(nfeatures);
    //Ciphertext ctemp;
    //for (int i = 0; i < dataplain.size(); i++) {
    //    encryptor.encrypt(dataplain[i], ctemp);
    //    dataenc.push_back(ctemp);
    //}
    //end = chrono::steady_clock::now();
    //diff = end - start;
    //cout << "Encrypting time = " << chrono::duration <double, milli>(diff).count() / 1000.0 << " s \n";

	return 0;
}
		
	
