#include "seal/seal.h"
#include "databasetools.h"
#include "ridgeregressiontools.h"
#include "threadpool.hpp"
#include <iostream>

int PP_Nesterov_Gradient_Descent_RR(double alpha, double lambda, bool bitsize) {
    
    // Set this to nr_threads -1 : it will run even on this thread
    thread_pool::thread_pool tp(9);

    unsigned int precision{};
    unsigned int modulus_size{};
    unsigned int default_value{};
    unsigned int sentenial_value{};
    unsigned int learning_iters{};

    if(bitsize) {
        // This choice corresponds to a 40-bit precision choice
        precision        = 40;
        modulus_size     = 20;
        default_value    = 40;
        sentenial_value  = 50;
        learning_iters   = 7;
    } else {
        precision = 30;
        modulus_size = 28;
        default_value = 30;
        sentenial_value = 40;
        learning_iters = 10;
    }

    EncryptionParameters parms(scheme_type::CKKS);
    size_t poly_modulus_degree = 32768;
    vector<int> mod(modulus_size, default_value);
    mod[0] = sentenial_value;
    mod[modulus_size-1] = sentenial_value;
    parms.set_poly_modulus_degree(poly_modulus_degree);

    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, mod));
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
    dVec weights;
    /*double alpha = 0.0009;
    double lambda = 1;*/

    //define the matrix C, which is zero apart from alphas on the diagonal
    dMat C;
    for (int i = 0; i < data[0].size(); i++) {
        dVec ctemp(data[0].size(), 0.0);
        ctemp[i] += alpha;
        C.push_back(ctemp);
    }
    //define the matrix I, which is zero apart from 1 on the diagonal:
    dMat I;
    for (int i = 0; i < data[0].size(); i++) {
        dVec ctemp(data[0].size(), 0.0);
        ctemp[i] += 1;
        I.push_back(ctemp);
    }
    //set the precision
    double scale = pow(2.0, precision);

    //encode the rows of the matrix C
    pVec Cplain;
    Plaintext ptemp, ptemp1;
    for (int i = 0; i < C.size(); i++) {
        encoder.encode(C[i], scale, ptemp);
        Cplain.push_back(ptemp);
    }
    pVec Iplain;
    for (int i = 0; i < C.size(); i++) {
        encoder.encode(I[i], scale, ptemp);
        Iplain.push_back(ptemp);
    }

    //cross validation
    dMatMat cvtrain, cvtest;
    dMat cvtrainresults, cvtestresults;

    CVrandomSampling(cvtrain, cvtrainresults, cvtest, cvtestresults, data, results);
    //cout << "fold sizes are:\n1 -- " << cvtrain[0].size() << "--" << cvtrainresults[0].size() << "\n";
    //cout << "2 -- " << cvtrain[1].size() << "--" << cvtrainresults[1].size() << "\n";
    //cout << "3 -- " << cvtrain[2].size() << "--" << cvtrainresults[2].size() << "\n";
    //cout << "4 -- " << cvtrain[3].size() << "--" << cvtrainresults[3].size() << "\n";
    //cout << "5 -- " << cvtrain[4].size() << "--" << cvtrainresults[4].size() << "\n";

    //cout << "test sizes are:\n1 -- " << cvtest[0].size() << "--" << cvtestresults[0].size() << "\n";
    //cout << "2 -- " << cvtest[1].size() << "--" << cvtestresults[1].size() << "\n";
    //cout << "3 -- " << cvtest[2].size() << "--" << cvtestresults[2].size() << "\n";
    //cout << "4 -- " << cvtest[3].size() << "--" << cvtestresults[3].size() << "\n";
    //cout << "5 -- " << cvtest[4].size() << "--" << cvtestresults[4].size() << "\n";
    dVec a, b, input;
    double mean, t, T, gamma;
    int nfeatures, n;
    pVec dataplain;
    Plaintext resultsplain, lambdap, alphap, gammap, mgammap;
    cVec dataenc;
    Ciphertext ctemp, ctemp1, resultsenc, Y, allsumtemp, Ytemp;

    input.push_back(0);
    for (int i = 1; i < cvtrain[0][0].size(); i++)input.push_back(lambda);
    encoder.encode(lambda, scale, lambdap);
    input.clear();
    encoder.encode(alpha, scale, alphap);

    double avgr2 = 0;

    for (int j = 0; j < 5; j++) {
        a.clear();
        b.clear();
        dataplain.clear();
        dataenc.clear();
        mean = 0;
        t = 1;
        cout << "starting fold " << j + 1 << "\n";
        //start by centering the results, and fitting the scale to the data
        center(cvtrainresults[j], mean);
        scale_fit(cvtrain[j], a, b);

        //extract dimensions
        nfeatures = cvtrain[j][0].size();
        n = cvtrain[j].size();
        //split the data into columns, leaving the first as it's a column of 1s. This would happen client side


        cout << "Encoding...";
        start = chrono::steady_clock::now();
        for (int i = 1; i < nfeatures; i++) {
            input.clear();
            for (int k = 0; k < n; k++)input.push_back(cvtrain[j][k][i]);
            encoder.encode(input, scale, ptemp);
            dataplain.push_back(ptemp);
        }

        encoder.encode(cvtrainresults[j], scale, resultsplain);
        end = chrono::steady_clock::now();
        diff = end - start;
        cout << "Encoding time = " << chrono::duration <double, milli>(diff).count() / 1000.0 << " s \n";

        cout << "Encrypting...";
        start = chrono::steady_clock::now();

        for (int i = 0; i < dataplain.size(); i++) {
            encryptor.encrypt(dataplain[i], ctemp);
            dataenc.push_back(ctemp);
        }

        encryptor.encrypt(resultsplain, resultsenc);
        end = chrono::steady_clock::now();
        diff = end - start;
        cout << "Encrypting time = " << chrono::duration <double, milli>(diff).count() / 1000.0 << " s \n";
        auto trainstart = chrono::steady_clock::now();
        cout << "Forming Y...";

        //first we create a vector with Y1 in the first slot --
        //multiply X1 and the results vector
        evaluator.multiply(dataenc[0], resultsenc, allsumtemp);
        evaluator.relinearize_inplace(allsumtemp, relin_keys);
        evaluator.rescale_to_next_inplace(allsumtemp);
        //allsum 
        for (int k = 0; k < log2(slot_count); k++) {
            ctemp = allsumtemp;
            evaluator.rotate_vector_inplace(ctemp, pow(2, k), gal_keys);
            evaluator.add_inplace(allsumtemp, ctemp);
        }

        //multiply by the 1st row of alpha C, only retaining sum in 1st slot
        ptemp = Cplain[1];
        evaluator.mod_switch_to_inplace(ptemp, allsumtemp.parms_id());
        evaluator.multiply_plain(allsumtemp, ptemp, Y);
        evaluator.rescale_to_next_inplace(Y);
        //now we add all the other entries
        
        std::vector<Ciphertext> allsums(nfeatures-2);
        std::mutex allsum_mutex;
        for (int l = 2; l < nfeatures; l++) {
            tp.push([l, &dataenc,&resultsenc, &evaluator, &relin_keys, &gal_keys, slot_count, &Cplain, &allsums, &allsum_mutex]() {
                Ciphertext ctemp;
                //multiply Xl and the results vector
                evaluator.multiply(dataenc[1. * l - 1], resultsenc, ctemp);
                evaluator.relinearize_inplace(ctemp, relin_keys);
                evaluator.rescale_to_next_inplace(ctemp);
                
                Ciphertext allsumtemp = ctemp;
                //allsum 
                for (int k = 0; k < log2(slot_count); k++) {
                    ctemp = allsumtemp;
                    evaluator.rotate_vector_inplace(ctemp, pow(2, k), gal_keys);
                    evaluator.add_inplace(allsumtemp, ctemp);
                }

                //multiply by the lth row of alpha C, only retaining sum in lth slot
                Plaintext ptemp = Cplain[l];
                evaluator.mod_switch_to_inplace(ptemp, allsumtemp.parms_id());
                evaluator.multiply_plain_inplace(allsumtemp, ptemp);
                evaluator.rescale_to_next_inplace(allsumtemp);
                std::lock_guard<std::mutex> lock(allsum_mutex);
                allsums[l-2] = allsumtemp;
            });
        }
        
        tp.wait_work();
        
        for(const auto& elem : allsums) {
            Y.scale() = elem.scale();
            evaluator.add_inplace(Y, elem);
        }

        //first iteration is Beta = Y,v=(1-gamma)Y=Y. We set this now: 
        Ciphertext Beta = Y;
        Ciphertext v = Beta;
        //and change t:
        t = (1. + sqrt(1. + 4 * t * t)) / 2.;

        Ciphertext Betaold;
        cout << "Forming the Matrix M...\n";
        cVec M(nfeatures);
        start = chrono::steady_clock::now();
        //the first entry of M should have first entry alpha*n. We can do this manually.
        input.clear();
        input.push_back(n);
        encoder.encode(n, scale, ptemp);
        encryptor.encrypt(ptemp, ctemp);
        evaluator.multiply_plain_inplace(ctemp, Cplain[0]);
        evaluator.rescale_to_next_inplace(ctemp);
        M[0] = ctemp;
        
        std::vector<Ciphertext> Xtemps(nfeatures-1);
        std::vector<Plaintext>  Ptemps(nfeatures-1);
        std::mutex X_mutex;
        std::mutex M_mutex;

        //now we go feature by feature, forming the rows of M
        for (int i = 1; i < nfeatures; i++) {
            tp.push([i, &dataenc, slot_count, &evaluator, &gal_keys, &relin_keys, &Cplain, &M, &Xtemps, &Ptemps, &X_mutex, &M_mutex](){
                Ciphertext Xtemp = dataenc[1. * i - 1];
                //first we allsum Xi, to find Mi0 = M0i
                Ciphertext allsumtemp = Xtemp;
                for (int k = 0; k < log2(slot_count); k++) {
                    Ciphertext ctemp = allsumtemp;
                    evaluator.rotate_vector_inplace(ctemp, pow(2, k), gal_keys);
                    evaluator.add_inplace(allsumtemp, ctemp);
                }
            
                //we need two copies of allsum temp -- one for Mi and one for M0.
                //first we eliminate all but the ith entry and add to M0:
                Ciphertext ctemp;
                Plaintext  ptemp;
                evaluator.multiply_plain(allsumtemp, Cplain[i], ctemp);
                evaluator.rescale_to_next_inplace(ctemp);
                {
                    std::lock_guard<std::mutex> lock(M_mutex);
                    evaluator.add_inplace(M[0], ctemp);
                }
            
                //now we eliminate all but the 0th entry
                evaluator.multiply_plain_inplace(allsumtemp, Cplain[0]);
                evaluator.rescale_to_next_inplace(allsumtemp);
                evaluator.mod_switch_to_next_inplace(allsumtemp);
                //form Mii:

                evaluator.square(Xtemp, ctemp);
                evaluator.relinearize_inplace(ctemp, relin_keys);
                evaluator.rescale_to_next_inplace(ctemp);
                for (int l = 0; l < log2(slot_count); l++) {
                    Ciphertext ctemp1 = ctemp;
                    evaluator.rotate_vector_inplace(ctemp1, pow(2, l), gal_keys);
                    evaluator.add_inplace(ctemp, ctemp1);
                }

                //switch down C[i]&multiply, eliminating all but the ith entry
                evaluator.mod_switch_to_next(Cplain[i], ptemp);
                evaluator.multiply_plain_inplace(ctemp, ptemp);
                evaluator.rescale_to_next_inplace(ctemp);
                //add to allsumtemp, generating a ciphertext Mi0,0,...,0,Mii,0,...,0
                allsumtemp.scale() = ctemp.scale();
                evaluator.add_inplace(allsumtemp, ctemp);
                {
                    std::lock_guard<std::mutex> lock(M_mutex);
                    M[i] = allsumtemp;
                }
                {
                    std::lock_guard<std::mutex> lock(X_mutex);
                    Xtemps[i-1] = Xtemp;
                    Ptemps[i-1] = ptemp;
                }
            });
        }

        tp.wait_work();
        
        for(int i = 1; i < nfeatures;i++) {
            auto Xtemp = Xtemps[i-1];
            auto ptemp = Ptemps[i-1];

            std::vector<Ciphertext> M_i_list(i-1);
            std::mutex M_i_mutex;
            //now we form Mik for k < i. This means we don't need to calculate Mik and Mki separately/twice
            for (int k = 1; k < i; k++) {
               tp.push([&evaluator, &dataenc, &Xtemp, &relin_keys, &gal_keys, i, k, slot_count, &M ,&ptemp, &Cplain, &M_i_list, &M_mutex, &M_i_mutex](){
                
                    //first multiply Xi by Xk
                    Ciphertext allsumtemp;
                    Ciphertext ctemp;

                    evaluator.multiply(dataenc[1. * k - 1], Xtemp, allsumtemp);
                    evaluator.relinearize_inplace(allsumtemp, relin_keys);
                    evaluator.rescale_to_next_inplace(allsumtemp);
                    
                    //now we need to allsum allsumtemp.
                    for (int l = 0; l < log2(slot_count); l++) {
                        ctemp = allsumtemp;
                        evaluator.rotate_vector_inplace(ctemp, pow(2, l), gal_keys);
                        evaluator.add_inplace(allsumtemp, ctemp);
                    }

                    //we need two copies of allsumtemp -- one for M[i] in the kth slot and one for
                    //M[k] in the ith slot...
                    //switch down a copy of C[k]             
                    Plaintext ptemp1;

                    evaluator.mod_switch_to_next(Cplain[k], ptemp1);
                    //M[k] in the ith slot                
                    evaluator.multiply_plain(allsumtemp, ptemp, ctemp);
                    evaluator.rescale_to_next_inplace(ctemp);
                
                    {
                        std::lock_guard<std::mutex> M_lock(M_mutex);
                        evaluator.add_inplace(M[k], ctemp);
                    }
                
                    //M[i] in the kth slot                
                    evaluator.multiply_plain_inplace(allsumtemp, ptemp1);
                    evaluator.rescale_to_next_inplace(allsumtemp);
                    {
                        std::lock_guard<std::mutex> M_i_lock(M_i_mutex);
                        M_i_list[k-1] = allsumtemp;
                    }
            });
        }
        
        tp.wait_work();

        for(const auto& allsumtemp : M_i_list) {
            evaluator.add_inplace(M[i], allsumtemp);
        }
    }


        for (int i = 0; i < nfeatures; i++) {
            evaluator.negate_inplace(M[i]);
        }
        end = chrono::steady_clock::now();
        diff = end - start;
        cout << "done. Time = " << chrono::duration <double, milli>(diff).count() / 1000.0 << " s \n";
        std::cout << "Starting iterations" << std::endl;
        allsums.clear();
        allsums.resize(nfeatures);
        //iterations:
        for (int k = 2; k < learning_iters; k++) {
            cout << "starting iteration " << k << "\n";
            T = (1. + sqrt(1. + 4 * t * t)) / 2.;
            gamma = (1. - t) / T;
            t = T;
            encoder.encode(gamma, scale, gammap);
            gamma = 1 - gamma;
            encoder.encode(gamma, scale, mgammap);
            //store a copy of Beta
            Betaold = Beta;
            Ytemp = Y;
            //go feature by feature
            
            for (int i = 0; i < nfeatures; i++) {
                tp.push([&M, &allsum_mutex, slot_count, &allsums, &relin_keys, &gal_keys, &evaluator, &Iplain, i, &v]() {
                    auto M_local = M[i];
                    Ciphertext ctemp;
                    Ciphertext allsumtemp;
                    

                    //first multiply beta by Mi:
                    evaluator.mod_switch_to_inplace(M_local, v.parms_id());
                    evaluator.multiply(M_local, v, allsumtemp);
                    evaluator.relinearize_inplace(allsumtemp, relin_keys);
                    evaluator.rescale_to_next_inplace(allsumtemp);
                    
                    //now all sum
                    for (int l = 0; l < log2(slot_count); l++) {
                        ctemp = allsumtemp;
                        evaluator.rotate_vector_inplace(ctemp, pow(2, l), gal_keys);
                        evaluator.add_inplace(allsumtemp, ctemp);
                    }

                    Plaintext ptemp;
                    //delete all but the ith entry by multiplying with Ii
                    evaluator.mod_switch_to(Iplain[i], allsumtemp.parms_id(), ptemp);
                    evaluator.multiply_plain_inplace(allsumtemp, ptemp);
                    evaluator.rescale_to_next_inplace(allsumtemp);
                    std::lock_guard<std::mutex> lock(allsum_mutex);
                    allsums[i] = allsumtemp;
                    });
            }

            tp.wait_work();
            for(const auto& allsumtemp : allsums) {
                    //adjust scale of Y (should only affect i=0) and add:
                    evaluator.mod_switch_to_inplace(Ytemp, allsumtemp.parms_id());
                    
                    Ytemp.scale() = allsumtemp.scale();
                    evaluator.add_inplace(Ytemp, allsumtemp);
            }

            //create (1-lambda*alpha)v:
            ctemp = v;
            evaluator.mod_switch_to(lambdap, v.parms_id(), ptemp);
            evaluator.multiply_plain_inplace(v, ptemp);
            evaluator.rescale_to_next_inplace(v);
            evaluator.mod_switch_to(alphap, v.parms_id(), ptemp);
            evaluator.multiply_plain_inplace(v, ptemp);
            evaluator.rescale_to_next_inplace(v);
            evaluator.negate_inplace(v);
            evaluator.mod_switch_to_next_inplace(ctemp);
            evaluator.mod_switch_to_next_inplace(ctemp);
            ctemp.scale() = v.scale();
            evaluator.add_inplace(v, ctemp);
            //and add to Y, forming Beta         
            Ytemp.scale() = v.scale();
            evaluator.add(v, Ytemp, Beta);
            //now update v = (1-gamma)beta+gamma*betaold

            evaluator.mod_switch_to_inplace(mgammap, Beta.parms_id());

            evaluator.multiply_plain(Beta, mgammap, ctemp);
            evaluator.rescale_to_next_inplace(ctemp);
            evaluator.mod_switch_to_inplace(gammap, Betaold.parms_id());
            evaluator.multiply_plain_inplace(Betaold, gammap);
            evaluator.rescale_to_next_inplace(Betaold);
            Betaold.scale() = ctemp.scale();
            evaluator.mod_switch_to_inplace(Betaold, ctemp.parms_id());
            evaluator.add(ctemp, Betaold, v);
            weights.clear();
            decryptor.decrypt(v, ptemp);
            encoder.decode(ptemp, input);
            for (int i = 0; i < nfeatures; i++)weights.push_back(input[i]);

            cout << "iteration " << k << " r^2: " << Rsquared(cvtrain[j], cvtrainresults[j], weights) << "\n";
            /*cout << "iteration " << k << " CV r^2: " << Rsquared(cvtest[j], cvtestresults[j], weights) << "\n";*/
        }
        
        auto trainend = chrono::steady_clock::now();
        diff = trainend - trainstart;
        cout << "fold " << j + 1 << "training done. Time = " << chrono::duration <double, milli>(diff).count() / 1000.0 << " s \n";
        weights.clear();
        decryptor.decrypt(v, ptemp);
        encoder.decode(ptemp, input);
        scale_columns(cvtest[j], a, b);
        shift_results(cvtestresults[j], mean);
        for (int i = 0; i < nfeatures; i++)weights.push_back(input[i]);
        cout << "fold " << j + 1 << " final r^2: " << Rsquared(cvtrain[j], cvtrainresults[j], weights) << "%\n";
        cout << "fold " << j + 1 << " cross validation r^2: " << Rsquared(cvtest[j], cvtestresults[j], weights) << "%\n";
        avgr2 += Rsquared(cvtest[j], cvtestresults[j], weights);
        std::cout << std::endl;
    }
    cout << "Average Cross Validation r^2: " << avgr2 / 5 << "%";
    return 0;
}
