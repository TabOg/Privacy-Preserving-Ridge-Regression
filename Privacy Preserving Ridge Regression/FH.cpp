#include "seal/seal.h"
#include "databasetools.h"
#include "ridgeregressiontools.h"
#include <iostream>
#include "threadpool.hpp"

int PP_Fixed_Hessian_RR(double lambda) {    

    thread_pool::thread_pool tp(10);

    EncryptionParameters parms(scheme_type::CKKS);
    size_t poly_modulus_degree = 32768;
    vector<int> mod;
    mod.push_back(50);
    for (int i = 0; i < 18; i++)mod.push_back(40);
    mod.push_back(50);
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
    

    
    //define the matrix I, which is zero apart from 1 on the diagonal:
    dMat I;
    for (int i = 0; i < data[0].size(); i++) {
        dVec ctemp(data[0].size(), 0.0);
        ctemp[i] += 1;
        I.push_back(ctemp);
    }


    //set 40 bits of precision
    double scale = pow(2.0, 40);
    Plaintext ptemp, ptemp1;
    //encode the rows of the matrix I
    pVec Iplain;
    for (int i = 0; i < I.size(); i++) {
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
    double mean;
    int nfeatures, n;
    pVec dataplain;
    Plaintext resultsplain, lambdap, alphap;
    cVec dataenc, H;
    Ciphertext ctemp, ctemp1, resultsenc, Y, allsumtemp, Ytemp, HVec;
    input.push_back(0);
    for (int i = 1; i < cvtrain[0][0].size(); i++)input.push_back(lambda);
    encoder.encode(lambda, scale, lambdap);
    input.clear();

    double avgr2 = 0;

    for (int j = 0; j < 5; j++) {
        a.clear();
        b.clear();
        dataplain.clear();
        dataenc.clear();
        H.clear();
        mean = 0;
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

        //multiply by the 1st row of I, only retaining sum in 1st slot
        ptemp = Iplain[1];
        evaluator.mod_switch_to_inplace(ptemp, allsumtemp.parms_id());
        evaluator.multiply_plain(allsumtemp, ptemp, Y);
        evaluator.rescale_to_next_inplace(Y);
        //now we add all the other entries
        std::mutex mutex;
        std::vector<Ciphertext> allsums(nfeatures-2);

        for (int l = 2; l < nfeatures; l++) {
            tp.push([&allsums, &mutex, &evaluator, l, &resultsenc, &relin_keys, &gal_keys, &dataenc, slot_count, &Iplain](){
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

                //multiply by the lth row of I, only retaining sum in lth slot
                Plaintext ptemp = Iplain[l];
                evaluator.mod_switch_to_inplace(ptemp, allsumtemp.parms_id());
                evaluator.multiply_plain_inplace(allsumtemp, ptemp);
                evaluator.rescale_to_next_inplace(allsumtemp);
                std::lock_guard<std::mutex> lock(mutex);
                allsums[l-2] = allsumtemp;
            });
        }

        tp.wait_work();
        // N.B In future iterations/releases this could be replaced 
        // by an addition tree if it is too slow (depending on the size of the allsums vector)
        for(const auto& allsumelem : allsums) {
            evaluator.add_inplace(Y, allsumelem);
        }

        cout << "Forming the Matrix M...\n";
        cVec M;

        //the first entry of M should have first entry n. We can do this manually.
        encoder.encode(n, scale, ptemp);
        encryptor.encrypt(ptemp, ctemp);
        //we need a copy of this for H[0]:
        H.resize(nfeatures);
        H[0] = ctemp;
        //but it will need to be one level lower
        evaluator.mod_switch_to_next_inplace(H[0]);
        //cout << "after feature 0:\n";
        //for (int i = 0; i < H.size(); i++) {
        //    decryptor.decrypt(H[i], ptemp);
        //    encoder.decode(ptemp, input);
        //    cout << "(";
        //    for (int k = 0; k < 15; k++)cout << input[k] << ",";
        //    cout << "\n";
        //}
        
        //back to M: retain 0th entry of (n,n,...,n)
        evaluator.multiply_plain_inplace(ctemp, Iplain[0]);
        evaluator.rescale_to_next_inplace(ctemp);
        M.resize(nfeatures);
        M[0] = ctemp;
        //now we go feature by feature, forming the rows of M and H
        // This loop is parallelised using the default strategy of pushing elements into the threadpool
        // and simply waiting for them to finish, writing back any underlying changes under the control of the mutex. 
        
        // To avoid congestion on a single mutex, here we just use a mutex for each H, M and I.
        std::mutex H_mutex;
        std::mutex M_mutex;
        std::mutex I_mutex; 
        std::vector<Ciphertext> Xtemps(nfeatures-1);

        for (int i = 1; i < nfeatures; i++) {
            tp.push([i, &dataenc, &Xtemps, slot_count, nfeatures, &gal_keys, &relin_keys, &evaluator, &H,&M, &Iplain, &H_mutex, &M_mutex, &I_mutex](){
                Ciphertext Xtemp = dataenc[1. * i - 1];
                //first we allsum Xi, to find Mi0 = M0i
                Ciphertext allsumtemp = Xtemp;
                Ciphertext ctemp;

                for (int k = 0; k < log2(slot_count); k++) {
                    ctemp = allsumtemp;
                    evaluator.rotate_vector_inplace(ctemp, pow(2, k), gal_keys);
                    evaluator.add_inplace(allsumtemp, ctemp);
                }
                //for H, we need this added to H[0] & forming the basis for H[i], both one level down:
                evaluator.mod_switch_to_next(allsumtemp, ctemp);
                
                // Note: because this lock occurs early in the loop
                // we can't wait for std::lock_guard to fall out of scope (it would deadlock H_mutex for all other threads until this loops again!)
                // So here instead we use unique_lock, which allows us better control
                // over when the lock is released
                std::unique_lock<std::mutex> H_lock(H_mutex);
                H_lock.lock();
                evaluator.add_inplace(H[0], ctemp);
                H_lock.unlock();

                // Rather than staying in a lock longer than we need to, we instead
                // store a locak copy of H[i] and use that: we will write back at the end.
                auto H_temp = allsumtemp;

                //cout << "adding M" << i << "0:\n";
                //for (int i = 0; i < 1; i++) {
                //    decryptor.decrypt(H[i], ptemp);
                //    encoder.decode(ptemp, input);
                //    cout << "(";
                //    for (int k = 0; k < 15; k++)cout << input[k] << ",";
                //    cout << "\n";
                //}
                //For M, we need two copies of allsum temp -- one for Mi and one for M0.
                //first we eliminate all but the ith entry and add to M0:
                //This does not need a lock because it's a read on Iplain!
                evaluator.multiply_plain(allsumtemp, Iplain[i], ctemp);
                evaluator.rescale_to_next_inplace(ctemp);

                std::unique_lock<std::mutex> M_lock;
                M_lock.lock();
                evaluator.add_inplace(M[0], ctemp);
                M_lock.unlock();
                //now we eliminate all but the 0th entry
                evaluator.multiply_plain_inplace(allsumtemp, Iplain[0]);
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

                //For H, we need this added to H[i]. H[i] is already at the correct level, just needs 
                //a modified scale:
                H_temp.scale() = ctemp.scale();
                evaluator.add_inplace(H_temp,ctemp);
                
                //For M, switch down I[i]&multiply, eliminating all but the ith entry
                Plaintext ptemp;
                std::unique_lock<std::mutex> I_lock(I_mutex);
                I_lock.lock();
                evaluator.mod_switch_to_next(Iplain[i], ptemp);
                I_lock.unlock();

                evaluator.multiply_plain_inplace(ctemp, ptemp);
                evaluator.rescale_to_next_inplace(ctemp);
                //add to allsumtemp, generating a ciphertext Mi0,0,...,0,Mii,0,...,0
                allsumtemp.scale() = ctemp.scale();
                evaluator.add_inplace(allsumtemp, ctemp);
                
                // Note: Here we need this lock, as otherwise we won't be synchronised for other loops.
                M_lock.lock();
                M[i] = allsumtemp;
                M_lock.unlock();

                H_lock.lock();
                H[i] = H_temp;
                Xtemps[i-1] = Xtemp;
                H_lock.unlock();
            });
        }

        tp.wait_work();

        // Need to add learning rate to H[d]
        Plaintext ptemp;
        evaluator.mod_switch_to(lambdap, H[nfeatures-1].parms_id(), ptemp);
        ptemp.scale() = H[nfeatures-1].scale();
        evaluator.add_plain_inplace(H[nfeatures-1], ptemp);
        
        for(int i = 1; i < nfeatures;i++) {
        //now we form Mik for k < i. This means we don't need to calculate Mik and Mki separately/twice
            for (int k = 1; k < i; k++) {
                //first multiply Xi by Xk
                Ciphertext Xtemp = Xtemps[i-1];
                evaluator.multiply(dataenc[1. * k - 1], Xtemp, allsumtemp);
                evaluator.relinearize_inplace(allsumtemp, relin_keys);
                evaluator.rescale_to_next_inplace(allsumtemp);
                //now we need to allsum allsumtemp.
                for (int l = 0; l < log2(slot_count); l++) {
                    ctemp = allsumtemp;
                    evaluator.rotate_vector_inplace(ctemp, pow(2, l), gal_keys);
                    evaluator.add_inplace(allsumtemp, ctemp);
                }

                //For H, we need this added to both H[i] and H[k]. Both should be at the right level
                //and scale:
                 
                evaluator.add_inplace(H[i], allsumtemp);
                evaluator.add_inplace(H[k], allsumtemp);
                
                //For M, we need two copies of allsumtemp -- one for M[i] in the kth slot and one for
                //M[k] in the ith slot...
                //switch down a copy of I[k]                
                evaluator.mod_switch_to_next(Iplain[k], ptemp1);
                //M[k] in the ith slot                
                evaluator.multiply_plain(allsumtemp, ptemp, ctemp);
                evaluator.rescale_to_next_inplace(ctemp);
                evaluator.add_inplace(M[k], ctemp);
                //M[i] in the kth slot                
                evaluator.multiply_plain_inplace(allsumtemp, ptemp1);
                evaluator.rescale_to_next_inplace(allsumtemp);
                evaluator.add_inplace(M[i], allsumtemp);
                //we need to add lambda to H[i] for i >0:
                if (i == nfeatures - 1) {
                    evaluator.mod_switch_to(lambdap, H[k].parms_id(), ptemp1);
                    ptemp1.scale() = H[k].scale();
                    evaluator.add_plain_inplace(H[k], ptemp1);
                }
            }
        }
        
        

            //cout << "after feature " << i << ": \n";
            //for (int i = 0; i < H.size(); i++) {
            //    decryptor.decrypt(H[i], ptemp);
            //    encoder.decode(ptemp, input);
            //    cout << "(";
            //    for (int k = 0; k < 15; k++)cout << input[k] << ",";
            //    cout << "\n";
            //}
            //cout << "\n";
        
        
    for (int i = 0; i < nfeatures; i++) {
            evaluator.negate_inplace(M[i]);
        }
        //for (int i = 0; i < H.size(); i++) {
        //    decryptor.decrypt(H[i], ptemp);
        //    encoder.decode(ptemp, input);
        //    cout << "(";
        //    for (int k = 0; k < 15; k++)cout << input[k];
        //    cout << "\n";
        //}
        cout << "\n";
        dataenc.clear();
        //time to calculate 1/H(i,i) -- first we need our starting point, T1 + T2D
        cout << "Calculating 1/H(i)...\n";
        start = chrono::steady_clock::now();
        Plaintext P_T1, P_T2;
        double t1, t2;
        t1 = T1(1. * nfeatures * n);
        t2 = T2(1. * nfeatures * n);
        encoder.encode(t1, scale, P_T1);
        encoder.encode(t2, scale, P_T2);

        //T2 needs to be multiplied by something 1 level down, so:
        evaluator.mod_switch_to_next_inplace(P_T2);

        //T1 will need to be added to something 2 levels down, so:
        evaluator.mod_switch_to_next_inplace(P_T1);
        evaluator.mod_switch_to_next_inplace(P_T1);
        
        
        for (int i = 0; i < nfeatures; i++) {
            tp.push([i, &H, &H_mutex, &evaluator, &P_T1, &P_T2, &relin_keys](){
                //negate and store a copy of H(i,i) for later:
                Ciphertext H_temp = H[i];
                Ciphertext ctemp1 = H[i];
                Ciphertext ctemp;

                evaluator.negate_inplace(ctemp1);

                //find v0 = T1 + T2H(i,i):
                evaluator.multiply_plain_inplace(H_temp, P_T2);
                evaluator.rescale_to_next_inplace(H_temp);
                P_T1.scale() = H_temp.scale();
                evaluator.add_plain_inplace(H_temp, P_T1);

                //now iterate Newton Raphson: each update is 2v - H[i,i]v^2
                for (int j = 0; j < 2; j++) {

                    //first double and store the result
                    evaluator.add(H_temp, H_temp, ctemp);

                    //now square the current value, relin and rescale

                    evaluator.square_inplace(H_temp);
                    evaluator.relinearize_inplace(H_temp, relin_keys);
                    evaluator.rescale_to_next_inplace(H_temp);

                    //now mod switch down our stored value Htemp, and multiply
                    evaluator.mod_switch_to_inplace(ctemp1, H_temp.parms_id());
                    evaluator.multiply_inplace(H_temp, ctemp1);
                    evaluator.relinearize_inplace(H_temp, relin_keys);
                    evaluator.rescale_to_next_inplace(H_temp);

                    //modify scale of ctemp = 2v, mod switch down, and add
                    ctemp.scale() = H_temp.scale();
                    evaluator.mod_switch_to_inplace(ctemp, H_temp.parms_id());
                    evaluator.add_inplace(H_temp, ctemp);
                }

                std::lock_guard<std::mutex> lock(H_mutex);
                H[i] = H_temp;
            });
        }

        tp.wait_work();

        //need to eliminate all but the ith entry from H[i]:
        evaluator.mod_switch_to(Iplain[0], H[0].parms_id(), ptemp);
        evaluator.multiply_plain_inplace(H[0], ptemp);
        evaluator.rescale_to_next_inplace(H[0]);
        
        //and we need to form HVec = (1/H00,1/H11,...,1/Hdd)
        HVec = H[0];
        for (int i = 1; i < H.size(); i++) {
            evaluator.mod_switch_to(Iplain[i], H[i].parms_id(), ptemp);
            evaluator.multiply_plain_inplace(H[i], ptemp);
            evaluator.rescale_to_next_inplace(H[i]);
            HVec.scale() = H[i].scale();
            evaluator.add_inplace(HVec, H[i]);
        
        }

        //decryptor.decrypt(HVec, ptemp);
        //encoder.decode(ptemp, input);
        //cout << "(";
        //for (int i = 0; i < 15; i++)cout << input[i] << ",";
        //cout << "\n";
        //first iteration is beta = HVecY. We set this now:
        evaluator.mod_switch_to_inplace(Y, HVec.parms_id());
        evaluator.multiply_inplace(Y, HVec);
        evaluator.relinearize_inplace(Y, relin_keys);
        evaluator.rescale_to_next_inplace(Y);
        Ciphertext Beta = Y;
        //we will also only need lambda in (0,lambda/H11,lambda/H22,...,lambda/Hdd)
        evaluator.mod_switch_to(lambdap, HVec.parms_id(), ptemp);
        evaluator.multiply_plain(HVec, ptemp, ctemp1);
        evaluator.rescale_to_next_inplace(ctemp1);

        std::vector<Ciphertext> betas(nfeatures);
        
        allsums.clear();
        allsums.resize(nfeatures);

        //iterations:
        for (int k = 2; k < 7; k++) {
            cout << "starting iteration " << k << "\n";
            Ytemp = Y;
    
            //go feature by feature
            for (int i = 0; i < nfeatures; i++) {
                tp.push([i, &M, &M_mutex, &H, &H_mutex, &Beta, &relin_keys, &gal_keys, &mutex, &evaluator, &allsums, slot_count](){
                        Ciphertext mtemp = M[i];
                        Ciphertext allsumtemp;
                        //first multiply beta by Mi:
                        evaluator.mod_switch_to_inplace(mtemp, Beta.parms_id());
                        evaluator.multiply(mtemp, Beta, allsumtemp);
                        
                        evaluator.relinearize_inplace(allsumtemp, relin_keys);
                        evaluator.rescale_to_next_inplace(allsumtemp);
                        //now all sum
                        for (int l = 0; l < log2(slot_count); l++) {
                            Ciphertext ctemp = allsumtemp;
                            evaluator.rotate_vector_inplace(ctemp, pow(2, l), gal_keys);
                            evaluator.add_inplace(allsumtemp, ctemp);
                        }

                        Ciphertext htemp = H[i];
                        //delete all but the ith entry & multiply by Hi
                        evaluator.mod_switch_to_inplace(htemp, allsumtemp.parms_id());
                        evaluator.multiply_inplace(allsumtemp, htemp);
                        evaluator.relinearize_inplace(allsumtemp, relin_keys);
                        evaluator.rescale_to_next_inplace(allsumtemp);
                
                        std::lock_guard<std::mutex> a_lock(mutex);
                        allsums[i] = allsumtemp;
                        
                        std::lock_guard<std::mutex> M_lock(M_mutex);
                        M[i] = mtemp;

                        std::lock_guard<std::mutex> H_lock(H_mutex);
                        H[i] = htemp;
                    });
            }

            tp.wait_work();
            for(const auto& allsumtemp : allsums) {
                    //adjust scale & modulus of Y (should only affect i=0) and add:
                    evaluator.mod_switch_to_inplace(Ytemp, allsumtemp.parms_id());
                    Ytemp.scale() = allsumtemp.scale();
                    evaluator.add_inplace(Ytemp, allsumtemp);
            }

            //create (1-lambda*H)beta:
            ctemp = Beta;
            // ctemp1 should contain Hvec * ptemp in it
            evaluator.mod_switch_to_inplace(ctemp1, Beta.parms_id());
            evaluator.multiply_inplace(Beta, ctemp1);
            evaluator.relinearize_inplace(Beta, relin_keys);
            evaluator.rescale_to_next_inplace(Beta);
            evaluator.negate_inplace(Beta);
            evaluator.mod_switch_to_next_inplace(ctemp);
            ctemp.scale() = Beta.scale();
            evaluator.add_inplace(Beta, ctemp);
            //and add to Y
            evaluator.mod_switch_to_inplace(Beta, Ytemp.parms_id());
            Ytemp.scale() = Beta.scale();
            evaluator.add_inplace(Beta, Ytemp);
            weights.clear();
            decryptor.decrypt(Beta, ptemp);
            encoder.decode(ptemp, input);
            for (int i = 0; i < nfeatures; i++)weights.push_back(input[i]);

            cout << "iteration " << k << " r^2: " << Rsquared(cvtrain[j], cvtrainresults[j], weights) << "\n";
            /*cout << "iteration " << k << " CV r^2: " << Rsquared(cvtest[j], cvtestresults[j], weights) << "\n";*/
        }

        auto trainend = chrono::steady_clock::now();
        diff = trainend - trainstart;
        cout << "fold " << j + 1 << "training done. Time = " << chrono::duration <double, milli>(diff).count() / 1000.0 << " s \n";
        weights.clear();
        decryptor.decrypt(Beta, ptemp);
        encoder.decode(ptemp, input);
        scale_columns(cvtest[j], a, b);
        shift_results(cvtestresults[j], mean);
        for (int i = 0; i < nfeatures; i++)weights.push_back(input[i]);
        cout << "fold " << j + 1 << " final r^2: " << Rsquared(cvtrain[j], cvtrainresults[j], weights) << "%\n";
        cout << "fold " << j + 1 << " cross validation r^2: " << Rsquared(cvtest[j], cvtestresults[j], weights) << "%\n";
        avgr2 += Rsquared(cvtest[j], cvtestresults[j], weights);
        cout << "Modulus chain index for Beta "
            << context->get_context_data(Beta.parms_id())->chain_index() << endl;
    }
    cout << "Average Cross Validation r^2: " << avgr2 / 5 << "%";

    return 0;
}
