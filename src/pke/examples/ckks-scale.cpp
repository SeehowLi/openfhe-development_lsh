/*
  Simple examples for CKKS scale and modulus
 */

#define PROFILE

#include "openfhe.h"

using namespace lbcrypto;

int main() {

    uint32_t multDepth = 4;


    
    uint32_t scaleModSize = 52;
    uint32_t firstmod = 59;


    uint32_t batchSize = 8;

    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetScalingTechnique(FIXEDAUTO);
    parameters.SetScalingModSize(scaleModSize);
    parameters.SetFirstModSize(firstmod);
    parameters.SetBatchSize(batchSize);
    std::cout << "Scale Mode is: " << parameters.GetScalingTechnique() << std::endl;

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    auto algo = cc->GetScheme();

    // Enable the features that you wish to use
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    std::cout << "CKKS scheme is using ring dimension " << cc->GetRingDimension() << std::endl << std::endl;

    auto keys = cc->KeyGen();

    cc->EvalMultKeyGen(keys.secretKey);

    cc->EvalRotateKeyGen(keys.secretKey, {1, -2});


    // Inputs
    std::vector<double> x1 = {0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> x2 = {5.0, 4.0, 3.0, 2.0, 1.0, 0.75, 0.5, 0.25};

    // Encoding as plaintexts
    Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(x1);
    Plaintext ptxt2 = cc->MakeCKKSPackedPlaintext(x2);

    std::cout << "Input x1: " << ptxt1 << std::endl;
    std::cout << "Input x2: " << ptxt2 << std::endl;

    // Encrypt the encoded vectors
    auto c1 = cc->Encrypt(keys.publicKey, ptxt1);
    auto c2 = cc->Encrypt(keys.publicKey, ptxt2);
    std::cout << "After encrypt level:" << c1->GetLevel() << " Encoding scale:" << log2(c1->GetScalingFactor()) << std::endl;
    std::cout << "Initial q: " << c1->GetElements()[0].GetModulus().GetMSB() << std::endl;
    std::cout << "scale_deg: " << c1->GetNoiseScaleDeg() << std::endl;
    
    auto paramsQ = cc->GetElementParams()->GetParams();
    std::cout << "\nModuli in Q:" << std::endl;
    for (uint32_t i = 0; i < paramsQ.size(); i++) {
      // q0 is a bit larger because its default size is 60 bits.
      // One can change this by supplying the firstModSize argument
      // in genCryptoContextCKKS.
      std::cout << "q" << i << ": " << log2(paramsQ[i]->GetModulus().ConvertToDouble()) << std::endl;
    }
    // auto paramsQP = cryptoParamsCKKS->GetParamsQP();
    // std::cout << "Moduli in P: " << std::endl;
    // BigInteger P = BigInteger(1);
    // for (uint32_t i = 0; i < paramsQP->GetParams().size(); i++) {
    //   if (i > paramsQ.size()) {
    //     P = P * BigInteger(paramsQP->GetParams()[i]->GetModulus());
    //     std::cout << "p" << i - paramsQ.size() << ": "
    //               << paramsQP->GetParams()[i]->GetModulus() << std::endl;
    //   }
    // }

    // Homomorphic addition
    auto cAdd = cc->EvalAdd(c1, c2);
    std::cout << "After add level:" << cAdd->GetLevel() << " Encoding scale:" << log2(cAdd->GetScalingFactor()) << std::endl;
    std::cout << "q: " << cAdd->GetElements()[0].GetModulus().GetMSB() << std::endl;
    std::cout << "scale_deg: " << cAdd->GetNoiseScaleDeg() << std::endl;

    // Homomorphic subtraction
    auto cSub = cc->EvalSub(c1, c2);

    // Homomorphic scalar multiplication
    auto cScalar = cc->EvalMult(c1, 4.0);
    std::cout << "After scalar level:" << cScalar->GetLevel() << " Encoding scale:" << log2(cScalar->GetScalingFactor()) << std::endl;
    std::cout << "q: " << cAdd->GetElements()[0].GetModulus().GetMSB() << std::endl;
    std::cout << "scale_deg: " << cScalar->GetNoiseScaleDeg() << std::endl;

    // Homomorphic multiplication
    auto cMul = cc->EvalMult(c1, c2);
    std::cout << "After mult1 level:" << cMul->GetLevel() << " Encoding scale:" << log2(cMul->GetScalingFactor()) << std::endl;
    std::cout << "q: " << cMul->GetElements()[0].GetModulus().GetMSB() << std::endl;
    std::cout << "scale_deg: " << cMul->GetNoiseScaleDeg() << std::endl;

    auto cMul2 = cc->EvalMult(cMul, c2);
    std::cout << "After mult2 level:" << cMul2->GetLevel() << " Encoding scale:" << log2(cMul2->GetScalingFactor()) << std::endl;
    std::cout << "q: " << cMul2->GetElements()[0].GetModulus().GetMSB() << std::endl;
    std::cout << "scale_deg: " << cMul2->GetNoiseScaleDeg() << std::endl;

    auto cMul3 = cc->EvalMult(cMul2, c2);
    std::cout << "After mult3 level:" << cMul3->GetLevel() << " Encoding scale:" << log2(cMul3->GetScalingFactor()) << std::endl;
    std::cout << "q: " << cMul3->GetElements()[0].GetModulus().GetMSB() << std::endl;
    std::cout << "scale_deg: " << cMul3->GetNoiseScaleDeg() << std::endl;

    auto cMul4 = cc->EvalMult(cMul3, c2);
    std::cout << "After mult4 level:" << cMul4->GetLevel() << " Encoding scale:" << log2(cMul4->GetScalingFactor()) << std::endl;
    std::cout << "q: " << cMul4->GetElements()[0].GetModulus().GetMSB() << std::endl;
    std::cout << "scale_deg: " << cMul4->GetNoiseScaleDeg() << std::endl;

    std::cout << "Composite Degree: " << std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cMul4->GetCryptoParameters())->GetCompositeDegree() << std::endl;
    // algo->ModReduceInternalInPlace(cMul4, 1);
    cc->ModReduceInPlace(cMul4);
    std::cout << "After mult4 level:" << cMul4->GetLevel() << " Encoding scale:" << log2(cMul4->GetScalingFactor()) << std::endl;
    std::cout << "q: " << cMul4->GetElements()[0].GetModulus().GetMSB() << std::endl;
    std::cout << "scale_deg: " << cMul4->GetNoiseScaleDeg() << std::endl;

    algo->ModReduceInternalInPlace(cMul4, 1);
    std::cout << "After mult4 level:" << cMul4->GetLevel() << " Encoding scale:" << log2(cMul4->GetScalingFactor()) << std::endl;
    std::cout << "q: " << cMul4->GetElements()[0].GetModulus().GetMSB() << std::endl;
    std::cout << "scale_deg: " << cMul4->GetNoiseScaleDeg() << std::endl;

    Plaintext result;
    std::cout.precision(8);

    std::cout << std::endl << "Results of homomorphic computations: " << std::endl;

    // Decrypt the result of multiplication
    cc->Decrypt(keys.secretKey, cMul4, &result);
    result->SetLength(batchSize);
    std::cout << "cMul4 = " << result << std::endl;

    return 0;
}
