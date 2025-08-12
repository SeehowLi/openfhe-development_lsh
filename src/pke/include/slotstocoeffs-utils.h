#ifndef SLOTS_TO_COEFFS_UTILS_H
#define SLOTS_TO_COEFFS_UTILS_H

#include <iostream>
#include <vector>
#include <complex>
#include <memory>
#include "openfhe.h"

using namespace lbcrypto;

// Function to convert slots to coefficients
Ciphertext<DCRTPoly> SlotsToCoeffs(const CryptoContext<DCRTPoly>& cryptoContext, const Ciphertext<DCRTPoly>& ciph);

// Function to convert coefficients to slots
Ciphertext<DCRTPoly> CoeffsToSlots(const CryptoContext<DCRTPoly>& cryptoContext, const Ciphertext<DCRTPoly>& ciph);

// Function to set up the SlotsToCoeffs operation
void EvalSlotsToCoeffsSetup(
    const CryptoContext<DCRTPoly>& cryptoContext,
    uint32_t levelBudget,
    uint32_t numSlots,
    uint32_t lDec);

// Function to set up the StC and CtS operation
void EvalStC_CtSSetup(
    const CryptoContext<DCRTPoly>& cryptoContext,
    uint32_t levelBudget,
    uint32_t numSlots,
    uint32_t lDec,
    uint32_t lEnc);

// Test function for SlotsToCoeffs
void testSlots2Coeffs();

void mytestSlots2Coeffs();

// Test function for SlotsToCoeffs and CoeffsToSlots--For Complex
// Sparse packing
void testStC_CtS();
void testCtS_StC();
// Full packing
void mytestStC_CtS();
void mytestCtS_StC();

// Test function for SlotsToCoeffs and CoeffsToSlots--For Real
// Sparse packing
void testStC_CtSReal();
void testCtS_StReal();
// Full packing
void mytestStC_CtSReal();
void mytestCtS_StCReal();

#endif // SLOTS_TO_COEFFS_UTILS_H
