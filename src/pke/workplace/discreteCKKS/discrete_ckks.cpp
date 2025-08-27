/*
 * @Author: SeehowLi lsh0126@nudt.edu.cn
 * @Date: 2025-07-03 21:55:24
 * @LastEditors: SeehowLi lsh0126@nudt.edu.cn
 * @LastEditTime: 2025-07-10 20:06:50
 * @FilePath: \openfhe-development\src\pke\examples\discrete_ckks.cpp
 * @Description: 用于实现离散CKKS的示例代码
 * 
 * Copyright (c) 2025 by $SeehowLi lsh0126@nudt.edu.cn, All Rights Reserved. 
 */

#include "homoencrypt-compute.h"
#include "openfhe.h"

using namespace lbcrypto;
HomoEncryptCompute hec;

void discrete_ckks_example();

int main() {
    discrete_ckks_example();

    return 0;
}

void discrete_ckks_example() {
    std::cout << "--------------------------------- DISCRETE CKKS EXAMPLE ---------------------------------"
              << std::endl;
    int ringdim = 1 << 12;
    int num_slots = ringdim / 2;
    int level = 8;
    int dcrtBits = 50;
    int q0bits = 50;
    std::vector<uint32_t> levelBudget = {3, 3};
    int lStC = 1;
    int lCtS = level - levelBudget[1] - levelBudget[0] - 1;
    // int t = 32;
    // int l =1;
    int slotsinuse = 8;
    hec.generate_context_bts(levelBudget, num_slots, level, ringdim, lStC, lCtS, dcrtBits, q0bits);
    
    std::vector<double> inputvalue(slotsinuse, 0.0);
    inputvalue = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
    std::cout << "Input value: ";
    auto pt = hec.encode(inputvalue, level - levelBudget[1], num_slots);
    std::cout << pt << std::endl;

    auto ct = hec.encrypt(pt);
    std::cout << "Current level: " << ct->GetLevel() << std::endl;
    std::cout << "Noise deg: " << ct->GetNoiseScaleDeg() << std::endl;
    // ct = hec.mult(ct, 1.0);
    // std::cout << "After square, current level: " << ct->GetLevel() << std::endl;
    // std::cout << "Noise deg: " << ct->GetNoiseScaleDeg() << std::endl;
    auto sk = hec.getsecretkey();
    ct = hec.SlotsToCoeffs(ct);
    hec.modreduceinternalinplace(ct, 1);
    std::cout << "After SlotsToCoeffs, current level: " << ct->GetLevel() << std::endl;
    std::cout << "Noise deg: " << ct->GetNoiseScaleDeg() << std::endl;
    Poly xx1=hec.myDecrypt(ct, sk);
    BigInteger q;
    BigInteger x;
    std::cout << "slots:" << ct->GetSlots() << std::endl;
    std::cout << "modulus size: " << ct->GetElements()[0].GetModulus().GetMSB() << std::endl;
    // for(int i=0; i<1000; i++) {
    //     std::cout << i << " ";
    //     std::cout << " ";
    //     q = xx1.GetModulus();
    //      if(xx1[i]>q/2)
    //         x = q - xx1[i];
    //     else
    //         x= xx1[i];
    //     std::cout << x.ConvertToDouble()/ct->GetScalingFactor() << std::endl;

    // }
    
    

    Plain result;
    result = hec.decrypt(ct);
    result->SetLength(slotsinuse);
    std::cout << "Decrypted result: " << result << std::endl;

}