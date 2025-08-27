#include "homoencrypt-compute.h"

int main(int argc, char* argv[]){
    // 检查参数数量
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <0|1>" << std::endl;
        std::cerr << "  0: Run StC first bootstrapping" << std::endl;
        std::cerr << "  1: Run CtS first bootstrapping" << std::endl;
        return 1;
    }

    // 解析参数
    int mode = std::atoi(argv[1]);
    if (mode != 0 && mode != 1 && mode != 2 && mode != 3) {
        std::cerr << "Error: Parameter must be 0 or 1 or 2 or 3" << std::endl;
        return 1;
    }

    HomoEncryptCompute hec;
    int depth = 30;
    std::vector<uint32_t> levelBudget = {3, 3};

    std::vector<double> x = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
    size_t encodedLength  = x.size();
    if(mode == 0) {
        std::cout << "Running StC first bootstrapping--full packing..." << std::endl;
        hec.generate_context_stcfirstbts_toy();
        std::cout << "StC first BTS Context Generate Scussfully!" << std::endl;

        // one extra level
        Plain ptxt = hec.encode(x, depth-levelBudget[1]-1);
        ptxt->SetLength(encodedLength);
        std::cout << "Input: " << ptxt << std::endl;

        Cipher ciph = hec.encrypt(ptxt);
        std::cout << "Initial number of levels remaining: " << depth - ciph->GetLevel() << std::endl;
        
        // Perform the bootstrapping operation. The goal is to increase the number of levels remaining
        // for HE computation.
        std::cout << "Start StC First Bootstrapping..." << std::endl;
        auto ciphertextAfter = hec.bootstrapStCfirst(ciph);
        std::cout << "StC First Bootstrapping done!" << std::endl;

        std::cout << "Number of levels remaining after bootstrapping: "
                << depth - ciphertextAfter->GetLevel() - (ciphertextAfter->GetNoiseScaleDeg() - 1) << std::endl
                << std::endl;

        Plain result;
        result = hec.decrypt(ciphertextAfter);
        result->SetLength(encodedLength);
        std::cout << "Output after bootstrapping \n\t" << result << std::endl;
    } else if(mode == 1) {
        std::cout << "Running CtS first bootstrapping--full packing..." << std::endl;
        hec.generate_context_ctsfirstbts_toy();
        std::cout << "CtS first BTS Context Generate Scussfully!" << std::endl;

        Plain ptxt = hec.encode(x, depth-1);
        ptxt->SetLength(encodedLength);
        std::cout << "Input: " << ptxt << std::endl;

        Cipher ciph = hec.encrypt(ptxt);
        std::cout << "Initial number of levels remaining: " << depth - ciph->GetLevel() << std::endl;

        // Perform the bootstrapping operation. The goal is to increase the number of levels remaining
        // for HE computation.
        std::cout << "Start CtS First Bootstrapping..." << std::endl;
        auto ciphertextAfter = hec.bootstrapCtSfirst(ciph);
        std::cout << "CtS First Bootstrapping done!" << std::endl;

        std::cout << "Number of levels remaining after bootstrapping: "
                << depth - ciphertextAfter->GetLevel() - (ciphertextAfter->GetNoiseScaleDeg() - 1) << std::endl
                << std::endl;

        Plain result;
        result = hec.decrypt(ciphertextAfter);
        result->SetLength(encodedLength);
        std::cout << "Output after bootstrapping \n\t" << result << std::endl;
    } else if(mode == 2){
        std::cout << "Running StC first bootstrapping--sparse packing..." << std::endl;
        int num_slots = 1 << 3;
        hec.generate_context_stcfirstbts_toy(num_slots);
        std::cout << "StC first BTS Context Generate Scussfully!" << std::endl;

        // one extra level
        Plain ptxt = hec.encode(x, depth-levelBudget[1]-1);
        ptxt->SetLength(encodedLength);
        std::cout << "Input: " << ptxt << std::endl;

        Cipher ciph = hec.encrypt(ptxt);
        std::cout << "Initial number of levels remaining: " << depth - ciph->GetLevel() << std::endl;
        
        // Perform the bootstrapping operation. The goal is to increase the number of levels remaining
        // for HE computation.
        std::cout << "Start StC First Bootstrapping..." << std::endl;
        auto ciphertextAfter = hec.bootstrapStCfirst(ciph);
        std::cout << "StC First Bootstrapping done!" << std::endl;

        std::cout << "Number of levels remaining after bootstrapping: "
                << depth - ciphertextAfter->GetLevel() - (ciphertextAfter->GetNoiseScaleDeg() - 1) << std::endl
                << std::endl;

        Plain result;
        result = hec.decrypt(ciphertextAfter);
        result->SetLength(encodedLength);
        std::cout << "Output after bootstrapping \n\t" << result << std::endl;

    } else{
        std::cout << "Running CtS first bootstrapping--sparsely packing..." << std::endl;
        int num_slots = 1 << 3;
        hec.generate_context_ctsfirstbts_toy(num_slots);
        std::cout << "CtS first BTS Context Generate Scussfully!" << std::endl;

        Plain ptxt = hec.encode(x, depth-1);
        ptxt->SetLength(encodedLength);
        std::cout << "Input: " << ptxt << std::endl;

        Cipher ciph = hec.encrypt(ptxt);
        std::cout << "Initial number of levels remaining: " << depth - ciph->GetLevel() << std::endl;

        // Perform the bootstrapping operation. The goal is to increase the number of levels remaining
        // for HE computation.
        std::cout << "Start CtS First Bootstrapping..." << std::endl;
        auto ciphertextAfter = hec.bootstrapCtSfirst(ciph);
        std::cout << "CtS First Bootstrapping done!" << std::endl;

        std::cout << "Number of levels remaining after bootstrapping: "
                << depth - ciphertextAfter->GetLevel() - (ciphertextAfter->GetNoiseScaleDeg() - 1) << std::endl
                << std::endl;

        Plain result;
        result = hec.decrypt(ciphertextAfter);
        result->SetLength(encodedLength);
        std::cout << "Output after bootstrapping \n\t" << result << std::endl;

    }
    return 0;
}