// #include "openfhe.h"
// #include "slotstocoeffs-utils.h"

// int main(){
//     // testSlots2Coeffs();
//     // mytestSlots2Coeffs();

//     testStC_CtS();
//     mytestStC_CtS();

//     testCtS_StC();
//     mytestCtS_StC();

//     testStC_CtSReal();
    
//     return 0;
// }

#include "openfhe.h"
#include "homoencrypt-compute.h"

int main(){
    // 创建 HomoEncryptCompute 实例
    HomoEncryptCompute hec;
    
    // 运行测试函数
    // 稀疏打包测试
    // hec.testSlots2Coeffs();
    hec.testStC_CtS();
    // hec.testCtS_StC();
    // hec.testStC_CtSReal();
    
    // 全打包测试
    // hec.mytestSlots2Coeffs();
    // hec.mytestStC_CtS();
    // hec.mytestCtS_StC();
    
    return 0;
}