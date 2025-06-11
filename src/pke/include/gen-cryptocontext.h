//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2014-2022, NJIT, Duality Technologies Inc. and other contributors
//
// All rights reserved.
//
// Author TPOC: contact@openfhe.org
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//==================================================================================

/*
  构造基于提供参数集的 CryptoContext
 */

/*
* 如何通过调用 GenCryptoContext() 生成 CryptoContext
*
* 1. 选择您想要使用的方案。在我们的教程示例中，我选择 CKKS。
* 2. 您的代码必须包含此头文件以及方案特定的上下文生成器头文件
*    (scheme/<scheme>/cryptocontext-<scheme>.h)：
*       #include "scheme/ckks/cryptocontext-ckks.h"
*       #include "gen-cryptocontext.h"
* 3. 创建一个参数对象作为 GenCryptoContext() 的参数。其通用形式如下：
*       CCParams<GeneratorName<Element>> parameters
*    其中：
*    - GeneratorName 是在 cryptocontext-<scheme>.h 中定义的类的名称。在我们的例子中，
*      它是 CryptoContextCKKS。
*    - Element 是一个模板参数，表示整数格子。因此，它可以保持为 Element 或替换为
*      Poly、NativePoly 或 DCRTPoly。在这里我们保持为 "Element"。
*      因此我们可以添加以下代码：
*       CCParams<CryptoContextCKKS<Element>> parameters;
* 4. 使用 CCParams<CryptoContextCKKS<Element>> 的 set 函数调整参数值，
*    因为对象是使用 scheme/cryptocontextparams-defaults.h 中的默认值创建的。
* 5. 调用 GenCryptoContext() 生成 CryptoContext。
*
* 现在您的代码应如下所示：
*       #include "scheme/ckks/cryptocontext-ckks.h"
*       #include "gen-cryptocontext.h"
*       ...........................................
*       CCParams<CryptoContextCKKS<Element>> parameters;
*       parameters.SetMultiplicativeDepth(1);
*       parameters.SetScalingModSize(50);
*       parameters.SetBatchSize(8);
*       parameters.SetSecurityLevel(HEStd_NotSet);
*       parameters.SetRingDim(16);
*
*       auto cryptoContext = GenCryptoContext(parameters);
*
*       cryptoContext->Enable(ENCRYPTION);
*       cryptoContext->Enable(KEYSWITCH);
*       cryptoContext->Enable(LEVELEDSHE);
*       ...........................................
*
* More examples can be found in src/pke/unittest/UnitTestAutomorphism.cpp or in
* src/pke/unittest/UnitTestEvalMult.cpp.
*/

#ifndef _GEN_CRYPTOCONTEXT_H_
#define _GEN_CRYPTOCONTEXT_H_

namespace lbcrypto {

// forward declarations (don't include headers as compilation fails when you do)
template <typename T>
class CCParams;

template <typename T>
typename T::ContextType GenCryptoContext(const CCParams<T>& params) {
    return T::genCryptoContext(params);
}

}  // namespace lbcrypto

#endif  // _GEN_CRYPTOCONTEXT_H_
