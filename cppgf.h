#ifndef INC_CPPGF_H_
#define INC_CPPGF_H_
/**
@file cppecc.h
@author t-sakai

# License
This software is distributed under two licenses, choose whichever you like.

## MIT License
Copyright (c) 2022 Takuro Sakai

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Public Domain
This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org>
*/
#include <cstdint>
#include <cstdlib>

namespace cppgf
{
/**
 * https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders
 */
template<uint32_t P>
class GaloisField
{
public:
    using this_type = GaloisField<P>;

    explicit GaloisField(uint32_t);
    ~GaloisField();

    void add(this_type& x0, this_type& x1);
    void mul(this_type& x0, this_type& x1);
    void div(this_type& x0, this_type& x1);

    uint32_t eval(uint8_t x) const;

private:
    GaloisField(const GaloisField&) = delete;
    GaloisField& operator=(const GaloisField&) = delete;

    static uint8_t add(uint8_t x0, uint8_t x1);
    static uint8_t sub(uint8_t x0, uint8_t x1);
    static uint8_t mul(uint8_t x0, uint8_t x1);
    static uint8_t div(uint8_t x0, uint8_t x1);

    void mulNoLUT(uint32_t x0, uint32_t x1, uint32_t prim);

    uint8_t exp_[256];
    uint8_t log_[256];
};

template<uint32_t P>
GaloisField<P>::GaloisField(uint32_t dim)
{
    ::memset(exp_, 0, sizeof(uint8_t) * 256);
    ::memset(log_, 0, sizeof(uint8_t) * 256);
    uint8_t x = 1;
    for(uint8_t i = 0; i < 256; ++i) {
        exp_[i] = x;
        log_[i] = i;
        x = mulNoLUT(x, 2, P);
    }
}

template<uint32_t P>
uint8_t GaloisField<P>::add(uint8_t x0, uint8_t x1)
{
    return x0 ^ x1;
}

template<uint32_t P>
uint8_t GaloisField<P>::sub(uint8_t x0, uint8_t x1)
{
    return x0 ^ x1;
}

template<uint32_t P>
uint8_t GaloisField<P>::mul(uint8_t x0, uint8_t x1)
{
    return x0 ^ x1;
}

template<uint32_t P>
uint8_t GaloisField<P>::div(uint8_t x0, uint8_t x1)
{
    if(0 == a || 0 == b) {
        return 0;
    }
    cppecc_u32 sum = CPPECC_STATIC_CAST(cppecc_s32)(gflog[a]) + gflog[b];
    if(CPPECC_GF_NW1 <= sum) {
        sum -= CPPECC_GF_NW1;
    }
    return gfexp[sum];
}

template<uint32_t P>
void GaloisField<P>::mulNoLUT(uint32_t x0, uint32_t x1, uint32_t prim)
{
    uint32_t r = 0;
    while(0 < x1) {
        if(0 != (x1 & 0x01U)) {
            r = r ^ x0;
        }
        x1 = x1 >> 1;
        x0 = x0 << 1;
        if(256U <= x0) {
            x0 = x0 ^ prim;
        }
    }
}
} // namespace cppgf
#endif //INC_CPPGF_H_