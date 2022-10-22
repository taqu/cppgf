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
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdlib>

namespace cppgf
{
template<uint32_t P>
class GaloisFieldTable
{
public:
    GaloisFieldTable();
    static uint8_t mulNoLUT(uint32_t x0, uint32_t x1);

    uint8_t exp_[256];
    uint8_t log_[256];
};

template<uint32_t P>
GaloisFieldTable<P>::GaloisFieldTable()
{
    ::memset(exp_, 0, sizeof(uint8_t) * 256);
    ::memset(log_, 0, sizeof(uint8_t) * 256);
    uint32_t x = 1;
    for(uint32_t i = 0; i < 256; ++i) {
        exp_[i] = x;
        log_[x] = i;
        x = mulNoLUT(x, 2);
    }
}

template<uint32_t P>
uint8_t GaloisFieldTable<P>::mulNoLUT(uint32_t x0, uint32_t x1)
{
    uint32_t r = 0;
    while(0 < x1) {
        if(0x01U == (x1 & 0x01U)) {
            r = r ^ x0;
        }
        x1 = x1 >> 1;
        x0 = x0 << 1;
        if(256U <= x0) {
            x0 = x0 ^ P;
        }
    }
    return static_cast<uint8_t>(r&0xFFU);
}

/**
 * https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders
 */
template<uint32_t P>
class GaloisField
{
public:
    using this_type = GaloisField<P>;

    explicit GaloisField(uint32_t size);
    ~GaloisField();

    uint32_t capacity() const;
    uint32_t size() const;
    void resize(uint32_t size);

    uint32_t eval(uint8_t x) const;

    const uint32_t& operator[](uint32_t x) const
    {
        return polynomial_[x];
    }

    uint32_t& operator[](uint32_t x)
    {
        return polynomial_[x];
    }

    static uint8_t add(uint8_t x0, uint8_t x1);
    static uint8_t sub(uint8_t x0, uint8_t x1);
    static uint8_t mul(uint8_t x0, uint8_t x1);
    static uint8_t div(uint8_t x0, uint8_t x1);
    static uint8_t pow(uint8_t x, uint8_t p);
    static uint8_t inverse(uint8_t x);

private:
    GaloisField(const GaloisField&) = delete;
    GaloisField& operator=(const GaloisField&) = delete;
    static GaloisFieldTable<P> table_;

    uint32_t capacity_;
    uint32_t size_;
    uint32_t* polynomial_;
};

template<uint32_t P>
GaloisField<P>::GaloisField(uint32_t size)
    : size_(size)
    , polynomial_(nullptr)
{
    assert(0 < size_);
    capacity_ = (size_ + 0x15U) & ~0x15U;
    polynomial_ = ::malloc(sizeof(uint32_t) * capacity_);
    ::memset(polynomial_, 0, sizeof(uint32_t) * capacity_);
}

template<uint32_t P>
GaloisField<P>::~GaloisField()
{
    ::free(polynomial_);
    polynomial_ = nullptr;
}

template<uint32_t P>
uint32_t GaloisField<P>::capacity() const
{
    return capacity_;
}

template<uint32_t P>
uint32_t GaloisField<P>::size() const
{
    return size_;
}

template<uint32_t P>
void GaloisField<P>::resize(uint32_t size)
{
    if(size <= capacity_) {
        return;
    }
    ::free(polynomial_);
    capacity_ = (size + 0x15U) & ~0x15U;
    size_ = size;
    polynomial_ = ::malloc(sizeof(uint32_t) * capacity_);
    ::memset(polynomial_, 0, sizeof(uint32_t) * capacity_);
}

template<uint32_t P>
uint32_t GaloisField<P>::eval(uint8_t x) const
{
    assert(0 < size_);
    uint8_t y = polynomial_[0];
    for(uint32_t i = 1; i < size_; ++i) {
        y = mul(y, x) ^ polynomial_[i];
    }
    return y;
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
    if(0 == x0 || 0 == x1) {
        return 0;
    }
    uint32_t sum = static_cast<uint32_t>(table_.log_[x0]) + table_.log_[x1];
    if(255 <= sum) {
        sum -= 255;
    }
    return table_.exp_[sum];
}

template<uint32_t P>
uint8_t GaloisField<P>::div(uint8_t x0, uint8_t x1)
{
    if(0 == x0) {
        return 0;
    }
    if(0 == x1) {
        return static_cast<uint8_t>(-1);
    }
    int32_t diff = static_cast<int32_t>(table_.log_[x0]) - table_.log_[x1];
    if(diff < 0) {
        diff += 255;
    }
    return table_.exp_[diff];
}

template<uint32_t P>
uint8_t GaloisField<P>::pow(uint8_t x, uint8_t p)
{
    return table_.exp_[(table_.log_[x] * p) % 255U];
}

template<uint32_t P>
uint8_t GaloisField<P>::inverse(uint8_t x)
{
    return table_.exp_[255U - table_.log_[x]];
}

template<uint32_t P>
void add(GaloisField<P>& result, const GaloisField<P>& x0, const GaloisField<P>& x1)
{
    uint32_t size = (std::max)(x0.size(), x1.size());
    result.resize(size);

    for(uint32_t i = 0; i < x0.size(); ++i) {
        result[i + size - x0.size()] = x0[i];
    }
    for(uint32_t i = 0; i < x1.size(); ++i) {
        result[i + size - x1.size()] ^= x1[i];
    }
}

template<uint32_t P>
void mul(GaloisField<P>& result, const GaloisField<P>& x0, const GaloisField<P>& x1)
{
    uint32_t total = x0.size() + x1.size() - 1;
    result.resize(total);
    for(uint32_t i = 0; i < total; ++i) {
        result[i] = 0;
    }
    for(uint32_t i = 0; i < x1.size(); ++i) {
        for(uint32_t j = 0; j < x0.size(); ++j) {
            result[i + j] ^= GaloisField<P>::mul(x0[j], x1[i]);
        }
    }
}
} // namespace cppgf
#endif //INC_CPPGF_H_