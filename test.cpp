#include "cppgf.h"

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <random>
#include <chrono>

int main(void)
{
    cppgf::GaloisFieldTable<0x11D> table;
    printf("exp\n");
    for(uint32_t i = 0; i < 256; ++i) {
        if(0 == (i % 16)) {
            printf("\n");
        }
        printf("0x%XU, ", table.exp_[i]);
    }
    printf("\n");
    printf("log\n");
    for(uint32_t i = 0; i < 256; ++i) {
        if(0 == (i % 16)) {
            printf("\n");
        }
        printf("0x%XU, ", table.log_[i]);
    }
    return 0;
}
