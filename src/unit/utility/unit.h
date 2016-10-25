#pragma once

#include<cstdlib>
#include<cstdio>

#define STR_VALUE(arg) #arg
#define CHECK(x) \
    if(!(x)) { \
        printf("Detected a false entry: %s\n",STR_VALUE(x)); \
        exit(EXIT_FAILURE);\
    }
