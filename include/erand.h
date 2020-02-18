//
// Created by Administrator on 2020/2/18 0018.
//

#ifndef SHADOW_ERAND_H
#define SHADOW_ERAND_H
#include <time.h>
#include <iostream>
#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h>
#include <stdio.h>
using namespace std;
double erand48(unsigned short xsubi[3])
{
    return (double)rand() / (double)RAND_MAX;
}
#endif //SHADOW_ERAND_H
