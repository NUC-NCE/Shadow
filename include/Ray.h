//
// Created by Administrator on 2020/2/18 0018.
//

#ifndef SHADOW_RAY_H
#define SHADOW_RAY_H

#include <vec.h>
struct Ray
{
    Vec o, d;
    Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};
#endif //SHADOW_RAY_H
