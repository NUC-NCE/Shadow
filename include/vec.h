//
// Created by Administrator on 2020/2/18 0018.
//

#ifndef SHADOW_VEC_H
#define SHADOW_VEC_H
#include <time.h>
#include <iostream>
#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h>
#include <stdio.h>
struct Vec{
    double x,y,z;
    Vec(double x_=0, double y_=0, double z_=0)
    {
        x=x_;
        y=y_;
        z=z_;
    }
    Vec operator+(const Vec &b) const
    {
        return Vec(x+b.x,y+b.y,z+b.z);
    }
    Vec operator-(const Vec &b) const
    {
        return Vec(x-b.x,y-b.y,z-b.z);
    }
    Vec operator*(double b) const
    {
        return Vec(x*b,y*b,z*b);
    }
    Vec mult(const Vec &b) const
    {
        return Vec(x*b.x,y*b.y,z*b.z);
    }
    Vec& norm()
    {
        return *this = *this * (1/sqrt(x*x+y*y+z*z));
    }
    double dot(const Vec &b) const
    {
        return x*b.x+y*b.y+z*b.z;    // cross:
    }
    Vec operator%(Vec&b)
    {
        return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);
    }

};

#endif //SHADOW_VEC_H
