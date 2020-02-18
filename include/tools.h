//
// Created by Administrator on 2020/2/18 0018.
//

#ifndef SHADOW_TOOLS_H
#define SHADOW_TOOLS_H
#include <time.h>
#include <iostream>
#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h>
#include <stdio.h>
#include <erand.h>
#include <vec.h>
#include <Ray.h>
#include <Sphere.h>

inline double clamp(double x)
{
    return x<0 ? 0 : x>1 ? 1 : x;
}
inline int toInt(double x)
{
    return int(pow(clamp(x),1/2.2)*255+.5);
}
Vec tc(0.0588, 0.361, 0.0941);
Vec sc = Vec(1,1,1)*.7;

Sphere spheres[] = {//Scene: radius, position, emission, color, material
        // center 50 40.8 62
        // floor 0
        // back  0
//  Sphere(1e5, Vec(50, 1e5+100, 0),  Vec(1,1,1)*1,Vec(),DIFF), //lite
//  Sphere(1e5, Vec(50, -1e5, 0),  Vec(),Vec(.3,.3,.1),DIFF), //grnd
//  Sphere(1e5, Vec(50, 1e5+100, 0),  Vec(0.761, 0.875, 1.00)*1.3,Vec(),DIFF),
//  //lite
        Sphere(1e5, Vec(50, 1e5+130, 0),  Vec(1,1,1)*1.3,Vec(),DIFF), //lite
        Sphere(1e2, Vec(50, -1e2+2, 47),  Vec(),Vec(1,1,1)*.7,DIFF), //grnd

        Sphere(1e4, Vec(50, -30, 300)+Vec(-sin(50*M_PI/180),0,cos(50*M_PI/180))*1e4, Vec(), Vec(1,1,1)*.99,SPEC),// mirr L
        Sphere(1e4, Vec(50, -30, 300)+Vec(sin(50*M_PI/180),0,cos(50*M_PI/180))*1e4,  Vec(), Vec(1,1,1)*.99,SPEC),// mirr R
        Sphere(1e4, Vec(50, -30, -50)+Vec(-sin(30*M_PI/180),0,-cos(30*M_PI/180))*1e4,Vec(), Vec(1,1,1)*.99,SPEC),// mirr FL
        Sphere(1e4, Vec(50, -30, -50)+Vec(sin(30*M_PI/180),0,-cos(30*M_PI/180))*1e4, Vec(), Vec(1,1,1)*.99,SPEC),// mirr


        Sphere(4, Vec(50,6*.6,47),   Vec(),Vec(.13,.066,.033), DIFF),//"tree"
        Sphere(16,Vec(50,6*2+16*.6,47),   Vec(), tc,  DIFF),//"tree"
        Sphere(11,Vec(50,6*2+16*.6*2+11*.6,47),   Vec(), tc,  DIFF),//"tree"
        Sphere(7, Vec(50,6*2+16*.6*2+11*.6*2+7*.6,47),   Vec(), tc,  DIFF),//"tree"

        Sphere(15.5,Vec(50,1.8+6*2+16*.6,47),   Vec(), sc,  DIFF),//"tree"
        Sphere(10.5,Vec(50,1.8+6*2+16*.6*2+11*.6,47),   Vec(), sc,  DIFF),//"tree"
        Sphere(6.5, Vec(50,1.8+6*2+16*.6*2+11*.6*2+7*.6,47),   Vec(), sc,  DIFF),//"tree"
};

inline bool intersect(const Ray &r, double &t, int &id)
{
    double n=sizeof(spheres)/sizeof(Sphere), d, inf=t=1e20;
    for(int i=int(n); i--;)
        if((d=spheres[i].intersect(r))&&d<t)
        {
            t=d;
            id=i;
        }
    return t<inf;
}
Vec radiance(const Ray &r, int depth, unsigned short *Xi)
{
    double t;                               // distance to intersection
    int id=0;                               // id of intersected object
    if (!intersect(r, t, id))
        return Vec(); // if miss, return black
    const Sphere &obj = spheres[id];        // the hit object
    Vec x=r.o+r.d*t, n=(x-obj.p).norm(), nl=n.dot(r.d)<0?n:n*-1, f=obj.c;
    double p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl
    if (++depth>5)
        if (erand48(Xi)<p)
            f=f*(1/p);
        else
            return obj.e; //R.R.
    if (obj.refl == DIFF)                   // Ideal DIFFUSE reflection
    {
        double r1=2*M_PI*erand48(Xi), r2=erand48(Xi), r2s=sqrt(r2);
        Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1):Vec(1))%w).norm(), v=w%u;
        Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();
        return obj.e + f.mult(radiance(Ray(x,d),depth,Xi));
    }
    else if (obj.refl == SPEC)              // Ideal SPECULAR reflection
        return obj.e + f.mult(radiance(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi));
    Ray reflRay(x, r.d-n*2*n.dot(r.d));     // Ideal dielectric REFRACTION
    bool into = n.dot(nl)>0;                // Ray from outside going in?
    double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=r.d.dot(nl), cos2t;
    if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    // Total internal reflection
        return obj.e + f.mult(radiance(reflRay,depth,Xi));
    Vec tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
    double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));
    double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);
    return obj.e + f.mult(depth>2 ? (erand48(Xi)<P ?   // Russian roulette
                                     radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) :
                          radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr);
}
#endif //SHADOW_TOOLS_H
