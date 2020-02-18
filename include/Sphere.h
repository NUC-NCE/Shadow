//
// Created by Administrator on 2020/2/18 0018.
//

#ifndef SHADOW_SPHERE_H
#define SHADOW_SPHERE_H

#include <vec.h>
#include <Ray.h>
enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()
struct Sphere
{
    double rad;       // radius
    Vec p, e, c;      // position, emission, color
    Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)
    Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_):
            rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}
    double intersect(const Ray &r) const   // returns distance, 0 if nohit
    {
        Vec op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        double t, eps=1e-4, b=op.dot(r.d), det=b*b-op.dot(op)+rad*rad;
        if (det<0)
            return 0;
        else
            det=sqrt(det);
        return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0);
    }
};
#endif //SHADOW_SPHERE_H
