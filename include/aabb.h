#ifndef __AABB_H__
#define __AABB_H__

#include "vec3.h"
#include "vec4.h"


struct aabb{
    bool empty;
    vec3d min_vec;
    vec3d max_vec;
    
    //methods
    aabb();
    aabb( const vec3d &min_vec, const vec3d &max_vec);
    void    grow(const vec3d &v);
    void    grow(const aabb &bb);
    void    scale(double factor);
    aabb    scaled(double factor) const;
    vec3d   get_center() const;
    vec3d   get_half_size() const ;
    bool    is_empty() const;
    bool    is_planar() const;
    void    set(const vec3d &min_vec, const vec3d &max_vec);
    double  get_surface_area() const;
    double  get_vol() const;   
};

std::ostream &operator<<(std::ostream &os, const aabb & bb);

#endif
