#include "aabb.h"
#include <cmath>

#define EPS_R 2.22044604925031308e-16;

static bool almost_equal(double x, double y){
    double max_x_y_one = fmax(fmax(1.0,fabs(x)),fabs(y)) ;

    return fabs(x - y) <= max_x_y_one*EPS_R ;
}

std::ostream &operator<<(std::ostream &os, const aabb & bb){
    return os << bb.min_vec  << "  " << bb.max_vec ;
}


aabb::aabb():empty(true),min_vec(make_vec3d(0,0,0)),max_vec(make_vec3d(0,0,0)) {};


aabb::aabb( const vec3d &min_vec_, const vec3d &max_vec_):
    empty(false),min_vec(min_vec_),max_vec(max_vec_){}


void  aabb::grow(const vec3d &v){
    if(!empty){
        min_vec = cwisemin(min_vec,v);
        max_vec = cwisemax(max_vec,v);
    }else{
        min_vec = v;
        max_vec = v;
        empty = false;
    }     
}


void  aabb::grow(const aabb &bb){
    grow(bb.min_vec); grow(bb.max_vec);
}


void  aabb::scale(double factor){
    *this = this->scaled(factor);
}


aabb  aabb::scaled(double factor) const{
    vec3d center = get_center();
    vec3d half_size = get_half_size();
    return aabb(center-half_size*factor,center+half_size*factor);
}


vec3d  aabb::get_center() const{
    return (max_vec+min_vec)*0.5;
}


vec3d  aabb::aabb::get_half_size() const{
    return (max_vec-min_vec)*0.5;
}




bool  aabb::is_empty() const{
    return empty;
}



bool  aabb::is_planar() const{
    vec3d s = max_vec-min_vec;
    return  almost_equal(s.x,0.0) || almost_equal(s.y,0.0) || almost_equal(s.z,0.0) ;
}



void  aabb::set(const vec3d &min_vec_, const vec3d &max_vec_){

    min_vec = min_vec_; 
    max_vec = max_vec_;    
}



double aabb::get_surface_area() const{
    vec3d side = max_vec-min_vec;
    return 2.0*(side.x*side.y+side.x*side.z +side.y*side.z);
}



double aabb::get_vol() const{
    vec3d side = max_vec-min_vec;
    return side.x*side.y*side.z;
}
