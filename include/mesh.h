#ifndef __MESH_H__
#define __MESH_H__


#include "vec3.h"
#include "vec4.h"
#include "aabb.h"
#include <vector>

struct mesh{
	std::vector<vec3d> vlist;
	std::vector<vec4i> flist;	
	aabb bb;
};

#endif
