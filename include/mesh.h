#ifndef __MESH_H__
#define __MESH_H__


#include "vec3.h"
#include "vec4.h"
#include "aabb.h"
#include <vector>

struct mesh{
	std::vector<vec3d> vlist;
	std::vector<vec4i> flist;	
    std::vector<vec3d> nlist; //list of normals
    std::vector<vec3d> x1list; //list of principal direction X1
    std::vector<double> k1_list; //first principal curvature
    std::vector<double> k2_list; //second principal curvature
	aabb bb;
};

#endif
