#ifndef __NORMALS_H__
#define __NORMALS_H__

#include "vec3.h"
#include <vector>

void compute_normals(
    const std::vector<vec3f> &vlist, //list of vertex
    const std::vector<vec3i> &flist, //list of faces
    std::vector<vec3f> &nlist //list of normals              
);

  



#endif
