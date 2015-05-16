#ifndef __NORMALS_AND_CURVATURES_H__
#define __NORMALS_AND_CURVATURES_H__

#include "vec3.h"
#include <vector>
#include <mesh.h>

void compute_normals(
    const std::vector<vec3f> &vlist, //list of vertex
    const std::vector<vec3i> &flist, //list of faces
    std::vector<vec3f> &nlist //list of normals              
);

void recompute_flat_verts_normals_and_curvature(mesh &global_mesh);

void write_normcurv_ssv_mt(std::ostream &os,const  mesh &m);

void write_normcurv_bin_mt(std::ostream &os,const  mesh &m);


#endif
