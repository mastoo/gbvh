#ifndef __TRIS_H__
#define __TRIS_H__

#include <iostream>
#include "mesh.h"

//write the mesh in Moller–Trumbore format
void write_tris_ssv_mt(std::ostream &os,const  mesh &m, bool normal_required);


void write_tris_bin_mt(std::ostream &os,const  mesh &m, bool normal_required);


#endif
