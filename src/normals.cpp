#include <normals.h>

void compute_normals(
    const std::vector<vec3f> &vlist, //list of vertex
    const std::vector<vec3i> &flist, //list of faces
    std::vector<vec3f> &nlist //list of normals              
){
    nlist.clear();
    nlist.reserve(vlist.size());
    std::vector< 
        std::vector<int> 
    >vert_to_face_map(vlist.size());
    
    std::vector<vec3f> face_normal_list;
    face_normal_list.reserve(flist.size()); 
    
    
    //preallocate memory space for each vertex 
    for (size_t  v = 0; v< vlist.size(); v++){
        //we expect on average that each vertex touches six faces 
        vert_to_face_map[v].reserve(6);
    }
    
    //for each vertex make a list of faces which share that particular vertex 
    //and compute the normal of each face
    for(size_t face = 0; face < flist.size(); face++){
       vert_to_face_map[flist[face].a].push_back(face);
       vert_to_face_map[flist[face].b].push_back(face);
       vert_to_face_map[flist[face].c].push_back(face);
       
       //compute normal for this face
       const vec3f &A = vlist[flist[face].a];
       const vec3f &B = vlist[flist[face].b];
       const vec3f &C = vlist[flist[face].c];
       
       face_normal_list.push_back(cross((B-A),C-A));
    }
    
    //from here 
    //we know for each vertex a list of faces that share that particular vertex
    //and the normal of each face

    //Compute  normal for each vertex by summing the normals of faces that use 
    // that vertex and normalizing.
    for (size_t v = 0; v< vlist.size(); v++){
    
        vec3f N = make_vec3f(0,0,0);
        for(size_t face = 0; face < vert_to_face_map[v].size(); face++){
            N = N + face_normal_list[vert_to_face_map[v][face]];
        }        

        nlist.push_back(normalize(N));
    }
}