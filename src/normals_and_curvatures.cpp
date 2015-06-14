#include <normals_and_curvatures.h>
#include <mesh.h>
#include <packed_struct.h>
#include <thrust/host_vector.h>
#include <iomanip>

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


void recompute_flat_verts_normals_and_curvature(mesh &global_mesh){
    
    
    {   //recompute faces and verts
        size_t new_num_verts = 3*global_mesh.flist.size();


        std::vector<vec3d> new_verts_list;
        std::vector<vec4i> new_faces_list;
        
        new_verts_list.reserve(new_num_verts);
        new_faces_list.reserve(global_mesh.flist.size());
        
        for(size_t f = 0; f < global_mesh.flist.size(); f++){
            vec4i  old_faces_face = global_mesh.flist[f];
            
            int ia = new_verts_list.size();        
            new_verts_list.push_back(global_mesh.vlist[old_faces_face.a]);
            int ib = new_verts_list.size();        
            new_verts_list.push_back(global_mesh.vlist[old_faces_face.b]);
            int ic = new_verts_list.size();        
            new_verts_list.push_back(global_mesh.vlist[old_faces_face.c]);
            
            vec4i new_face = make_vec4i(ia,ib,ic,old_faces_face.d);       
            
            new_faces_list.push_back(new_face);    
            
        }
        
        global_mesh.vlist.swap(new_verts_list);
        global_mesh.flist.swap(new_faces_list);
    }
    
    
    //compute nornals and principal curvature vectors
    
    global_mesh.nlist.resize(global_mesh.vlist.size());
    global_mesh.x1list.resize(global_mesh.vlist.size());
    global_mesh.k1_list.resize(global_mesh.vlist.size());
    global_mesh.k2_list.resize(global_mesh.vlist.size());
    
    for(size_t f = 0; f < global_mesh.flist.size(); f++){
        
        const vec3d A = global_mesh.vlist[global_mesh.flist[f].a];
        const vec3d B = global_mesh.vlist[global_mesh.flist[f].b];
        const vec3d C = global_mesh.vlist[global_mesh.flist[f].c];
        //compute faces normal
        const vec3d n_v =  normalize(cross(B-A,C-A));
        //compute the principal curvature vector
        const vec3d X1_v = normalize(B-A);
        
        //here we assume flat triangles
        const double k1 = 0;
        const double k2 = 0;

        global_mesh.nlist[global_mesh.flist[f].a] = n_v;
        global_mesh.nlist[global_mesh.flist[f].b] = n_v;
        global_mesh.nlist[global_mesh.flist[f].c] = n_v;
        
        global_mesh.x1list[global_mesh.flist[f].a] = X1_v;
        global_mesh.x1list[global_mesh.flist[f].b] = X1_v;
        global_mesh.x1list[global_mesh.flist[f].c] = X1_v;
        

        global_mesh.k1_list[global_mesh.flist[f].a] = k1;
        global_mesh.k1_list[global_mesh.flist[f].b] = k1;
        global_mesh.k1_list[global_mesh.flist[f].c] = k1;

        
        global_mesh.k2_list[global_mesh.flist[f].a] = k2;
        global_mesh.k2_list[global_mesh.flist[f].b] = k2;
        global_mesh.k2_list[global_mesh.flist[f].c] = k2;        
        
    }
    
    
}




void write_normcurv_ssv_mt(std::ostream &os,const  mesh &m){
    os << "#number of vertices"<< std::endl;
    os << m.vlist.size() <<std::endl;
    os << "##start normal list" << std::endl;
    os << "# x y z" << std::endl;
    for(size_t i = 0; i < m.nlist.size(); i++){
        os << std::setprecision(15) << m.nlist[i].x << " " << m.nlist[i].y << " "  << m.nlist[i].z << " " << std::endl;
    }   
    os << "##end normals list" << std::endl;
    os << "##start first principal direction  list" << std::endl;
    for(size_t i = 0; i < m.x1list.size(); i++){
        os << std::setprecision(15) << m.x1list[i].x << " " << m.x1list[i].y << " "  << m.x1list[i].z << " " << std::endl;
    }
    os << "##end normals list" << std::endl;
    
    os << "##start curvature list" << std::endl;
    for(size_t i = 0; i < m.x1list.size(); i++){
        os << std::setprecision(15) << m.k1_list[i] << " " << m.k2_list[i] << std::endl;
    }  
    os << "##end curvature list" << std::endl;

    os << std::endl;
    
}


//write the mesh in Moller–Trumbore format
void write_normcurv_bin_mt(std::ostream &os,const  mesh &m){
    os << "#number of vertices"<< std::endl;
    os << m.vlist.size() <<std::endl;
    //create a packed list of normals
    {
        thrust::host_vector<vertex4f> packed_normals_k1(m.nlist.size());
        for(size_t i =0; i < m.nlist.size(); i++){
            packed_normals_k1[i].x = m.nlist[i].x;
            packed_normals_k1[i].y = m.nlist[i].y;            
            packed_normals_k1[i].z = m.nlist[i].z;
            packed_normals_k1[i].w = m.k1_list[i];
        }
            
        vertex4f * packed_normals_k1_ptr = thrust::raw_pointer_cast(packed_normals_k1.data());
        os.write(reinterpret_cast<char*>(packed_normals_k1_ptr), sizeof(packed_normals_k1[0])*packed_normals_k1.size());
    }

    {
        thrust::host_vector<vertex4f> packed_x1_k2(m.x1list.size());
        for(size_t i =0; i < m.x1list.size(); i++){
            packed_x1_k2[i].x = m.x1list[i].x;
            packed_x1_k2[i].y = m.x1list[i].y;            
            packed_x1_k2[i].z = m.x1list[i].z;
            packed_x1_k2[i].w = m.k2_list[i];
        }
            
        vertex4f * packed_x1_k2_ptr = thrust::raw_pointer_cast(packed_x1_k2.data());
        os.write(reinterpret_cast<char*>(packed_x1_k2_ptr), sizeof(packed_x1_k2[0])*packed_x1_k2.size());
    }
    

    
    
    os << std::endl;
    
}

