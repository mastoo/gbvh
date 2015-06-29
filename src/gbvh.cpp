#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <iomanip>
#include <utility>
#include <map>
#include <set>
#include <cassert>
#include "vec3.h"
#include "vec4.h"
#include "aabb.h"
#include "bvh.h"
#include "tris.h"
#include "normals_and_curvatures.h"

enum INPUT_FORMAT {STL,PLY,OBJ};

struct mesh_file{
	std::string filename;
	std::vector<vec3f> vlist; //list of vertex
	std::vector<vec3i> flist; //list of faces
    std::vector<vec3f> nlist; //list of normals
    std::vector<vec3f> x1list; //list of the first principal direction
    std::vector<float> k1    ; //list of the first curvature
    std::vector<float> k2    ; //list of the first curvature
    
	INPUT_FORMAT frmt;	
    
    bool normals_available;
    bool curvature_available;
};




void print_usage();

bool compare_extension(const std::string &filename,const std::string &extension);

std::string trim(const std::string& str, const std::string& whitespace = " \t");

bool get_valid_line(std::istream &ifile,std::string &line);
				
bool load_ply( mesh_file &ply_mesh);
				
bool load_obj( mesh_file &ply_mesh);
					  
int main(int argc, char** argv){
	
	std::vector<mesh_file> mesh_list;
	mesh  global_mesh;
    bool write_normal_curvature = false;   
    bool assume_flat  = false;
	OUTPUT_FORMAT oformat = BIN;
	
	std::cerr << "gbvh " << __TIME__ << " " << __DATE__ << std::endl;
	
	if(argc <= 1){
		print_usage();
		return 2;
	}
	
	for(int i = 1; i < argc ; i++ ){
		
		std::string arg = std::string(argv[i]);
		if(arg[0] != '-'){
			mesh_file m;
			
			m.filename = std::string(argv[i]);		
			if(compare_extension(m.filename,"ply")){
				m.frmt = PLY;
			}else if(compare_extension(m.filename,"stl")){
				m.frmt = STL;
			}else if(compare_extension(m.filename,"obj")){
				m.frmt = OBJ;
			}else{
				std::cerr << "file format not recognized:  " <<m.filename << std::endl;
			}		
			mesh_list.push_back(m);			
		}else{
			if(arg == "-ssv"){
				oformat = SSV;
			}else if(arg == "-bin"){
				oformat = BIN;
			}else if(arg == "-nbin"){
				oformat = NBIN;
			}
			
            if(arg == "-flat"){
                assume_flat = true;
            }

            if(arg == "-wnc"){
                write_normal_curvature = true;
            }
            

		}
	}
 
	// load meshes
	for(size_t i = 0; i < mesh_list.size() ; i++ ){
		switch(mesh_list[i].frmt){
			case STL:{
				assert(false);
			}break;
			case PLY:{
				if(!load_ply(mesh_list[i])){
								 
					std::cerr << "Error opening file: " << mesh_list[i].filename << std::endl;
					return -1;			 
				}
			}break;
			case OBJ:{
				if(!load_obj(mesh_list[i])){
								 
					std::cerr << "Error opening file: " << mesh_list[i].filename << std::endl;
					return -1;			 
				}
			}break;
			
		}
	}
	

	if(write_normal_curvature && !assume_flat){
        for(size_t i = 0; i < mesh_list.size() ; i++ ){
            if(!mesh_list[i].normals_available){
                compute_normals(
                    mesh_list[i].vlist,
                    mesh_list[i].flist,
                    mesh_list[i].nlist
                );
            }
            
            if(!mesh_list[i].curvature_available){
                compute_curvatures(
                    mesh_list[i].vlist,
                    mesh_list[i].flist,
                    mesh_list[i].nlist,
                    mesh_list[i].x1list,
                    mesh_list[i].k1,
                    mesh_list[i].k2
                );
            }            
        }

    }
	

	
	{//unify the meshes
		size_t nverts  = 0;
		size_t nfaces  = 0; 
		//count the total number of meshes
		for(size_t i = 0; i < mesh_list.size() ; i++ ){
			nverts += mesh_list[i].vlist.size();
			nfaces += mesh_list[i].flist.size();			
		}
		
		global_mesh.vlist.resize(nverts);
		global_mesh.flist.resize(nfaces);
        
        
        if(write_normal_curvature && !assume_flat){
            global_mesh.nlist.resize(nverts);
            global_mesh.x1list.resize(nverts);
            global_mesh.k1_list.resize(nverts);
            global_mesh.k2_list.resize(nverts);
        }
        
        
		//fill the mesh
		size_t gvidx = 0;
		size_t gfidx = 0;
		for(size_t i =0; i < mesh_list.size(); i++){
    		int face_offset = gvidx;
			for(size_t j = 0; j < mesh_list[i].vlist.size(); j++){
				const vec3f v   = mesh_list[i].vlist[j];

                
				global_mesh.vlist[gvidx] = make_vec3d(v.x,v.y,v.z);
                
                if(write_normal_curvature && ! assume_flat){

                    const vec3f n   = mesh_list[i].nlist[j];
                    const vec3f x1  = mesh_list[i].x1list[j];
                    const float k1  = mesh_list[i].k1[j];
                    const float k2  = mesh_list[i].k2[j];
                    global_mesh.nlist[gvidx] = make_vec3d(n.x,n.y,n.z);
                    global_mesh.x1list[gvidx] = make_vec3d(x1.x,x1.y,x1.z);
                    
                    global_mesh.k1_list[gvidx] = k1;
                    global_mesh.k2_list[gvidx] = k2;
                
                }
				gvidx++;
			}
			
			for(size_t j = 0; j < mesh_list[i].flist.size(); j++){
				const vec3i face = mesh_list[i].flist[j];
				//set index of the mesh in the field d
				global_mesh.flist[gfidx] = make_vec4i(face.a+face_offset,
													  face.b+face_offset,
													  face.c+face_offset,i);
				gfidx++;
			}			
		}
	}

	
	global_mesh.bb = compute_aabb(global_mesh.vlist,global_mesh.flist);
	
	std::cerr << "bounding box: " << global_mesh.bb << std::endl;
	
	bvh_builder bbuild(bvh_platform(),
					   global_mesh);
	
	bbuild.build();
	
    
    if(assume_flat && write_normal_curvature){
        recompute_flat_verts_normals_and_curvature(global_mesh);
    }
    
    
	switch(oformat){
		case SSV:{
			write_tris_ssv_mt(std::cout, global_mesh);
			bbuild.write_bvh_ssv(std::cout,oformat);
            if(write_normal_curvature){
                write_normcurv_ssv_mt(std::cout, global_mesh);
            }
            
		}break;
		case BIN:case NBIN:{
            
			write_tris_bin_mt(std::cout, global_mesh);
			bbuild.write_bvh_ssv(std::cout,oformat);
            if(write_normal_curvature){
                write_normcurv_bin_mt(std::cout,global_mesh);
            }
		}break;
		
	};
	

    
    
	return 0;
}


void print_usage(){
	std::cout << "usage: gbvh [option] inputfile1 inputfile2 ..." << std::endl;
    std::cout << "\t-ssv:\t produce ascii output" << std::endl;
    std::cout << "\t-bin:\t produce binary format" << std::endl;    
    std::cout << "\t-nbin:\t produce packed binary format" << std::endl;
    std::cout << "\t-flat:\t assume flat triangles" << std::endl;
    std::cout << "\t-wnc:\t write normals and curvatures" << std::endl;

}


//TODO: this doesn't handle correctly the filepath
bool compare_extension(const std::string &filename, 
					   const std::string &extension){
						   
	return filename.substr(filename.find_last_of(".") + 1) == extension;
}

bool load_ply( mesh_file &ply_mesh ){
						
	enum DATA_FMT {ASCII,BINARY_LITTLE_ENDIAN,BINARY_BIG_ENDIAN} fmt;
    std::string line;
    std::string token;
    size_t num_verts = 0;
    size_t num_faces = 0;
    ply_mesh.normals_available = false;
    ply_mesh.curvature_available = false;
    
    std::ifstream ifile(ply_mesh.filename.c_str());
		
	if (!ifile.is_open()) {		
		return false;
	} 
    
	
	{   //read the magic number
	    get_valid_line(ifile,line);
	    std::stringstream ss(line);
	    ss >> token;
	    assert(token == "ply");
	}
	
	{	//read the format
		get_valid_line(ifile,line);
	    std::stringstream ss(line);
		ss >> token;
	    assert(token == "format");
	    ss >> token;
	    if(token == "ascii"){
			fmt = ASCII;
		}else if(token == "binary_little_endian"){
			fmt = BINARY_LITTLE_ENDIAN;
			assert(false);
		}else if(token == "binary_little_endian"){
			fmt = BINARY_BIG_ENDIAN;
			assert(false);
		}else{
			assert(false);
		}
	}
	
	//read the header
	do{
    	get_valid_line(ifile,line);
    	std::stringstream ss(line);
    	ss >> token ; 
    	
    	if(token == "element"){
            ss >> token; //read the next token
            if(token == "vertex"){
                ss >> num_verts; //read the number of vertex
                
                
            }else if(token == "face"){
                ss >> num_faces; //read the number of vertex
            }                	        	    
    	}else if(token == "property"){
            ss >> token; //read the next token
            ss >> token;
            if(token == "nx" || token == "ny" || token == "nz") {
                ply_mesh.normals_available = true;
            }else if(token == "x1x" || token == "x1y" || token == "x1z" ||
                     token == "k1"  || token == "k2"
            ){
                ply_mesh.curvature_available = true;
            }
        }
	}while(token != "end_header");
	
    if(ply_mesh.curvature_available == true && ply_mesh.normals_available == false){
        assert(false);
    }
    
	ply_mesh.vlist.resize(0);
	ply_mesh.flist.resize(0);
	ply_mesh.nlist.resize(0);
    ply_mesh.x1list.resize(0);
    ply_mesh.k1.resize(0);
    ply_mesh.k2.resize(0);
    
    
	ply_mesh.vlist.reserve(num_verts);
	ply_mesh.flist.reserve(num_faces);
    if(ply_mesh.normals_available){
        ply_mesh.nlist.reserve(num_verts);
    }
    
    if(ply_mesh.curvature_available){
        ply_mesh.x1list.reserve(num_verts);
        ply_mesh.k1.reserve(num_verts);
        ply_mesh.k2.reserve(num_verts);
    }
    
	if(fmt == ASCII){
		//read verts
		for(size_t i = 0; i< num_verts; i++){
			vec3f v;
            
			get_valid_line(ifile,line);
			std::stringstream ss(line);
			ss >> v.x >> v.y >> v.z;
			ply_mesh.vlist.push_back(v);
            
            if(ply_mesh.normals_available){
                vec3f n;
                ss >> n.x >> n.y >> n.z;
                ply_mesh.nlist.push_back(n);                
            }

            if(ply_mesh.curvature_available){
                vec3f x1;
                float k1;
                float k2;
                
                ss >> x1.x >> x1.y >> x1.z >> k1 >> k2;
                ply_mesh.x1list.push_back(x1);
                ply_mesh.k1.push_back(k1);
                ply_mesh.k2.push_back(k2);
            }
            
		}
		assert(ply_mesh.vlist.size()==num_verts);
		//read faces
		for(size_t i = 0; i< num_faces; i++){
			get_valid_line(ifile,line);
			std::stringstream ss(line);
			int list_len,a,b,c;
			ss >> list_len >> a >> b >> c;
			assert(list_len == 3);
			ply_mesh.flist.push_back(make_vec3i(a,b,c));	    
		}
		assert(ply_mesh.flist.size()==num_faces);
	}
	
	

	
	return true;
}


bool load_obj(mesh_file &ply_mesh){
						
    std::string line;
    std::string prefix;
    
    ply_mesh.normals_available = false;
    ply_mesh.curvature_available = false;
    
    std::ifstream ifile(ply_mesh.filename.c_str());
		
	if (!ifile.is_open()) {		
		return false;
	}     
	
	ply_mesh.vlist.resize(0);
	ply_mesh.flist.resize(0);
	
	while(get_valid_line(ifile,line)){
		std::stringstream ss(line);
		ss >> prefix;
		
		if(prefix == "v"){
			float x,y,z;
			ss >> x >> y >> z;   
			ply_mesh.vlist.push_back(make_vec3f(x,y,z));	         
		}else if(prefix == "f"){
			int a,b,c;
			ss >> a >> b >> c;
			ply_mesh.flist.push_back(make_vec3i(a-1,b-1,c-1));	  
		}
		
	}

	
	
	return true;
}



/***** TRIM ****/
std::string trim(const std::string& str, const std::string& whitespace){
    const size_t strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const size_t strEnd = str.find_last_not_of(whitespace);
    const size_t strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}


bool get_valid_line(std::istream &ifile,std::string &line){
    do{
        std::getline(ifile, line);
        line = trim(line); // trim leading and trailing white spaces
    }while(ifile.good() && (line[0]=='#'  || line.size()==0));
    
    return ifile.good();
}
