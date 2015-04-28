#include "tris.h"
#include <iomanip>
#include <thrust/host_vector.h>



//write the mesh in Moller–Trumbore format
void write_tris_ssv_mt(std::ostream &os,const  mesh &m){
	os << "#number of vertices  number of triangles" << std::endl;
	os << m.vlist.size() << " " << m.flist.size() << std::endl;
	os << "##start vertex list" << std::endl;
	os << "# x y z" << std::endl;
	for(size_t i = 0; i < m.vlist.size(); i++){
		os << std::setprecision(15) << m.vlist[i].x << " " << m.vlist[i].y << " "  << m.vlist[i].z << std::endl;
	}	
	os << "##end vertex list" << std::endl;

	os << "##start face list" << std::endl;
	os << "# a b c mesh_index" << std::endl;
	for(size_t i = 0; i < m.flist.size(); i++){
		os << m.flist[i].a << " " << m.flist[i].b << " "  << m.flist[i].c << " " << m.flist[i].d << std::endl;
	}	
	os << "##end face list" << std::endl;
	os << std::endl;
	
}


//Queste sono le strutture dati che vengono salvata nel file binario
//Gli operatori uguale servono per effettuare le conversioni dalle strutture 
//dati utilizzate nel programma

struct __attribute__((__packed__)) vertex3f{
	float x,y,z;
	
	void operator =(const vec3d &v){
		x = v.x;
		y = v.y;
		z = v.z;
	}
};

struct __attribute__((__packed__)) face4i{
	int a,b,c,d;
	
	void operator =(const vec4i &f){
		a = f.a;
		b = f.b;
		c = f.c;
		d = f.d;		
	}
};



//write the mesh in Moller–Trumbore format
void write_tris_bin_mt(std::ostream &os,const  mesh &m){
	os << "#number of vertices  number of triangles" << std::endl;
	os << m.vlist.size() << " " << m.flist.size() << std::endl;
	//create a packed list of vertices
	thrust::host_vector<vertex3f> packed_verts(m.vlist.begin(),m.vlist.end());
		
	vertex3f * packed_verts_ptr = thrust::raw_pointer_cast(packed_verts.data());
	os.write(reinterpret_cast<char*>(packed_verts_ptr), sizeof(packed_verts[0])*packed_verts.size());
	
	thrust::host_vector<face4i> packed_faces(m.flist.begin(),m.flist.end());
	face4i * packed_faces_ptr = thrust::raw_pointer_cast(packed_faces.data());

	os.write(reinterpret_cast<char*>(packed_faces_ptr), sizeof(packed_faces[0])*packed_faces.size());

	os << std::endl;
	
}