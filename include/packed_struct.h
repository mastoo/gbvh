#ifndef __PACKED_STRUCT__H__
#define __PACKED_STRUCT__H__



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


struct __attribute__((__packed__)) vertex4f{
    float x,y,z,w;
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


#endif
