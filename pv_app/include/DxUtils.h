#ifndef DX_UTILS_HPP
#define DX_UTILS_HPP


struct D3DXVECTOR3
{
    float x,y,z;
    D3DXVECTOR3(float x,float y,float z)
    {
        this->x=x;
        this->y=y;
        this->z=z;
    }
};

struct D3DXVECTOR4
{
    float x,y,z,w;
    D3DXVECTOR4(float x,float y,float z,float w)
    {
        this->x=x;
        this->y=y;
        this->z=z;
        this->w=w;
    }
};


struct D3DXMATRIX
{
   float a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44;
   D3DXMATRIX(float a11,float a12,float a13,float a14,float a21,float a22,float a23,float a24,float a31,float a32,float a33,float a34,float a41,float a42,float a43,float a44)
   {
       this->a11=a11;
       this->a12=a12;
       this->a13=a13;
       this->a14=a14;

       this->a21=a21;
       this->a22=a22;
       this->a23=a23;
       this->a24=a24;

       this->a31=a31;
       this->a32=a32;
       this->a33=a33;
       this->a34=a34;

       this->a41=a41;
       this->a42=a42;
       this->a43=a43;
       this->a44=a44;
   }
};


inline float D3DXVec4Dot(D3DXVECTOR4 *first,D3DXVECTOR4 *second)
{
    return first->x*second->x+first->y*second->y+first->z*second->z+first->w*second->w;
}

inline float D3DXVec3Dot(D3DXVECTOR3 *first,D3DXVECTOR3 *second)
{
    return first->x*second->x+first->y*second->y+first->z*second->z;
}


inline void mulMatrixVec(D3DXVECTOR4 *result,D3DXVECTOR4 *vector,D3DXMATRIX *matrix)
{
    /*result->x=matrix->a11*vector->x+matrix->a12*vector->y+matrix->a13*vector->z+matrix->a14*vector->w;
    result->y=matrix->a21*vector->x+matrix->a22*vector->y+matrix->a23*vector->z+matrix->a24*vector->w;
    result->z=matrix->a31*vector->x+matrix->a32*vector->y+matrix->a33*vector->z+matrix->a34*vector->w;
    result->w=matrix->a41*vector->x+matrix->a42*vector->y+matrix->a43*vector->z+matrix->a44*vector->w;*/


    result->x=matrix->a11*vector->x+matrix->a21*vector->y+matrix->a31*vector->z+matrix->a41*vector->w;
    result->y=matrix->a12*vector->x+matrix->a22*vector->y+matrix->a32*vector->z+matrix->a42*vector->w;
    result->z=matrix->a13*vector->x+matrix->a23*vector->y+matrix->a33*vector->z+matrix->a43*vector->w;
    result->w=matrix->a14*vector->x+matrix->a24*vector->y+matrix->a34*vector->z+matrix->a44*vector->w;

}

inline void Normalize(D3DXVECTOR3 *vector,D3DXVECTOR3 *result)
{
   D3DXVECTOR3 temp(0,0,0);

   float len=sqrt(vector->x*vector->x+vector->y*vector->y+vector->z*vector->z);
   temp.x=(float)vector->x/(float)len;
   temp.y=(float)vector->y/(float)len;
   temp.z=(float)vector->z/(float)len;

   *result=temp;
}

inline D3DXVECTOR3 D3DAdd(D3DXVECTOR3 *first,D3DXVECTOR3 *second)
{
    return D3DXVECTOR3(first->x+second->x,first->y+second->y,first->z+second->z);
}

inline D3DXVECTOR3 D3DSubtract(D3DXVECTOR3 *first,D3DXVECTOR3 *second)
{
    return D3DXVECTOR3(first->x-second->x,first->y-second->y,first->z-second->z);
}

inline D3DXVECTOR3 D3DXCross(D3DXVECTOR3 *first,D3DXVECTOR3 *second)
{
    float x=first->y*second->z-second->y*first->z;
    float y=first->x*second->z-second->x*first->z;
    float z=first->x*second->y-second->x*first->y;

    return D3DXVECTOR3(x,-y,z);
}

#endif
