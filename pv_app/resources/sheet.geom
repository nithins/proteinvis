#version 330 compatibility
#extension GL_EXT_gpu_shader4: enable
#extension GL_EXT_geometry_shader4: enable

layout( triangles_adjacency ) in;
layout( triangle_strip, max_vertices=14 ) out;

varying in float ids[];

#define noOfSegments 15

varying in vec3 colorsOut[];

varying out vec4 Color;
varying out vec3 sideNormal;
varying out vec3 upNormal;
varying out vec3 Position;
varying out float EdgeNumber;


void main()
{

float width=0.6;
float height=0.16f;
int vertIndex=int(ids[1]);

vec3 v12=normalize(gl_PositionIn[2].xyz);
vec3 v13=gl_PositionIn[3].xyz-gl_PositionIn[1].xyz;
vec3 v34=normalize(gl_PositionIn[4].xyz);
vec3 v35=gl_PositionIn[5].xyz-gl_PositionIn[3].xyz;

v12*=width;
v34*=width;

vec3 v34xv35=normalize(cross(v34,v35));
v34xv35*=height;
vec3 v12xv13=normalize(cross(v12,v13));
v12xv13*=height;


vec3 side[2];
side[0]=v12;
side[1]=v34;

vec3 up[2];
up[0]=v12xv13;
up[1]=v34xv35;

vec3 pos[2];
pos[0]=gl_PositionIn[1].xyz;
pos[1]=gl_PositionIn[3].xyz;


float[2] sign;
sign[0] = -1.0;
sign[1] = 1.0;

Color=vec4(colorsOut[1].xyz/256,1);


const int sArr[10] = int[10](1, 1, 0, 0, 0, 0, 1, 1, 1,1);
const int uArr[10] =int[10](1, 1, 1, 1, 0, 0, 0, 0, 1,1);

int index=0;
int sideSign=0;
int upSign=0;
int flag=0;
vec3 tempPosVar;
mat4 MVI = transpose(inverse(gl_ModelViewMatrix));

float[2] arrowScales;
arrowScales[0]=1;
arrowScales[1]=1;
if(vertIndex>=0)
{
	if(vertIndex!=0){
	arrowScales[0]=1-(float(int(vertIndex-1)%noOfSegments)/float(noOfSegments-1))+(float(1)/float(noOfSegments-1))+0.3f;
	arrowScales[1]=1-(float(int(vertIndex-1)%noOfSegments)/float(noOfSegments-1))+0.3f;
	side[0]+=(width*1.5f)*side[0];
	side[1]+=(width*1.5f)*side[1];
	flag=1;
	}
	index=0;
	EdgeNumber=6;
	if(vertIndex==1 || vertIndex==0 || vertIndex==noOfSegments)
	{
		if(vertIndex==noOfSegments)
			index=1;

		for(int j=0;j<4;j++)
		{
			sideSign=int(index*(1-2*int(j%2))+int(j%2));
			tempPosVar=pos[index]+sign[sideSign]*(side[index])*arrowScales[index]+sign[int(j/2)]*(up[index]);		
			gl_Position = gl_ModelViewProjectionMatrix * vec4(tempPosVar,1);
			Position=(gl_ModelViewMatrix*vec4(tempPosVar,1)).xyz;
			EmitVertex();
		}
		EndPrimitive();
	}
}

for(int i=0;i<10;i++)
{

	index=int(i%2);
	sideSign=sArr[i];
	upSign=uArr[i];
	tempPosVar=pos[index]+sign[sideSign]*side[index]*arrowScales[index]+sign[upSign]*up[index];		
	gl_Position = gl_ModelViewProjectionMatrix * vec4(tempPosVar,1);
	Position=(gl_ModelViewMatrix*vec4(tempPosVar,1)).xyz;	
	upNormal=mat3(MVI)*(sign[upSign]*up[index]);
	sideNormal=mat3(MVI)*(sign[sideSign]*side[index]+flag*cross(sign[sideSign]*side[index],sign[upSign]*up[index]));
	EdgeNumber=i/2;
	EmitVertex();
}

}
