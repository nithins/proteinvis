#version 330 compatibility
#extension GL_EXT_gpu_shader4: enable
#extension GL_EXT_geometry_shader4: enable


layout( lines_adjacency ) in;
layout( triangle_strip, max_vertices=18 ) out;


varying out vec3 sideNormal;
varying out vec3 upNormal;
varying out vec3 Position;
varying out float EdgeNumber;

void main()
{

float width=1f;
float height=0.1;

vec3 p0=gl_PositionIn[0].xyz;
vec3 p1=gl_PositionIn[1].xyz;
vec3 p2=gl_PositionIn[2].xyz;


vec3 helixUp=normalize(gl_PositionIn[3].xyz);
vec3 segmentVector=p1-p0;
vec3 nextSegmentVector=p2-p1;
vec3 helixSide1=normalize(cross(segmentVector,helixUp));
vec3 helixSide2=normalize(cross(nextSegmentVector,helixUp));

helixUp*=width;
helixSide1*=height;
helixSide2*=height;


vec3 side[2];
side[0]=helixSide1;
side[1]=helixSide2;

vec3 pos[2];
pos[0]=p0;
pos[1]=p1;


float[2] sign;
sign[0] = -1.0;
sign[1] = 1.0;



const int sArr[10] = int[10](1,1,1,1,0,0,0,0,1,1);
const int uArr[10] =int[10](0,0,1,1,1,1,0,0,0,0);

int index=0;
int sideSign=0;
int upSign=0;
vec3 tempPosVar;
mat4 MVI = transpose(inverse(gl_ModelViewMatrix));






for(int i=0;i<10;i++)
{

	index=int(i%2);
	sideSign=sArr[i];
	upSign=uArr[i];
	vec3 initNormal=sign[sideSign]*side[index]+sign[upSign]*helixUp;		
	tempPosVar=pos[index]+initNormal;
	gl_Position = gl_ModelViewProjectionMatrix * vec4(tempPosVar,1);
	Position=(gl_ModelViewMatrix*vec4(tempPosVar,1)).xyz;
	
	EdgeNumber=int(i/2);
	sideNormal=mat3(MVI)*sign[sideSign]*side[index];
	upNormal=mat3(MVI)*sign[upSign]*helixUp;
	EmitVertex();
}



}
