#version 330 compatibility
#extension GL_EXT_gpu_shader4: enable
#extension GL_EXT_geometry_shader4: enable

layout( lines ) in;
layout( triangle_strip, max_vertices=18 ) out;

varying out vec3 upperAxisPosition;
varying out vec3 lowerAxisPosition;
varying out vec3 Position;
varying out mat4 withTransMat;
varying out mat4 onlyRotMat;


uniform float ug_cylinder_radius;


void main()
{


vec3 pos0=gl_PositionIn[0].xyz;
vec3 pos1=gl_PositionIn[1].xyz;

/*if(pos0.y>pos1.y)
{
vec3 temp=pos0;
pos0=pos1;
pos1=temp;
}*/


vec3 axisWrd   = normalize(pos1-pos0);
vec3 upVec   = vec3(0,0,1);
if(abs(dot(axisWrd,upVec)) > 0.99)
	upVec = vec3(0,1,0);
vec3 sideVec = normalize(cross(axisWrd,upVec));
upVec = normalize(cross(sideVec,axisWrd));



vec3 side[2];
side[0]=sideVec*ug_cylinder_radius;
side[1]=sideVec*ug_cylinder_radius;

vec3 up[2];
up[0]=upVec*ug_cylinder_radius;
up[1]=upVec*ug_cylinder_radius;

vec3 pos[2];
pos[0]=pos0;
pos[1]=pos1;


float[2] sign;
sign[0] = -1.0;
sign[1] = 1.0;



const int sArr[10] = int[10](1, 1, 0, 0, 0, 0, 1, 1, 1,1);
const int uArr[10] =int[10](1, 1, 1, 1, 0, 0, 0, 0, 1,1);

int index=0;
int sideSign=0;
int upSign=0;
int flag=0;
vec3 tempPosVar;

vec3 positionCameraSpace0=(gl_ModelViewMatrix* vec4(pos0,1)).xyz;
vec3 positionCameraSpace1=(gl_ModelViewMatrix* vec4(pos1,1)).xyz;

vec3 centre=((positionCameraSpace0+positionCameraSpace1)/2);
mat4 trans=transpose(mat4(1,0,0,-centre.x,0,1,0,-centre.y,0,0,1,-centre.z,0,0,0,1));


vec4 n1Trans=trans*vec4(positionCameraSpace0,1);
vec4 n2Trans=trans*vec4(positionCameraSpace1,1);
	

vec4 axis = normalize(n2Trans - n1Trans);

float angleFrRot=float(acos(axis.y / sqrt(axis.x * axis.x + axis.y * axis.y)));

if (n2Trans.x < 0)
	angleFrRot *= -1;

float cosAngle = float(cos(angleFrRot));
float sinAngle = float(sin(angleFrRot));
mat4 rotZ = transpose(mat4(cosAngle,-sinAngle,0,0,sinAngle,cosAngle,0,0,0,0,1,0,0,0,0,1));
//mat4 inverseRotZ = {cosAngle,sinAngle,0,0,-sinAngle,cosAngle,0,0,0,0,1,0,0,0,0,1};


angleFrRot = float(acos(sqrt(axis.x * axis.x + axis.y * axis.y) / length(axis)));
if (n2Trans.z > 0)
	angleFrRot *= -1;
cosAngle = float(cos(angleFrRot));
sinAngle = float(sin(angleFrRot));

mat4 rotX = transpose(mat4(1,0,0,0,0,cosAngle,-sinAngle,0,0,sinAngle,cosAngle,0,0,0,0,1));
//vec4x4 inverseRotX = {1,0,0,0,0,cosAngle,sinAngle,0,0,-sinAngle,cosAngle,0,0,0,0,1};



mat4 PositionTransMatrix=rotX*rotZ*trans;
mat4 NormalTransMatrix=rotX*rotZ;
withTransMat=PositionTransMatrix;
onlyRotMat=NormalTransMatrix;
upperAxisPosition=(PositionTransMatrix*vec4(positionCameraSpace1,1)).xyz;
lowerAxisPosition=(PositionTransMatrix*vec4(positionCameraSpace0,1)).xyz;

vec4 positionTrans;

for(int i=0;i<10;i++)
{

	index=int(i%2);
	sideSign=sArr[i];
	upSign=uArr[i];
	tempPosVar=pos[index]+sign[sideSign]*side[index]+sign[upSign]*up[index];		
	gl_Position = gl_ModelViewProjectionMatrix * vec4(tempPosVar,1);
	//Position=(gl_ModelViewMatrix* vec4(tempPosVar,1)).xyz;
	//Position=tempPosVar;	
	positionTrans=gl_ModelViewMatrix* vec4(tempPosVar,1);
	positionTrans=PositionTransMatrix*positionTrans;
	Position=positionTrans.xyz;
	EmitVertex();
}
EndPrimitive();



for(int i=0;i<2;i++)
{
	index=i;
	for(int j=0;j<4;j++)
		{
			sideSign=int(index*(1-2*int(j%2))+int(j%2));
			tempPosVar=pos[index]+sign[sideSign]*side[index]+sign[int(j/2)]*(up[index]);		
			gl_Position = gl_ModelViewProjectionMatrix * vec4(tempPosVar,1);
			//Position=(gl_ModelViewMatrix*vec4(tempPosVar,1)).xyz;
			//Position=tempPosVar;	
			positionTrans=gl_ModelViewMatrix* vec4(tempPosVar,1);
			positionTrans=PositionTransMatrix*positionTrans;
			Position=positionTrans.xyz;		
			EmitVertex();
		}
		EndPrimitive();
}


}
