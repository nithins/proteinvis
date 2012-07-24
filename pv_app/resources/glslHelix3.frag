#extension GL_ARB_gpu_shader5 : enable

varying in vec3 upperAxisPosition;
varying in vec3 lowerAxisPosition;
varying in vec3 Position;
varying in mat4 withTransMat;
varying in mat4 onlyRotMat;


uniform float ug_cylinder_radius;


#define NUM_LIGHTS 2

void DirectionalLight( in int i,
		       in vec3 normal,
		       inout vec4 ambient,
		       inout vec4 diffuse,
		       inout vec4 specular)
{
  float nDotVP;
  float nDotHV;
  float pf;
  
  nDotVP = max(0.0,dot(normal,normalize(vec3(gl_LightSource[i].position))));
  nDotHV = max(0.0,dot(normal,normalize(vec3(gl_LightSource[i].halfVector))));
  
  if(nDotVP == 0.0)
    pf = 0.0;
  else
    pf = pow(nDotHV,gl_FrontMaterial.shininess);
  
  ambient  += gl_LightSource[i].ambient;
  diffuse  += gl_LightSource[i].diffuse*nDotVP;
  specular += gl_LightSource[i].specular*pf;
}


void PointLight(in int i,
		in vec3 eye,
		in vec3 ecPosition3,
		in vec3 normal,
		inout vec4 ambient,
		inout vec4 diffuse,
		inout vec4 specular)
{
  float nDotVP;
  float nDotHV;
  float pf;
  float attenuation;
  float d;
  vec3  VP;
  vec3  halfVector;
  
  VP = vec3(gl_LightSource[i].position)-ecPosition3;
  
  d  = length(VP);
  
  VP = normalize(VP);
  
  attenuation = 1.0/(gl_LightSource[i].constantAttenuation +
		     gl_LightSource[i].linearAttenuation*d+
		     gl_LightSource[i].quadraticAttenuation*d*d);
		     
  halfVector = normalize(VP+eye);
  
  nDotVP = max(0.0,dot(normal,VP));
  nDotHV = max(0.0,dot(normal,halfVector));
  
  if(nDotVP == 0.0)
    pf = 0.0;
  else
    pf = pow(nDotHV,gl_FrontMaterial.shininess);
  
  ambient  += gl_LightSource[i].ambient*attenuation;
  diffuse  += gl_LightSource[i].diffuse*nDotVP*attenuation;
  specular += gl_LightSource[i].specular*pf*attenuation;  
}



void main()
{
float thickness=length(gl_ModelViewMatrix*vec4(4082,0,0,0));
float pitch=length(gl_ModelViewMatrix*vec4(0.9f,0,0,0));;
float halfAxisLength=length((upperAxisPosition-lowerAxisPosition)/2);
vec4 cameraPositionTrans=(withTransMat*vec4(0,0,0,1));

vec3 cameraPosition=cameraPositionTrans.xyz;


vec3 v=normalize(Position-cameraPosition);
vec3 cVec=cameraPosition;
float r= length(gl_ModelViewMatrix*vec4(ug_cylinder_radius,0,0,0));
float rinner= length(gl_ModelViewMatrix*vec4(ug_cylinder_radius-.2f,0,0,0));
vec4 Color=vec4(0,1,0,1);

float a=v.x*v.x+v.z*v.z;
	float b=2*(cVec.x*v.x+cVec.z*v.z);
	float c=cVec.x*cVec.x+cVec.z*cVec.z-r*r;
	
	
	float det=b*b-4*a*c;
	
	if(det<0)
	discard;
	
	
	float tin=(-b+sqrt(det))/(2*a);
	float tout=(-b-sqrt(det))/(2*a);
	

	float t=min(tin,tout);
	bool isFirst=true;
	vec3 position=cameraPosition+t*v;


c=cVec.x*cVec.x+cVec.z*cVec.z-rinner*rinner;
	
	
det=b*b-4*a*c;
	
	
float tinInner=(-b+sqrt(det))/(2*a);
float toutInner=(-b-sqrt(det))/(2*a);
	


bool isFront=true;
	vec3 normal;
bool isNormalSet=false;

	if(abs(position.y)>halfAxisLength)
		isFront=false;
		
	if(isFront)
	{
		float len=sqrt(length(position)*length(position)-r*r);
		vec3 onS=vec3(0,1,0)*len;
	
		if(position.y<0)
		onS=vec3(0,-1,0)*len;
	
		normal=normalize(position-onS);
	

		float helixParam=position.y/pitch;
		vec3 helixPosition=vec3(cos(helixParam),position.y,sin(helixParam));

		vec3 helixNormal=normalize(helixPosition-onS);

		float angle=acos(dot(normal,helixNormal));
		angle*=(180/3.1412);
	

		if(angle>thickness)
		{
		
			vec3 position2=cameraPosition+min(tinInner,toutInner)*v;
			len=sqrt(length(position2)*length(position2)-rinner*rinner);
			vec3 onS=vec3(0,1,0)*len;
	
			if(position2.y<0)
			onS=vec3(0,-1,0)*len;
	
			vec3 prevNorm=normal;
			normal=normalize(position2-onS);
	

			helixParam=position2.y/pitch;
			vec3 helixPosition=vec3(cos(helixParam),position2.y,sin(helixParam));

			vec3 helixNormal=normalize(helixPosition-onS);
			vec3 myCross=normalize(cross(normal,helixNormal));
			angle=acos(dot(normal,helixNormal));
			angle*=(180/3.1412);	
		
			if(angle>thickness || det<0)
			{
				isFront=false;
				normal=prevNorm;
			}		
			else
				normal=vec3(0,1,0);
		}
		
	}
	if(!isFront)
	{


			vec3 position2=cameraPosition+max(tinInner,toutInner)*v;
			float len=sqrt(length(position2)*length(position2)-rinner*rinner);
			vec3 onS=vec3(0,1,0)*len;
	
			if(position2.y<0)
			onS=vec3(0,-1,0)*len;
	
			vec3 prevNorm=normal;
			normal=normalize(position2-onS);
	

			float helixParam=position2.y/pitch;
			vec3 helixPosition=vec3(cos(helixParam),position2.y,sin(helixParam));

			vec3 helixNormal=normalize(helixPosition-onS);
			float angle=acos(dot(normal,helixNormal));
			angle*=(180/3.1412);	
		
			bool isOnInnerHelix=true;
			if(angle>thickness || det<0)
			{
				isOnInnerHelix=false;
			}		


		t=max(tin,tout);

		position=cameraPosition+t*v;

		helixParam=position.y/pitch;
		helixPosition=vec3(cos(helixParam),position.y,sin(helixParam));


		len=sqrt(length(position)*length(position)-r*r);
		onS=vec3(0,1,0)*len;
	
		if(position.y<0)
		onS=vec3(0,-1,0)*len;
	


		helixNormal=normalize(helixPosition-onS);

		normal=normalize(position-onS);
		angle=acos(dot(normal,helixNormal));
		angle*=(180/3.1412);

		Color=vec4(1,0,0,1);
		bool isOnBackFace=true;
		if(angle>thickness)
		{
			isOnBackFace=false;
		}

		if(isOnInnerHelix && !isOnBackFace)
		{
			normal=vec3(0,-1,0);
			Color=vec4(0,0,1,0);
		}

		if(!isOnInnerHelix && isOnBackFace)
		{
			normal=vec3(0,1,0);
			Color=vec4(0,0,1,0);
		}

		if(!isOnInnerHelix && !isOnBackFace)
			discard;

		if(abs(position.y)>halfAxisLength)
		discard;
	
	}



	



	vec4 positionTransformed=inverse(withTransMat)*vec4(position,1);
	position=positionTransformed.xyz;

	vec4 normalTransformation=inverse(onlyRotMat)*vec4(normal,1);
	normal=normalTransformation.xyz;
	normal=normalize(normal);



  vec4 amb  = vec4(0.0);
  vec4 diff = vec4(0.0);
  vec4 spec = vec4(0.0);




  vec4 color = vec4(0);
  vec4 inputCol=Color;
  

  for(int i = 0 ; i < NUM_LIGHTS; ++i)
  {
    if(gl_LightSource[i].position.w == 0.0)
      DirectionalLight(i,normalize(normal),amb,diff,spec);
    else
      PointLight(i,vec3(0,0,0),position,normalize(normal),amb,diff,spec);
  }
  


  
  color +=  amb*inputCol;
  color += diff*inputCol;
  color += spec*gl_FrontMaterial.specular;


  gl_FragColor   =  color;
}
