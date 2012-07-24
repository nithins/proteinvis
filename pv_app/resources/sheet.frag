#version 120

varying in vec4 Color;
varying in vec3 sideNormal;
varying in vec3 upNormal;
varying in vec3 Position;
varying in float EdgeNumber;

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
	
  vec4 amb  = vec4(0.0);
  vec4 diff = vec4(0.0);
  vec4 spec = vec4(0.0);



  vec4 color = vec4(0);
  vec3 Normal;
  vec4 inputCol;
  
  if(EdgeNumber >=0 && EdgeNumber<=1){
	inputCol=Color;
	Normal=upNormal;}

  if(EdgeNumber>1 && EdgeNumber<2){
	inputCol=vec4(0.2f,0.2f,0.2f,1);
	Normal=sideNormal;}

  if(EdgeNumber>2 && EdgeNumber<3){
	inputCol=Color;
	Normal=upNormal;}

  if(EdgeNumber>3 && EdgeNumber<4){
	inputCol=vec4(0.2f,0.2f,0.2f,1);
	Normal=sideNormal;}
  
  for(int i = 0 ; i < NUM_LIGHTS; ++i)
  {
    if(gl_LightSource[i].position.w == 0.0)
      DirectionalLight(i,normalize(Normal),amb,diff,spec);
    else
      PointLight(i,vec3(0,0,0),Position,normalize(Normal),amb,diff,spec);
  }
  


  
  color +=  amb*inputCol;
  color += diff*inputCol;
  color += spec*gl_FrontMaterial.specular;


  gl_FragColor   =  color;	
}
