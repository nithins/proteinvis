#version 330 compatibility

in vec3  x_dir;
in vec3  q;
in vec3  r;
in vec3  mc_pos;

const float pitch  = 5.4;
const float width  = 0.6;

#define NUM_LIGHTS 2
#define PI     3.14159265358979323846264338
#define TWOPI (2*PI)

vec4 perform_lighting(vec4 color,vec3 pt,vec3 N)
{
  vec4 amb  = vec4(0.0);
  vec4 diff = vec4(0.0);
  vec4 spec = vec4(0.0);

  
  for(int i = 0 ; i < NUM_LIGHTS; ++i)
  {
    if(gl_LightSource[i].position.w == 0.0)
    {
      
      float nDotVP = max(0.0,dot(N,normalize(vec3(gl_LightSource[i].position))));
      float nDotHV = max(0.0,dot(N,normalize(vec3(gl_LightSource[i].halfVector))));

      float pf = 0.0;
      
      if(nDotVP != 0.0)
        pf = pow(nDotHV,gl_FrontMaterial.shininess);
      
      amb   += gl_LightSource[i].ambient;
      diff  += gl_LightSource[i].diffuse*nDotVP;
      spec  += gl_LightSource[i].specular*pf;

    }
      
    else
    {
      vec3 VP = vec3(gl_LightSource[i].position)-pt;
  
      float d  = length(VP);
      
      VP = normalize(VP);
      
      float attenuation = 1.0/(gl_LightSource[i].constantAttenuation +
                gl_LightSource[i].linearAttenuation*d+
                gl_LightSource[i].quadraticAttenuation*d*d);
                
      vec3 halfVector = normalize(VP);
      
      float nDotVP = max(0.0,dot(N,VP));
      float nDotHV = max(0.0,dot(N,halfVector));

      float pf = 0.0;
      
      if(nDotVP != 0.0)
        pf = pow(nDotHV,gl_FrontMaterial.shininess);
      
      amb   += gl_LightSource[i].ambient*attenuation;
      diff  += gl_LightSource[i].diffuse*nDotVP*attenuation;
      spec  += gl_LightSource[i].specular*pf*attenuation;  

    }
  }
  
  vec4 colorout = vec4(0);
  
  colorout +=  amb*color;
  colorout += diff*color;
  colorout += spec*gl_FrontMaterial.specular;

  return vec4(colorout.xyz,1);

}


bool ray_cylinder_ixn(vec3 c, vec3 cdir, float R, vec3 l, vec3 ldir, 
                      inout vec3 pnear, inout vec3 pfar)
{

  // solve this quadratic eqn for t
  // 
  // ||                 ldir.cdir                         (l-c).cdir     ||2     2
  // || [ ldir - cdir ---------------- ] t +  l-c - cdir --------------- ||  = r
  // ||                 cdir.cdir                         cdir.cdir      ||


  vec3 u = ldir  - cdir*dot(ldir ,cdir)/dot(cdir,cdir);
  vec3 v = (l-c) - cdir*dot((l-c),cdir)/dot(cdir,cdir); 


  float _a = dot(u,u);
  float _b = 2*dot(u,v);
  float _c = dot(v,v) - R*R;

  float _d = _b*_b - 4*_a*_c;

  if (_d <0 ) 
    return false;

  _d =sqrt(_d);

  float tnear = (-_b - _d)/(2*_a);
  float tfar  = (-_b + _d)/(2*_a);

  pnear = l + tnear*ldir;
  pfar  = l + tfar*ldir;

  return true;
}


// line: l(t) = u + t(v-u) .. 
// returns the coeff (t) of the point on l to which p is closest
float pt_line_proj_coeff(vec3 u,vec3 v,vec3 p)
{
  return dot(p-u,v-u)/dot(v-u,v-u);
}

void main()
{
  float radius = length(x_dir);
  vec3    e    = (gl_ModelViewMatrixInverse*vec4(0,0,0,1)).xyz;
  vec3 edir    = mc_pos-e;

  vec3 pnear   =vec3(0,0,0),pfar=vec3(0,0,0);

  if(ray_cylinder_ixn(q,r-q,radius,e,edir,pnear,pfar) == false)
    discard;

  vec3 qr    =  normalize(r-q);
  vec3 y_dir =  cross(qr,x_dir);
  
  
  vec3 pt       = pnear;
  bool isFar    = false;
  float t       = pt_line_proj_coeff(q,r,pt);  
  vec3 apt      = q + t*(r-q);
  vec3 normal   = normalize(pt - apt);
  float theta   = TWOPI*length(apt-q)/pitch;
  vec3  hdir    = (x_dir*cos(theta) + y_dir*sin(theta))/radius;

  
  if( t <0.0 || t >= 1.0 || dot(hdir,normal)  < 1.0-width)
  {
    pt       = pfar;    
    isFar    = true;
    t        = pt_line_proj_coeff(q,r,pt);
    apt      = q + t*(r-q);
    normal   = normalize(pt - apt);
    theta    = TWOPI*length(apt-q)/pitch;
    hdir     = (x_dir*cos(theta) + y_dir*sin(theta))/radius;    
  }
    
  if( t <0.0 || t >= 1.0 || dot(hdir,normal)  < 1.0-width)
    discard;
    
  if(isFar)
    normal *=-1;
  
  normal        = normalize(gl_NormalMatrix*normal);

  vec3 color   = gl_Color.xyz;

  pt = (gl_ModelViewMatrix*vec4(pt,1)).xyz;
  gl_FragColor   = perform_lighting(vec4(color,1),pt,normal);

  vec4  dc_pt     = gl_ProjectionMatrix*vec4(pt,1); 
  gl_FragDepth    = (dc_pt.z/dc_pt.w +1.0)/2.0;
}

