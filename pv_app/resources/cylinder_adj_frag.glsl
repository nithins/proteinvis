#version 330 compatibility

in vec3  p;
in vec3  q;
in vec3  r;
in vec3  s;
in vec3  mc_pos;

const float radius = 0.2;

float side_of_plane(vec3 pt,vec3 n,vec3 ppt)
{
  return dot(n,pt) - dot(n,ppt);
}

#define NUM_LIGHTS 2

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

vec3 closest_line_pt(vec3 l,vec3 ldir,vec3 p)
{
  return l+ldir*dot(p-l,ldir)/dot(ldir,ldir);
}

const float plane_shift_eps = 0.0001;

void main()
{
  vec3    e = (gl_ModelViewMatrixInverse*vec4(0,0,0,1)).xyz;
  vec3 edir = mc_pos-e;

  vec3 pnear,pfar;

  if(ray_cylinder_ixn(q,r-q,radius,e,edir,pnear,pfar) == false)
    discard;

  vec3 pt = pnear;

  vec3 qr = normalize(r-q);

  vec3 pqr = (normalize(q-p) + qr)/2;
  vec3 qrs = (qr + normalize(s-r))/2;

  if(side_of_plane(pt,pqr,q-qr*plane_shift_eps) < 0 || 
     side_of_plane(pt,qrs,r+qr*plane_shift_eps) > 0 )
    discard;

  vec3 normal  = normalize(gl_NormalMatrix*(pt - closest_line_pt(q,r-q,pt)));
  vec3 color   = gl_Color.xyz;

  pt = (gl_ModelViewMatrix*vec4(pt,1)).xyz;
  gl_FragColor   = perform_lighting(vec4(color,1),pt,normal);

  vec4  dc_pt     = gl_ProjectionMatrix*vec4(pt,1); 
  gl_FragDepth    = (dc_pt.z/dc_pt.w +1.0)/2.0;
}

