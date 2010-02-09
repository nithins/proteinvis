varying in vec4 f_diffuse,f_ambient;
varying in vec3 f_lightDir,f_halfVector;
varying in vec3 f_mc_prism_coords;
varying in mat4 f_mvm;
varying in mat4 f_mvm_inv;


void main()
{

  vec3 mc_p0  = (f_mvm_inv*vec4(0,0,0,1)).xyz;
  vec3 mc_p1  = f_mc_prism_coords;
  vec3 mc_v   = mc_p1 - mc_p0;

  float a     = dot(mc_v,mc_v);
  float b     = 2.0*(dot(mc_v,mc_p0));
  float c     = dot(mc_p0,mc_p0)-1.0;
  float d     = b*b - 4.0*a*c;

  if(d <0.0) 
  {
/*    gl_FragDepth   = gl_FragCoord.z;
    gl_FragColor   = gl_Color;
    return;*/
    discard;
  }

  d = sqrt(d);

  vec3 mc_sp_pt   = mc_p0 + ((-b - d)/(2.0*a))*mc_v;
  vec3 mc_normal  = mc_sp_pt;

  vec3  n = normalize((f_mvm*vec4(mc_normal,1.0)).xyz-(f_mvm*vec4(0,0,0,1.0)).xyz);
  float n_dot_l,n_dot_hv;
  vec4  lightcolor = f_ambient;

  n_dot_l = max(dot(n,f_lightDir),0.0);

  if (n_dot_l > 0.0) 
  {
    n_dot_hv    = max(dot(n,normalize(f_halfVector)),0.0);
    if(n_dot_hv > 0.0)
    {
      lightcolor += gl_FrontMaterial.specular * gl_LightSource[0].specular * pow(n_dot_hv,gl_FrontMaterial.shininess);
    }
    lightcolor += f_diffuse * n_dot_l;
  }

  vec4 dc_sp_pt = gl_ProjectionMatrix*f_mvm*vec4(mc_sp_pt,1.0); 

  gl_FragColor = gl_Color*lightcolor;
  gl_FragDepth = (dc_sp_pt.z/dc_sp_pt.w +1.0)/2.0 ;

}