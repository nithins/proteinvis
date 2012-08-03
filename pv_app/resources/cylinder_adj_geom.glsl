#version 330 compatibility
#extension GL_EXT_gpu_shader4: enable
#extension GL_EXT_geometry_shader4 : enable

//#define ENABLE_CAPS

#ifndef ENABLE_CAPS
layout(lines_adjacency) in;
#else
layout(triangles) in;
#endif
layout(triangle_strip, max_vertices=24) out;


#ifndef ENABLE_CAPS
out vec3  p;
#endif
out vec3  q;
out vec3  r;
out vec3  s;
out vec3  mc_pos;

const float radius = 0.4;

void draw_quad(vec3 a, vec3 b,vec3 c,vec3 d)
{
  mc_pos = a;
  gl_Position  = gl_ModelViewProjectionMatrix*vec4(mc_pos,1); 
  EmitVertex();  

  mc_pos = b;
  gl_Position  = gl_ModelViewProjectionMatrix*vec4(mc_pos,1); 
  EmitVertex();  

  mc_pos = c;
  gl_Position  = gl_ModelViewProjectionMatrix*vec4(mc_pos,1); 
  EmitVertex();  

  mc_pos = d;
  gl_Position  = gl_ModelViewProjectionMatrix*vec4(mc_pos,1); 
  EmitVertex();  

  EndPrimitive();
}

vec3 line_plane_ixn(vec3 pn, vec3 pp, vec3 ld, vec3 lp)
{
  float t = (dot(pn,pp) - dot(pn,lp))/(dot(pn,ld));
  return lp + t*ld;
}

void main()
{
#ifndef ENABLE_CAPS
  p = gl_PositionIn[0].xyz;
  q = gl_PositionIn[1].xyz;
  r = gl_PositionIn[2].xyz;
  s = gl_PositionIn[3].xyz;
#else
  q = gl_PositionIn[0].xyz;
  r = gl_PositionIn[1].xyz;
  s = gl_PositionIn[2].xyz;
#endif

  vec3  qr = normalize(r-q);
  vec3   c = (q + r) /2;  
#ifdef ENABLE_CAPS
  vec3   p = q - qr;
#endif  
  vec3 pqr = (normalize(q-p) + qr)/2;
  vec3 qrs = (qr +  normalize(s-r))/2;

  vec3 oqr1 = vec3(0,0,1);

  if(dot(oqr1,qr) > 0.99)
    oqr1 = vec3(0,1,0);
  
  vec3 oqr2 = normalize(cross(qr,oqr1));
  oqr1 = normalize(cross(oqr2,qr));

  vec3 c0q0 = line_plane_ixn(pqr,q,qr,q+radius*(-oqr1-oqr2));
  vec3 c1q0 = line_plane_ixn(pqr,q,qr,q+radius*(+oqr1-oqr2));
  vec3 c1q1 = line_plane_ixn(pqr,q,qr,q+radius*(+oqr1+oqr2));
  vec3 c0q1 = line_plane_ixn(pqr,q,qr,q+radius*(-oqr1+oqr2));

  vec3 c0r0 = line_plane_ixn(qrs,r,qr,r+radius*(-oqr1-oqr2));
  vec3 c1r0 = line_plane_ixn(qrs,r,qr,r+radius*(+oqr1-oqr2));
  vec3 c1r1 = line_plane_ixn(qrs,r,qr,r+radius*(+oqr1+oqr2));
  vec3 c0r1 = line_plane_ixn(qrs,r,qr,r+radius*(-oqr1+oqr2));

  gl_FrontColor = gl_FrontColorIn[1];

  draw_quad(c1r0,c1r1,c0r0,c0r1);
  draw_quad(c1q0,c0q0,c1q1,c0q1);

  draw_quad(c0q0,c1q0,c0r0,c1r0);
  draw_quad(c0q1,c0r1,c1q1,c1r1);

  draw_quad(c0q0,c0r0,c0q1,c0r1);
  draw_quad(c1q0,c1q1,c1r0,c1r1);
}