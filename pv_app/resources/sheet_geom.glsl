#version 330 compatibility
#extension GL_EXT_gpu_shader4: enable
#extension GL_EXT_geometry_shader4: enable

//#define ENABLE_TIPS

layout( lines_adjacency ) in;
#ifndef ENABLE_TIPS
layout( triangle_strip, max_vertices=16 ) out;
#else
layout( triangle_strip, max_vertices=24 ) out;
#endif

in vec3  NormalIn[];
in float WidthIn[];

out vec3 normal;
out vec3 wc_pos;

const float height = 0.41;

void draw_quad(vec3 a, vec3 b,vec3 c,vec3 d)
{
  wc_pos      = (gl_ModelViewMatrix*vec4(a,1)).xyz; 
  gl_Position = gl_ProjectionMatrix * vec4(wc_pos,1); EmitVertex();

  wc_pos = (gl_ModelViewMatrix*vec4(b,1)).xyz; 
  gl_Position = gl_ProjectionMatrix * vec4(wc_pos,1); EmitVertex();

  wc_pos = (gl_ModelViewMatrix*vec4(c,1)).xyz; 
  gl_Position = gl_ProjectionMatrix * vec4(wc_pos,1); EmitVertex();

  wc_pos = (gl_ModelViewMatrix*vec4(d,1)).xyz; 
  gl_Position = gl_ProjectionMatrix * vec4(wc_pos,1); EmitVertex();

  EndPrimitive();
}

void draw_quad(vec3 a,vec3 an, vec3 b,vec3 bn,vec3 c,vec3 cn,vec3 d,vec3 dn)
{
  normal      = an;
  wc_pos      = (gl_ModelViewMatrix*vec4(a,1)).xyz; 
  gl_Position = gl_ProjectionMatrix * vec4(wc_pos,1); EmitVertex();

  normal      = bn;
  wc_pos = (gl_ModelViewMatrix*vec4(b,1)).xyz; 
  gl_Position = gl_ProjectionMatrix * vec4(wc_pos,1); EmitVertex();

  normal      = cn;
  wc_pos = (gl_ModelViewMatrix*vec4(c,1)).xyz; 
  gl_Position = gl_ProjectionMatrix * vec4(wc_pos,1); EmitVertex();

  normal      = dn;
  wc_pos = (gl_ModelViewMatrix*vec4(d,1)).xyz; 
  gl_Position = gl_ProjectionMatrix * vec4(wc_pos,1); EmitVertex();

  EndPrimitive();
}

void main()
{
  vec3 p = gl_PositionIn[0].xyz;
  vec3 q = gl_PositionIn[1].xyz;
  vec3 r = gl_PositionIn[2].xyz;

  vec3 pq     = q-p;
  vec3 qr     = r-q;
  
#ifdef ENABLE_TIPS
  if(length(qr) < 0.0001)
    qr = pq;
#endif 

  float pwidth =  WidthIn[0];
  float qwidth =  WidthIn[1];

  vec3 ps     = NormalIn[0]*pwidth;
  vec3 qs     = NormalIn[1]*qwidth;

  vec3 pup    = normalize(cross(pq,ps))*height;
  vec3 qup    = normalize(cross(qr,qs))*height;

  vec3 pup_wc = normalize(gl_NormalMatrix*pup);
  vec3 qup_wc = normalize(gl_NormalMatrix*qup);
  vec3 ps_wc  = normalize(gl_NormalMatrix*ps);
  vec3 qs_wc  = normalize(gl_NormalMatrix*qs);

  gl_FrontColor = gl_FrontColorIn[0];

  draw_quad(p+pup-ps,pup_wc,
            q+qup-qs,qup_wc,
            p+pup+ps,pup_wc,
            q+qup+qs,qup_wc);

  draw_quad(p-pup-ps,-pup_wc,
            p-pup+ps,-pup_wc,
            q-qup-qs,-qup_wc,
            q-qup+qs,-qup_wc);

  gl_FrontColor = vec4(0.25,0,0,1);

  draw_quad(p+pup+ps,ps_wc,
            q+qup+qs,qs_wc,
            p-pup+ps,ps_wc,
            q-qup+qs,qs_wc);

  draw_quad(p+pup-ps,-ps_wc,
            p-pup-ps,-ps_wc,
            q+qup-qs,-qs_wc,
            q-qup-qs,-qs_wc);

#ifdef ENABLE_TIPS
  gl_FrontColor = vec4(0.25,0,0,1);

  normal = normalize(gl_NormalMatrix*qr);
  
  draw_quad(q+qup+qs,
            q+qup-qs,
            q-qup+qs,
            q-qup-qs);

  normal = normalize(-gl_NormalMatrix*pq);

  draw_quad(p+pup+ps,
            p-pup+ps,
            p+pup-ps,
            p-pup-ps);

#endif
}
