#version 330 compatibility
#extension GL_EXT_gpu_shader4: enable


in  float Width;

out float WidthIn;
out vec3  NormalIn;

void main()
{  
  gl_FrontColor = gl_Color;
  gl_Position   = gl_Vertex;
  NormalIn      = normalize(gl_Normal - gl_Vertex.xyz);
  WidthIn       = Width;
}
