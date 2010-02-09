attribute in float radius;
varying out float g_radius;

void main()
{
  gl_FrontColor = gl_Color;
  g_radius      = radius;
  gl_Position   = gl_Vertex;
}
