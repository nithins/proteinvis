#version 150

varying in float indices;
varying out float ids;
varying in vec3 colorsInp;
varying out vec3 colorsOut;
varying in vec4 position;
void main()
{
  ids=indices;
  colorsOut=colorsInp;
  gl_Position   = position;
}
