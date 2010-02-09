#version 120
#extension GL_EXT_geometry_shader4 : enable

varying in float  g_radius[];

varying out vec4 f_diffuse,f_ambient;
varying out vec3 f_lightDir,f_halfVector;
varying out vec3 f_mc_prism_coords;
varying out mat4 f_mvm;
varying out mat4 f_mvm_inv;

const vec2 bin01   = vec2(-1,1);

void main()
{
  // lighting stuff 

  f_lightDir    = normalize(vec3(gl_LightSource[0].position));
  f_halfVector  = normalize(gl_LightSource[0].halfVector.xyz);
  f_diffuse     = gl_LightSource[0].diffuse;
  f_ambient     = gl_LightSource[0].ambient;
  f_ambient    += gl_LightModel.ambient ;

  

  for(int i=0; i< gl_VerticesIn; i++)
  { 
    mat4 scale_mat      = mat4(g_radius[i],0,0,0,0,g_radius[i],0,0,0,0,g_radius[i],0,0,0,0,1);
    mat4 scale_mat_inv  = mat4(1.0/g_radius[i],0,0,0,0,1.0/g_radius[i],0,0,0,0,1.0/g_radius[i],0,0,0,0,1.0);

    mat4 trans_mat      = mat4(1,0,0,0,0,1,0,0,0,0,1,0, gl_PositionIn[i].xyz,1);
    mat4 trans_mat_inv  = mat4(1,0,0,0,0,1,0,0,0,0,1,0,-gl_PositionIn[i].xyz,1);
    f_mvm               = gl_ModelViewMatrix*trans_mat*scale_mat;
    f_mvm_inv           = scale_mat_inv*trans_mat_inv*gl_ModelViewMatrixInverse;

    gl_FrontColor       = gl_FrontColorIn[i]; 

    mat3 mcCylinderAxes = mat3(1,0,0,0,1,0,0,0,1);
    int[3] valid_face_version;

    for(int j = 0 ; j < 3 ; j ++)
    {
      vec4 cp0 = f_mvm*vec4(bin01[0]*mcCylinderAxes[(j+2)%3],1.0);
      vec4 cp1 = f_mvm*vec4(bin01[1]*mcCylinderAxes[(j+2)%3],1.0);

      if(length(cp0) < length(cp1)) 
	valid_face_version[j] = 0; 
      else 
	valid_face_version[j] = 1;
    }


    for(int j = 0 ; j < 3 ; j ++)
    {
      int k = valid_face_version[j];

      vec3 r_dir = mcCylinderAxes[j];
      vec3 u_dir = bin01[k]*mcCylinderAxes[(j+1)%3];
      vec3 n_dir = bin01[k]*mcCylinderAxes[(j+2)%3];	

      for(int l = 0 ; l < 4; l++)
      {
	f_mc_prism_coords = n_dir + bin01[l%2]*r_dir + bin01[(l/2)%2]*u_dir;
	gl_Position       = gl_ProjectionMatrix*f_mvm*vec4(f_mc_prism_coords,1.0); 
	EmitVertex();
      }
      EndPrimitive();
    }
  }
}

