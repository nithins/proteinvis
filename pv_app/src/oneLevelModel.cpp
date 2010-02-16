/***************************************************************************
 *   Copyright (C) 2009 by Nithin Shivashankar,   *
 *   nithin@gauss   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <cfloat>

#include <QFile>
#include <QDir>
#include <QResource>

#include <GL/glew.h>

#include <glutils.h>
#include <GLSLProgram.h>

#include <oneLevelModel.h>


using namespace std;

onelevel_model_t::onelevel_model_t ( )
{

  if ( s_sphereShader == NULL )
  {

    string sphere_log;

    QFile sp_vert ( ":/shaders/sphere.vert" );
    QFile sp_geom ( ":/shaders/sphere.geom" );
    QFile sp_frag ( ":/shaders/sphere.frag" );

    sp_vert.open ( QIODevice::ReadOnly );
    sp_geom.open ( QIODevice::ReadOnly );
    sp_frag.open ( QIODevice::ReadOnly );

    s_sphereShader = GLSLProgram::createFromSourceStrings
                     (
                       string ( sp_vert.readAll().constData() ),
                       string ( sp_geom.readAll().constData() ),
                       string ( sp_frag.readAll().constData() ),
                       GL_POINTS,
                       GL_QUADS
                     );

    sp_vert.close();
    sp_geom.close();
    sp_frag.close();

    s_sphereShader->GetProgramLog ( sphere_log );

    _LOG_VAR ( sphere_log );
  }

  if ( s_cylinderShader == NULL )
  {
    string cyl_log;

    QFile cyl_vert ( ":/shaders/cylinder.vert" );
    QFile cyl_geom ( ":/shaders/cylinder.geom" );
    QFile cyl_frag ( ":/shaders/cylinder.frag" );

    cyl_vert.open ( QIODevice::ReadOnly );
    cyl_geom.open ( QIODevice::ReadOnly );
    cyl_frag.open ( QIODevice::ReadOnly );

    s_cylinderShader = GLSLProgram::createFromSourceStrings
                       (
                         string ( cyl_vert.readAll().constData() ),
                         string ( cyl_geom.readAll().constData() ),
                         string ( cyl_frag.readAll().constData() ),
                         GL_LINES,
                         GL_QUADS
                       );

    cyl_vert.close();
    cyl_geom.close();
    cyl_frag.close();

    s_cylinderShader->GetProgramLog ( cyl_log );

    _LOG_VAR ( cyl_log );
  }
}

onelevel_model_t::~onelevel_model_t ()
{
}

void onelevel_model_t::render_spheres
( glutils::bufobj_ptr_t center_coord ,
  glutils::bufobj_ptr_t sphere_color,
  glutils::bufobj_ptr_t sphere_radii ,
  const double & add_radius,
  const double & alpha
) const
{

  glPushAttrib ( GL_ENABLE_BIT );

  s_sphereShader->use();

  s_sphereShader->sendUniform ( "ug_add_radius", ( float ) add_radius );

  s_sphereShader->sendUniform ( "ug_alpha_value", ( float ) alpha );

  GLuint radius_attrib = s_sphereShader->getAttributeLocation ( "radius" );

  sphere_radii->bind_to_vertex_attrib_pointer ( radius_attrib );

  sphere_color->bind_to_color_pointer();

  center_coord->bind_to_vertex_pointer();

  glDrawArrays ( GL_POINTS, 0, center_coord->get_num_items() );

  center_coord->unbind_from_vertex_pointer();

  sphere_color->unbind_from_color_pointer();

  sphere_radii->unbind_from_vertex_attrib_pointer ( radius_attrib );

  s_sphereShader->disable();

  glPopAttrib();
}

void onelevel_model_t::render_spheres
( glutils::bufobj_ptr_t center_coord,
  glutils::bufobj_ptr_t sphere_color,
  glutils::bufobj_ptr_t sphere_radii,
  glutils::bufobj_ptr_t sphere_indxs,
  const double & add_radius,
  const double & alpha
) const
{

  glPushAttrib ( GL_ENABLE_BIT );

  s_sphereShader->use();

  s_sphereShader->sendUniform ( "ug_add_radius", ( float ) add_radius );

  s_sphereShader->sendUniform ( "ug_alpha_value", ( float ) alpha );

  GLuint radius_attrib = s_sphereShader->getAttributeLocation ( "radius" );

  sphere_radii->bind_to_vertex_attrib_pointer ( radius_attrib );

  sphere_color->bind_to_color_pointer();

  center_coord->bind_to_vertex_pointer();

  glBindBuffer ( sphere_indxs->target(), sphere_indxs->id() );

  glDrawElements ( GL_POINTS, sphere_indxs->get_num_items(), sphere_indxs->src_type(), 0 );

  glBindBuffer ( sphere_indxs->target(), 0 );

  center_coord->unbind_from_vertex_pointer();

  sphere_color->unbind_from_color_pointer();

  sphere_radii->unbind_from_vertex_attrib_pointer ( radius_attrib );

  s_sphereShader->disable();

  glPopAttrib();
}

void onelevel_model_t::render_spheres
( glutils::bufobj_ptr_t center_coord ,
  glutils::bufobj_ptr_t sphere_color,
  double sphere_radius ) const

{

  glPushAttrib ( GL_ENABLE_BIT );

  s_sphereShader->use();

  s_sphereShader->sendUniform ( "ug_add_radius", ( float ) 0.0 );

  s_sphereShader->sendUniform ( "ug_alpha_value", ( float ) 0.0 );

  GLuint radius_attrib = s_sphereShader->getAttributeLocation ( "radius" );

  glVertexAttrib1f ( radius_attrib, ( float ) sphere_radius );

  sphere_color->bind_to_color_pointer();

  center_coord->bind_to_vertex_pointer();

  glDrawArrays ( GL_POINTS, 0, center_coord->get_num_items() );

  center_coord->unbind_from_vertex_pointer();

  sphere_color->unbind_from_color_pointer();

  s_sphereShader->disable();

  glPopAttrib();
}

void onelevel_model_t::render_spheres
( glutils::bufobj_ptr_t center_coord ,
  glutils::bufobj_ptr_t sphere_color,
  double sphere_radius,
  glutils::bufobj_ptr_t sphere_indxs ) const

{

  glPushAttrib ( GL_ENABLE_BIT );

  s_sphereShader->use();

  s_sphereShader->sendUniform ( "ug_add_radius", ( float ) 0.0 );

  s_sphereShader->sendUniform ( "ug_alpha_value", ( float ) 0.0 );

  GLuint radius_attrib = s_sphereShader->getAttributeLocation ( "radius" );

  glVertexAttrib1f ( radius_attrib, ( float ) sphere_radius );

  sphere_color->bind_to_color_pointer();

  center_coord->bind_to_vertex_pointer();

  glBindBuffer ( sphere_indxs->target(), sphere_indxs->id() );

  glDrawElements ( GL_POINTS, sphere_indxs->get_num_items(), sphere_indxs->src_type(), 0 );

  glBindBuffer ( sphere_indxs->target(), 0 );

  center_coord->unbind_from_vertex_pointer();

  sphere_color->unbind_from_color_pointer();

  s_sphereShader->disable();

  glPopAttrib();
}

void onelevel_model_t::render_cylinders
( glutils::bufobj_ptr_t endpts_coord,
  glutils::bufobj_ptr_t endpts_indxs,
  double cyl_radius ) const
{

  glPushAttrib ( GL_ENABLE_BIT );

  s_cylinderShader->use();

  s_cylinderShader->sendUniform ( "ug_cylinder_radius", ( float ) cyl_radius );

  endpts_coord->bind_to_vertex_pointer();

  glBindBuffer ( endpts_indxs->target(), endpts_indxs->id() );

  glDrawElements ( GL_LINES, endpts_indxs->get_num_items() *2, endpts_indxs->src_type(), 0 );

  glBindBuffer ( endpts_indxs->target(), 0 );

  endpts_coord->unbind_from_vertex_pointer();

  s_cylinderShader->disable();

  glPopAttrib();

}

void onelevel_model_t::render_sphere
( const double &x,
  const double &y,
  const double &z,
  const double &r
)
{
  s_sphereShader->use();

  GLuint radius_attrib = s_sphereShader->getAttributeLocation ( "radius" );

  glBegin ( GL_POINTS );
  glVertexAttrib1f ( radius_attrib, ( float ) r );
  glVertex3d ( x, y, z );
  glEnd();

  s_sphereShader->disable();
}

void onelevel_model_t::render_cylinder
( const double &x1,
  const double &y1,
  const double &z1,
  const double &x2,
  const double &y2,
  const double &z2,
  const double &radius
)
{

  s_cylinderShader->use();

  s_cylinderShader->sendUniform ( "ug_cylinder_radius", ( float ) radius );

  glColor3ub ( 128, 128, 128 );

  glBegin ( GL_LINES );
  glVertex3d ( x1, y1, z1 );
  glVertex3d ( x2, y2, z2 );
  glEnd();

  s_cylinderShader->disable();
}

GLSLProgram*  onelevel_model_t::s_cylinderShader = NULL;
GLSLProgram*  onelevel_model_t::s_sphereShader   = NULL;





