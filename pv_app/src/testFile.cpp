/***************************************************************************
 *   Copyright (C) 2012 by Pranav D Bagur,
 *
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

#include <fstream>
#include <string>
#include <sstream>
#include <set>
#include <GL/glew.h>
#include <boost/bind.hpp>

#include <logutil.h>
#include <cpputils.h>

#include <protein.h>
#include <testFile.h>
#include <malloc.h>
#include <math.h>

using namespace glutils;

test_model_t::test_model_t()
{
  InitShaders();

  vertex_list_t tubeVerts;


  tubeVerts.push_back(vertex_t(0,0,0));
  tubeVerts.push_back(vertex_t(1,0,0));
  tubeVerts.push_back(vertex_t(1,1,0));
  tubeVerts.push_back(vertex_t(0,1,0));
  tubeVerts.push_back(vertex_t(0,0,0));
  tubeVerts.push_back(vertex_t(1,0,0));
  tubeVerts.push_back(vertex_t(1,1,0));




  m_tube_Pts_bo=glutils::make_buf_obj(tubeVerts);

  _LOG_VAR(m_tube_Pts_bo->get_num_items());

}

void test_model_t::InitShaders()
{
  //initialize tubes shader
  string test_log;

  QFile cyl_vert ( "/home/nithin/projects/proteinvis/pv_app/resources/cylinder_adj_vert.glsl" );
  QFile cyl_geom ( "/home/nithin/projects/proteinvis/pv_app/resources/cylinder_adj_geom.glsl" );
  QFile cyl_frag ( "/home/nithin/projects/proteinvis/pv_app/resources/cylinder_adj_frag.glsl" );


  cyl_vert.open ( QIODevice::ReadOnly );
  cyl_geom.open ( QIODevice::ReadOnly );
  cyl_frag.open ( QIODevice::ReadOnly );

  s_tubeShader = GLSLProgram::createFromSourceStrings
      (
        string ( cyl_vert.readAll().constData() ),
        string ( cyl_geom.readAll().constData() ),
        string ( cyl_frag.readAll().constData() ),
        GL_LINE_STRIP_ADJACENCY,
        GL_TRIANGLE_STRIP
        );

  cyl_vert.close();
  cyl_geom.close();
  cyl_frag.close();

  s_tubeShader->GetProgramLog ( test_log );

  _LOG_VAR ( test_log );
}

void test_model_t::Render()
{

  glPushMatrix();
  glScalef(10,10,10);
  glPushAttrib ( GL_ENABLE_BIT );
  glDisable(GL_LIGHTING);
  glColor3f(0.2,0.2,0.2);

  m_tube_Pts_bo->bind_to_vertex_pointer();
  glDrawArrays(GL_LINE_STRIP,0,m_tube_Pts_bo->get_num_items());
  m_tube_Pts_bo->unbind_from_vertex_pointer();


  s_tubeShader->use();
  m_tube_Pts_bo->bind_to_vertex_pointer();
  glDrawArrays(GL_LINE_STRIP_ADJACENCY,0,m_tube_Pts_bo->get_num_items());
  m_tube_Pts_bo->unbind_from_vertex_pointer();
  s_tubeShader->disable();

  glPopAttrib();
  glPopMatrix();
}

