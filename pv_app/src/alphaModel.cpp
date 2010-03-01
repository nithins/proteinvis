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


#include <fstream>
#include <string>
#include <sstream>
#include <cfloat>
#include <set>

#include <logutil.h>
#include <cpputils.h>

#include <protein.h>
#include <alphaModel.h>

using namespace std;

typedef three_tuple_t<unsigned int,false> tri_t;
typedef two_tuple_t<unsigned int,false> ege_t;

typedef std::set < tri_t> tri_set_t;
typedef std::set < ege_t> ege_set_t;

const char * g_num_tets_line_prefix = "REMARK  Number of tetrahedra :";

const char * g_tet_line_prefix = "Tetra";

const char * g_num_tris_line_prefix = "REMARK  Number of triangles :";

const char * g_tri_line_prefix = "Trig";

const char * g_num_eges_line_prefix = "REMARK  Number of edges :";

const char * g_ege_line_prefix = "Edge";

void alpha_complex_model_t::read_file ( const char * filename )
{
  fstream alpFile ( filename, ios::in );

  // get the number of tetras

  string num_tet_pfx ( g_num_tets_line_prefix );

  string tet_pfx ( g_tet_line_prefix );

  uint num_tets = 0;

  while ( !alpFile.eof() )
  {
    string line;

    getline ( alpFile, line );

    if ( string ( line.begin(), line.begin() + num_tet_pfx.size() ) == num_tet_pfx )
    {
      stringstream linestream ( string ( line.begin() + num_tet_pfx.size(), line.end() ) );

      linestream >> num_tets;

      break;
    }
  }

  uint num_tets_read = 0;

  uint * tets = new uint[num_tets*4];

  while ( !alpFile.eof() )
  {
    if ( num_tets == num_tets_read )
      break;

    string line;

    getline ( alpFile, line );

    vector<string> tokens;

    tokenize_string ( line, tokens );

    if ( tokens[0] != tet_pfx )
      continue;

    tets[4*num_tets_read+0] = atoi ( tokens[2].c_str() ) - 1;

    tets[4*num_tets_read+1] = atoi ( tokens[3].c_str() ) - 1;

    tets[4*num_tets_read+2] = atoi ( tokens[4].c_str() ) - 1;

    tets[4*num_tets_read+3] = atoi ( tokens[5].c_str() ) - 1;

    int flip = atoi ( tokens[6].c_str() );

    if ( flip == 1 )
      swap ( tets[4*num_tets_read+0], tets[4*num_tets_read+1] );

    tokens.clear();

    ++num_tets_read;
  }


  // get the number of triangles

  string num_t_pfx ( g_num_tris_line_prefix );

  string trig_pfx ( g_tri_line_prefix );

  uint num_tris = 0;

  while ( !alpFile.eof() )
  {
    string line;

    getline ( alpFile, line );

    if ( string ( line.begin(), line.begin() + num_t_pfx.size() ) == num_t_pfx )
    {
      stringstream linestream ( string ( line.begin() + num_t_pfx.size(), line.end() ) );

      linestream >> num_tris;

      break;
    }
  }

  uint num_tris_read = 0;

  tri_set_t tri_set;

  while ( !alpFile.eof() )
  {
    if ( num_tris == num_tris_read )
      break;

    string line;

    getline ( alpFile, line );

    vector<string> tokens;

    tokenize_string ( line, tokens );

    if ( tokens[0] != trig_pfx )
      continue;

    uint v1 = atoi ( tokens[2].c_str() ) - 1;

    uint v2 = atoi ( tokens[3].c_str() ) - 1;

    uint v3 = atoi ( tokens[4].c_str() ) - 1;

    tri_set.insert ( tri_t ( v1, v2, v3 ) );

    tokens.clear();

    ++num_tris_read;
  }

  // get the number of edges

  string num_e_pfx ( g_num_eges_line_prefix );

  string ege_pfx ( g_ege_line_prefix );

  uint num_eges = 0;

  while ( !alpFile.eof() )
  {
    string line;

    getline ( alpFile, line );

    if ( string ( line.begin(), line.begin() + num_e_pfx.size() ) == num_e_pfx )
    {
      stringstream linestream
          ( string ( line.begin() + num_e_pfx.size(), line.end() ) );

      linestream >> num_eges;

      break;
    }

  }

  ege_set_t ege_set;

  uint num_eges_read = 0;

  while ( !alpFile.eof() )
  {
    if ( num_eges == num_eges_read )
      break;

    string line;

    getline ( alpFile, line );

    vector<string> tokens;

    tokenize_string ( line, tokens );

    if ( tokens[0] != ege_pfx )
      continue;

    uint v1 = atoi ( tokens[2].c_str() ) - 1;

    uint v2 = atoi ( tokens[3].c_str() ) - 1;

    ege_set.insert ( ege_t( v1, v2 ) );

    tokens.clear();

    ++num_eges_read;
  }

  alpFile.close();

  if ( num_tets_read != num_tets )
    throw genericException ( "Not all Tetras were read " );

  if ( num_tris_read != num_tris )
    throw genericException ( "Not all trigs were read " );

  if ( num_eges_read != num_eges )
    throw genericException ( "Not all edges were read" );

  // throw away all edges present in the tri set

  for ( tri_set_t::iterator it = tri_set.begin();it != tri_set.end() ;++it )
  {
    uint v1 = (*it)[0];
    uint v2 = (*it)[1];
    uint v3 = (*it)[2];

    ege_set.erase ( ege_t( v1, v2 ) );
    ege_set.erase ( ege_t( v2, v3 ) );
    ege_set.erase ( ege_t( v3, v1 ) );
  }

  num_eges = ege_set.size();

  uint * eges = new uint[num_eges *2];

  uint num_eges_added = 0;

  for ( ege_set_t::iterator it = ege_set.begin();it != ege_set.end();++it, ++num_eges_added )
  {
    eges[2*num_eges_added+0] = (*it)[0];
    eges[2*num_eges_added+1] = (*it)[1];
  }

  // throw away all triangles present in the tetra set

  for ( uint i = 0 ; i < num_tets; ++i )
  {
    uint v1 = tets[4*i+0];
    uint v2 = tets[4*i+1];
    uint v3 = tets[4*i+2];
    uint v4 = tets[4*i+3];

    tri_set.erase ( tri_t( v2, v3, v4 ) );
    tri_set.erase ( tri_t( v1, v3, v4 ) );
    tri_set.erase ( tri_t( v1, v2, v4 ) );
    tri_set.erase ( tri_t( v1, v3, v3 ) );
  }

  num_tris = tri_set.size();

  uint * tris = new uint[num_tris *3];

  uint num_tris_added = 0;

  for ( tri_set_t::iterator it = tri_set.begin();it != tri_set.end();++it, ++num_tris_added )
  {
    tris[3*num_tris_added+0] = (*it)[0];
    tris[3*num_tris_added+1] = (*it)[1];
    tris[3*num_tris_added+2] = (*it)[2];
  }

  glutils::bufobj_ptr_t tet_bo,tri_bo,ege_bo,vrt_bo,col_bo;

  tet_bo = glutils::buf_obj_t::create_bo( tets, GL_UNSIGNED_INT, 4,
                                          GL_ELEMENT_ARRAY_BUFFER,
                                          sizeof ( uint )  *4*num_tets, 0 );

  tri_bo = glutils::buf_obj_t::create_bo( tris, GL_UNSIGNED_INT, 3,
                                          GL_ELEMENT_ARRAY_BUFFER,
                                          sizeof ( uint )  *3*num_tris, 0 );

  ege_bo = glutils::buf_obj_t::create_bo( tris, GL_UNSIGNED_INT, 2,
                                          GL_ELEMENT_ARRAY_BUFFER,
                                          sizeof ( uint )  *2*num_eges, 0 );

  vrt_bo = m_protein_rd->get_coord_bo();

  col_bo = glutils::buf_obj_t::create_bo();

  m_tet_ren = glutils::create_buffered_flat_tetrahedrons_ren ( vrt_bo, tet_bo, col_bo );

  m_tri_ren = glutils::create_buffered_flat_triangles_ren ( vrt_bo, tri_bo, col_bo );

  m_ege_ren = glutils::create_buffered_lines_ren ( vrt_bo, ege_bo, col_bo );

  delete []tets;

  delete []eges;

  delete []tris;
}

alpha_complex_model_t::alpha_complex_model_t
(
  const char *filename,
  protein_rd_t * protein_rd
)    :
    m_protein_rd ( protein_rd )
{
  read_file ( filename );
}

alpha_complex_model_t::~alpha_complex_model_t ()
{
  delete m_tet_ren;
  delete m_tri_ren;
  delete m_ege_ren;
}

int alpha_complex_model_t::render_tetrahedrons() const
{
  glColor3f ( 0.8, 0.2, 0.8 );

  glPushAttrib ( GL_ENABLE_BIT );

  glEnable ( GL_RESCALE_NORMAL );

  glDisable ( GL_CULL_FACE );

  m_tet_ren->render();

  glPopAttrib ();

  return 0;
}

int alpha_complex_model_t::render_triangles() const
{
  glColor3f ( 0.8, 0.2, 0.8 );

  glPushAttrib ( GL_ENABLE_BIT );

  glEnable ( GL_RESCALE_NORMAL );

  glDisable ( GL_CULL_FACE );

  m_tri_ren->render();

  glPopAttrib ();

  return 0;
}

int alpha_complex_model_t::render_edges() const
{
  glColor3f ( 0.8, 0.2, 0.8 );

  glPushAttrib ( GL_ENABLE_BIT );

  glEnable ( GL_RESCALE_NORMAL );

  m_ege_ren->render();

  glPopAttrib ();

  return 0;
}
