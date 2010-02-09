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
#include <set>

#include <boost/bind.hpp>

#include <logutil.h>
#include <cpputils.h>

#include <protein.h>
#include <pocketModel.h>


using namespace std;


template <typename list_t, typename ind_t>

class compare_list_items
{
    const list_t &list;

  public:

    compare_list_items ( const list_t &l ) : list ( l ) {}

    bool operator() ( ind_t i, ind_t j )
    {
      return list[i] < list[j];
    }
};

ostream& operator<< ( ostream &o, const pair<uint, uint> &p )
{
  o << "(" << p.first << "," << p.second << ")";
  return o;
}



void pocket_model_t::read_file ( const char * tepoc_filename, const char * tet_filename )
{

  using boost::bind;

  fstream tepocFile ( tepoc_filename, ios::in );

  vector<uint> alpha_list;

  vector<uint> tetno_list;

  vector<uint> pocno_list;

  while ( !tepocFile.eof() )
  {
    string line;

    getline ( tepocFile, line );

    vector<string> tokens;

    line = stripWS ( stripLineComments ( line ) );

    if ( line.size() == 0 )
      continue;

    tokenize_string ( line, tokens );

    alpha_list.push_back ( atoi ( tokens[0].c_str() ) );

    tetno_list.push_back ( atoi ( tokens[1].c_str() ) - 1 );

    pocno_list.push_back ( atoi ( tokens[2].c_str() ) );
  }

  fstream tetFile ( tet_filename, ios::in );

  string num_tets_line;

  getline ( tetFile, num_tets_line );

  vector<uint> tetra_list;

  while ( !tetFile.eof() )
  {
    string line;

    getline ( tetFile, line );

    line = stripWS ( stripLineComments ( line ) );

    if ( line.size() == 0 )
      continue;

    vector<string> tokens;

    tokenize_string ( line, tokens );

    tetra_list.push_back ( atoi ( tokens[1].c_str() ) - 1 );

    tetra_list.push_back ( atoi ( tokens[2].c_str() ) - 1 );

    tetra_list.push_back ( atoi ( tokens[3].c_str() ) - 1 );

    tetra_list.push_back ( atoi ( tokens[4].c_str() ) - 1 );
  }

  vector < pair<uint, uint> >  alpha_ranges;

  uint range_start = 0;

  for ( uint i = 1 ; i < alpha_list.size();i++ )
  {
    if ( alpha_list[i] != alpha_list[i-1] )
    {
      alpha_ranges.push_back ( make_pair ( range_start, i ) );
      range_start = i;
    }
  }

  alpha_ranges.push_back ( make_pair ( range_start, alpha_list.size() ) );

  vector<uint> pocno_remapping;

  pocno_remapping.resize ( pocno_list.size() );

  generate ( pocno_remapping.begin(), pocno_remapping.end(), num_generator_t<uint>() );

  compare_list_items<vector<uint>, uint > comp ( pocno_list );


  for ( uint i = 0 ; i < alpha_ranges.size();i++ )
  {
    stable_sort ( pocno_remapping.begin() + alpha_ranges[i].first,
                  pocno_remapping.begin() + alpha_ranges[i].second,
                  comp
                );
  }

  vector < pair<uint, uint> >  pocno_ranges;

  for ( uint i = 0 ; i < alpha_ranges.size();i++ )
  {
    uint pocno_range_start = alpha_ranges[i].first;

    uint alpha_range_start = pocno_ranges.size();

    for ( uint j = alpha_ranges[i].first + 1 ; j < alpha_ranges[i].second;j++ )
    {
      if ( pocno_list[pocno_remapping[j]] != pocno_list[pocno_remapping[j-1]] )
      {
        pocno_ranges.push_back ( make_pair ( pocno_range_start, j ) );
        pocno_range_start = j;
      }
    }

    pocno_ranges.push_back ( make_pair ( pocno_range_start, alpha_ranges[i].second ) );

    alpha_ranges[i] = make_pair ( alpha_range_start, pocno_ranges.size() );
  }

  m_pockets.alpha_pocket_ranges     = new alpha_pocket_range_t[alpha_ranges.size() ];

  m_pockets.num_alpha_pocket_ranges = alpha_ranges.size() ;

  m_pockets.pocket_tet_ranges       = new pocket_tet_range_t[pocno_ranges.size() ];

  m_pockets.num_pocket_tet_ranges   = pocno_ranges.size();

  m_pockets.tet_idx                 = new uint [pocno_list.size() *4];

  m_pockets.num_tet_idx             = pocno_list.size() * 4 ;

  for ( uint i = 0 ; i < alpha_ranges.size();i++ )
  {
    m_pockets.alpha_pocket_ranges[i].start = alpha_ranges[i].first;
    m_pockets.alpha_pocket_ranges[i].end   = alpha_ranges[i].second;
  }

  for ( uint i = 0 ; i < pocno_ranges.size();i++ )
  {
    m_pockets.pocket_tet_ranges[i].start = pocno_ranges[i].first * 4;
    m_pockets.pocket_tet_ranges[i].end   = pocno_ranges[i].second * 4;
  }


  for ( uint i = 0 ; i < pocno_remapping.size();i++ )
  {
    uint v1 = tetra_list[tetno_list[pocno_remapping[i]] * 4 + 0];
    uint v2 = tetra_list[tetno_list[pocno_remapping[i]] * 4 + 1];
    uint v3 = tetra_list[tetno_list[pocno_remapping[i]] * 4 + 2];
    uint v4 = tetra_list[tetno_list[pocno_remapping[i]] * 4 + 3];

    m_pockets.tet_idx[4*i + 0 ] = v1;
    m_pockets.tet_idx[4*i + 1 ] = v2;
    m_pockets.tet_idx[4*i + 2 ] = v3;
    m_pockets.tet_idx[4*i + 3 ] = v4;

  }
}

void pocket_model_t::pockets_t::init()
{
  alpha_pocket_ranges     = NULL;
  pocket_tet_ranges       = NULL;
  tet_idx                 = NULL;

  num_alpha_pocket_ranges = 0;
  num_pocket_tet_ranges   = 0;
  num_tet_idx             = 0;

}

void pocket_model_t::pockets_t::destroy()
{
  if ( num_alpha_pocket_ranges != 0 )
    delete []alpha_pocket_ranges;

  if ( num_pocket_tet_ranges != 0 )
    delete []pocket_tet_ranges;

  if ( num_tet_idx != 0 )
    delete []tet_idx;

}

void pocket_model_t::clear_render()
{
  if ( m_ren != NULL )
  {
    delete m_ren;
    m_ren = NULL;
  }

  m_atom_indxs.clear();

  m_atom_bonds.clear();

}


pocket_model_t::pocket_model_t
( const char *tepoc_fn,  const char *tet_fn,  protein_rd_t * protein_rd )
    : m_protein_rd ( protein_rd ), m_ren ( NULL )
{
  m_pockets.init();

  read_file ( tepoc_fn, tet_fn );
}

pocket_model_t::~pocket_model_t ()
{
  clear_render();

  m_pockets.destroy();
}

bool pocket_model_t::check_alpha_num ( const uint &alphanum )
{
  if ( alphanum >= m_pockets.num_alpha_pocket_ranges )
  {
    _LOG ( "invalid alpha num" );
    _LOG_VAR ( alphanum );
    _LOG_VAR ( m_pockets.num_alpha_pocket_ranges );
    return false;
  }

  return true;
}

void pocket_model_t::setup_render ( const uint &alphanum , const uint &pocno )
{
  clear_render();

  if ( !check_alpha_num ( alphanum ) )
    return;

  uint pockets_begin = m_pockets.alpha_pocket_ranges[alphanum].start;

  uint pockets_end   = m_pockets.alpha_pocket_ranges[alphanum].end;

  uint tet_begin     = m_pockets.pocket_tet_ranges[pockets_begin].start;

  uint tet_end       = m_pockets.pocket_tet_ranges[pockets_end-1].end;

  if ( pocno != ( uint ) - 1 )
  {
    if ( pockets_begin + pocno >= pockets_end )
    {
      _LOG ( "invalid poc no" );
      _LOG_VAR ( alphanum );
      _LOG_VAR ( pocno );
      _LOG_VAR ( pockets_begin );
      _LOG_VAR ( pockets_end );
      return;
    }

    tet_begin     = m_pockets.pocket_tet_ranges[pockets_begin+pocno].start;

    tet_end       = m_pockets.pocket_tet_ranges[pockets_begin+pocno].end;
  }


  uint num_tet       = ( tet_end - tet_begin ) / 4;

  glutils::buf_obj_t ver_bo = m_protein_rd->get_coord_bo();

  glutils::buf_obj_t tet_bo
  ( m_pockets.tet_idx + tet_begin,
    GL_UNSIGNED_INT,
    4,
    GL_ELEMENT_ARRAY_BUFFER,
    num_tet*4*sizeof ( uint ),
    0 );

  glutils::buf_obj_t col_bo;

  m_ren = glutils::create_buffered_flat_tetrahedrons_ren ( ver_bo, tet_bo, col_bo );

  set<uint> atomset_set;

  for ( uint i = tet_begin;i < tet_end;i += 4 )
  {
    atomset_set.insert ( m_pockets.tet_idx[i+0] );
    atomset_set.insert ( m_pockets.tet_idx[i+1] );
    atomset_set.insert ( m_pockets.tet_idx[i+2] );
    atomset_set.insert ( m_pockets.tet_idx[i+3] );

    if ( ( m_pockets.tet_idx[i+0] >= m_protein_rd->get_protein()->get_num_atoms() ) ||
         ( m_pockets.tet_idx[i+1] >= m_protein_rd->get_protein()->get_num_atoms() ) ||
         ( m_pockets.tet_idx[i+2] >= m_protein_rd->get_protein()->get_num_atoms() ) ||
         ( m_pockets.tet_idx[i+3] >= m_protein_rd->get_protein()->get_num_atoms() ) )
    {
      _LOG ( "Screw up" );
    }
  }

  uint num_atomset = atomset_set.size();

  uint *atomset = new uint[num_atomset];

  copy ( atomset_set.begin(), atomset_set.end(), atomset );

  uint *bondset;
  uint num_bondset;

  m_protein_rd->get_protein()->get_subset_atom_bonds ( atomset, num_atomset, bondset, num_bondset );

  m_atom_indxs.src_ptr    = atomset;
  m_atom_indxs.src_type   = GL_UNSIGNED_INT ;
  m_atom_indxs.src_comp   = 1 ;
  m_atom_indxs.target     = GL_ELEMENT_ARRAY_BUFFER;
  m_atom_indxs.stride     = 0;
  m_atom_indxs.size       = num_atomset * sizeof ( uint );
  m_atom_indxs.upload();

  m_atom_bonds.src_ptr    = bondset;
  m_atom_bonds.src_type   = GL_UNSIGNED_INT ;
  m_atom_bonds.src_comp   = 2 ;
  m_atom_bonds.target     = GL_ELEMENT_ARRAY_BUFFER;
  m_atom_bonds.stride     = 0;
  m_atom_bonds.size       = num_bondset * 2 * sizeof ( uint );
  m_atom_bonds.upload();

  delete []bondset;
  delete []atomset;

}

int pocket_model_t::render() const
{
  if ( m_ren != NULL )
  {

    glColor3f ( 0.2, 0.8, 0.5 );

    glPushAttrib ( GL_ENABLE_BIT );

    glDisable ( GL_CULL_FACE );

    glEnable ( GL_RESCALE_NORMAL );

    m_ren->render();

    glPopAttrib ();
  }

  return 0;
}


