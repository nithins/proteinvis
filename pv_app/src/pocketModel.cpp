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



void pocket_model_t::read_file ( const std::string &  tepoc_filename, const std::string & tet_filename )
{

  using boost::bind;

  fstream tepocFile ( tepoc_filename.c_str(), ios::in );

  vector<uint> alpha_list;

  vector<uint> tetno_list;

  vector<uint> pocno_list;

  if(tepocFile.is_open() == false)
    throw std::runtime_error("pocket file does not exist or is empty");

  while ( !tepocFile.eof() )
  {
    string line;

    getline ( tepocFile, line );

    vector<string> tokens;

    line = stripWS ( stripLineComments ( line ) );

    if ( line.size() == 0 )
      continue;

    tokenize_string ( line, tokens );

    if(tokens.size() < 3) continue;

    alpha_list.push_back ( atoi ( tokens[0].c_str() ) );

    tetno_list.push_back ( atoi ( tokens[1].c_str() ) - 1 );

    pocno_list.push_back ( atoi ( tokens[2].c_str() ) );
  }

  fstream tetFile ( tet_filename.c_str(), ios::in );

  if(tetFile.is_open() == false)
    throw std::runtime_error("tet file does not exist or is empty");

  string num_tets_line;

  getline ( tetFile, num_tets_line );

  glutils::quad_idx_list_t tetra_list;

  while ( !tetFile.eof() )
  {
    string line;

    getline ( tetFile, line );

    line = stripWS ( stripLineComments ( line ) );

    if ( line.size() == 0 )
      continue;

    vector<string> tokens;

    tokenize_string ( line, tokens );

    glutils::quad_idx_t tetra
        (atoi ( tokens[1].c_str() ) - 1 ,atoi ( tokens[2].c_str() ) - 1,
         atoi ( tokens[3].c_str() ) - 1 ,atoi ( tokens[4].c_str() ) - 1 );

    tetra_list.push_back (tetra);
  }

  std::vector<alpha_pocket_range_t>  alpha_ranges;

  uint range_start = 0;

  for ( uint i = 1 ; i < alpha_list.size();i++ )
  {
    if ( alpha_list[i] != alpha_list[i-1] )
    {
      alpha_ranges.push_back ( alpha_pocket_range_t ( range_start, i ) );
      range_start = i;
    }
  }

  alpha_ranges.push_back ( alpha_pocket_range_t ( range_start, alpha_list.size() ) );

  vector<uint> pocno_remapping;

  pocno_remapping.resize ( pocno_list.size() );

  generate ( pocno_remapping.begin(), pocno_remapping.end(), num_generator_t<uint>() );

  compare_list_items<vector<uint>, uint > comp ( pocno_list );


  for ( uint i = 0 ; i < alpha_ranges.size();i++ )
  {
    stable_sort ( pocno_remapping.begin() + alpha_ranges[i][0],
                  pocno_remapping.begin() + alpha_ranges[i][1],
                  comp
                );
  }

  std::vector<pocket_tet_range_t>  pocno_ranges;

  for ( uint i = 0 ; i < alpha_ranges.size();i++ )
  {
    uint pocno_range_start = alpha_ranges[i][0];

    uint alpha_range_start = pocno_ranges.size();

    for ( uint j = alpha_ranges[i][0] + 1 ; j < alpha_ranges[i][1];j++ )
    {
      if ( pocno_list[pocno_remapping[j]] != pocno_list[pocno_remapping[j-1]] )
      {
        pocno_ranges.push_back ( pocket_tet_range_t ( pocno_range_start, j ) );
        pocno_range_start = j;
      }
    }

    pocno_ranges.push_back ( pocket_tet_range_t ( pocno_range_start, alpha_ranges[i][1]) );

    alpha_ranges[i] = alpha_pocket_range_t ( alpha_range_start, pocno_ranges.size() );
  }

  m_pockets.alpha_pocket_ranges.resize(alpha_ranges.size() );

  m_pockets.pocket_tet_ranges.resize(pocno_ranges.size() );

  m_pockets.tet_idx.resize(pocno_list.size());

  std::copy(alpha_ranges.begin(),alpha_ranges.end(),
            m_pockets.alpha_pocket_ranges.begin());

  std::copy(pocno_ranges.begin(),pocno_ranges.end(),
            m_pockets.pocket_tet_ranges.begin());

  for ( uint i = 0 ; i < pocno_remapping.size();i++ )
    m_pockets.tet_idx[i] = tetra_list[tetno_list[pocno_remapping[i]]];
}


void pocket_model_t::pockets_t::destroy()
{
}

void pocket_model_t::clear_render()
{
}


pocket_model_t::pocket_model_t
    ( const std::string & tepoc_fn,
      const std::string &tet_fn,
      boost::shared_ptr<protein_rd_t> protein_rd)
    : m_protein_rd ( protein_rd )
{
  read_file ( tepoc_fn, tet_fn );
}

pocket_model_t::~pocket_model_t ()
{
}

bool pocket_model_t::check_alpha_num ( const uint &alphanum )
{
  if ( alphanum >= m_pockets.alpha_pocket_ranges.size())
    throw std::logic_error("invalid alpha num");

  return true;
}

void pocket_model_t::setup_render ( const uint &alphanum , const uint &pocno )
{
  clear_render();

  uint pockets_begin = m_pockets.alpha_pocket_ranges[alphanum][0];

  uint pockets_end   = m_pockets.alpha_pocket_ranges[alphanum][1];

  uint tet_begin     = m_pockets.pocket_tet_ranges[pockets_begin][0];

  uint tet_end       = m_pockets.pocket_tet_ranges[pockets_end-1][1];

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

    tet_begin     = m_pockets.pocket_tet_ranges[pockets_begin+pocno][0];

    tet_end       = m_pockets.pocket_tet_ranges[pockets_begin+pocno][1];
  }


  glutils::bufobj_ptr_t ver_bo = m_protein_rd->get_coord_bo();

  glutils::quad_idx_list_t tets(tet_end - tet_begin );

  std::copy(m_pockets.tet_idx.begin()+tet_begin,m_pockets.tet_idx.begin()+tet_end,tets.begin());

  glutils::bufobj_ptr_t tet_bo = glutils::make_buf_obj(tets);

  m_ren.reset(glutils::create_buffered_flat_tetrahedrons_ren ( ver_bo, tet_bo));

  set<uint> atomset_set;

  for ( uint i = tet_begin;i < tet_end;i++ )
  {
    atomset_set.insert ( m_pockets.tet_idx[i][0] );
    atomset_set.insert ( m_pockets.tet_idx[i][1] );
    atomset_set.insert ( m_pockets.tet_idx[i][2] );
    atomset_set.insert ( m_pockets.tet_idx[i][3] );

    if ( ( m_pockets.tet_idx[i][0] >= m_protein_rd->get_protein()->get_num_atoms() ) ||
         ( m_pockets.tet_idx[i][1] >= m_protein_rd->get_protein()->get_num_atoms() ) ||
         ( m_pockets.tet_idx[i][2] >= m_protein_rd->get_protein()->get_num_atoms() ) ||
         ( m_pockets.tet_idx[i][3] >= m_protein_rd->get_protein()->get_num_atoms() ) )
    {
      _LOG ( "Screw up" );
    }
  }

  glutils::point_idx_list_t atom_idx_list(atomset_set.size());

  copy ( atomset_set.begin(), atomset_set.end(), atom_idx_list.begin());

  uint * bondlist_c;
  uint   num_bonds;

  m_protein_rd->get_protein()->get_subset_atom_bonds
      ( atom_idx_list.data(),atom_idx_list.size(), bondlist_c,num_bonds );

  glutils::line_idx_list_t bondlist(num_bonds);

  std::copy(bondlist_c,  bondlist_c + 2*num_bonds, (uint*)bondlist.data());

  delete []bondlist_c;

  m_atom_indxs = glutils::make_buf_obj(atom_idx_list);

  m_atom_bonds = glutils::make_buf_obj(bondlist);


}

int pocket_model_t::render() const
{
  if ( m_ren != NULL )
  {

    glColor3f ( 0.2, 0.8, 0.5 );

    glPushAttrib ( GL_ENABLE_BIT|GL_POLYGON_BIT );

    glDisable ( GL_CULL_FACE );

    glEnable ( GL_RESCALE_NORMAL );

    glPolygonMode ( GL_FRONT_AND_BACK, GL_FILL );

    m_ren->render();

    glPopAttrib ();
  }

  return 0;
}
