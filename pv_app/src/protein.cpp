
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
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>
#include <list>
#include <set>
#include <cmath>
#include <climits>
#include <cstring>

#include <boost/regex.hpp>

#include <logutil.h>
#include <cpputils.h>

#include <protein.h>

#include <limits>

uchar g_standard_colors[][3] =
{
  {0x00, 0x00, 0x00}, //Black
  {0x80, 0x00, 0x00}, //Maroon
  {0x00, 0x80, 0x00}, //Green
  {0x00, 0x00, 0x80}, //Navy
  {0xC0, 0xC0, 0xC0}, //Silver
  {0xFF, 0x00, 0x00}, //Red
  {0x00, 0xFF, 0x00}, //Lime
  {0x00, 0x00, 0xFF}, //Blue
  {0x80, 0x80, 0x80}, //Gray
  {0x80, 0x00, 0x80}, //Purple
  {0x80, 0x80, 0x00}, //Olive
  {0x00, 0x80, 0x80}, //Teal
  {0xFF, 0xFF, 0xFF}, //White
  {0xFF, 0x00, 0xFF}, //Fuchsia
  {0xFF, 0xFF, 0x00}, //Yellow
  {0x00, 0xFF, 0xFF}, //Aqua
  {0xFF, 0x7F, 0x00}, //Orange
};

#define RGB_BLACK    g_standard_colors[0]
#define RGB_MAROON   g_standard_colors[1]
#define RGB_GREEN    g_standard_colors[2]
#define RGB_NAVY     g_standard_colors[3]
#define RGB_SILVER   g_standard_colors[4]
#define RGB_RED      g_standard_colors[5]
#define RGB_LIME     g_standard_colors[6]
#define RGB_BLUE     g_standard_colors[7]
#define RGB_GRAY     g_standard_colors[8]
#define RGB_PURPL    g_standard_colors[9]
#define RGB_OLIVE    g_standard_colors[10]
#define RGB_TEAL     g_standard_colors[11]
#define RGB_WHITE    g_standard_colors[12]
#define RGB_FUCHSIA  g_standard_colors[13]
#define RGB_YELLOW   g_standard_colors[14]
#define RGB_AQUA     g_standard_colors[15]
#define RGB_ORANGE   g_standard_colors[16]

std_atom_types_t protein_t::std_atom_types[] =
{
  {"Uk", "Unknown", RGB_AQUA, 2.0, 2.0},
  {"H", "Hydrogen", RGB_WHITE, 0.320, 1.100},
  {"O", "Oxygen", RGB_RED, 0.680, 1.348},
  {"Cl", "Chlorine", RGB_GREEN, 2.0, 2.0},
  {"N", "Nitrogen", RGB_BLUE, 0.680, 1.400},
  {"C", "Carbon", RGB_GRAY, 0.720, 1.548},
  {"S", "Sulpher", RGB_YELLOW, 1.020, 1.808},
  {"P", "Phosphorus", RGB_ORANGE, 1.9, 1.9}
};

const double g_water_moleclue_ball_radius = 1.4;

using namespace std;

bool read_crd_file ( const char *, protein_t & );
bool read_pdb_file ( const char *, protein_t & );

uint get_std_atom_type_idx ( string name )
{
  uint num_std_atom_types = sizeof ( protein_t::std_atom_types ) / sizeof ( std_atom_types_t );

  if ( num_std_atom_types == 0 )
    return 0;

  static string atomnames_re_str;

  if ( atomnames_re_str.size() == 0 )
  {

    stringstream atomnames_re_sstream;

    atomnames_re_sstream << "^(" << protein_t::std_atom_types[0].name;

    for ( uint i = 1 ; i < num_std_atom_types ;i++ )
      atomnames_re_sstream << "|" << protein_t::std_atom_types[i].name;

    atomnames_re_sstream << ")";

    atomnames_re_str = atomnames_re_sstream.str();
  }

  static map<string, uint> name_idx_map;

  if ( name_idx_map.size() == 0 )
  {
    for ( uint i = 0 ; i < num_std_atom_types;i++ )
      name_idx_map[protein_t::std_atom_types[i].name] = i;
  }

  boost::regex atomnames_re ( atomnames_re_str );

  boost::smatch atomnames_matches;

  if ( boost::regex_search ( name, atomnames_matches, atomnames_re ) )
    return name_idx_map[string ( atomnames_matches[0].first,
                                 atomnames_matches[0].second ) ];
  else
    return 0;
}

string chompLine ( const string& line, const char &comment_char = '#' )
{
  return stripWS ( stripLineComments ( line, comment_char ) );
}

protein_t::protein_t ( const char *filename )
{
  init();

  string filename_str ( filename );

  string file_extn ( filename_str.begin() + ( filename_str.size() - 4 ), filename_str.end() );

  if ( file_extn == ".crd" )
    read_crd_file ( filename, *this );
  else
  {
    if ( file_extn == ".pdb" )
      read_pdb_file ( filename, *this );
    else
    {
      _LOG ( "cannot read extn type " << file_extn );
    }
  }

  collect_bb_atoms();

  collect_bb_bonds();
}

protein_t::~protein_t()
{
  destroy();
}

void protein_t::destroy()
{
  if ( num_atoms != 0 )
    delete []atoms;

  if ( num_bonds != 0 )
    delete []bonds;

  if ( num_acids != 0 )
    delete []acids;

  if ( num_chains != 0 )
    delete []chains;

  if ( num_atom_types != 0 )
    delete []atom_types;

  if ( num_acid_types != 0 )
    delete []acid_types;

  if ( num_helices != 0 )
    delete []helices;

  if ( num_sheets!= 0 )
    delete []sheets;

  init();
}

void protein_t::init()
{
  atoms          = NULL;
  bonds          = NULL;
  acids          = NULL;
  chains         = NULL;
  atom_types     = NULL;
  acid_types     = NULL;
  helices        = NULL;
  sheets         = NULL;

  num_atoms      = 0;
  num_bonds      = 0;
  num_acids      = 0;
  num_chains     = 0;
  num_atom_types = 0;
  num_acid_types = 0;
  num_helices    = 0;
  num_sheets     = 0;
}

void protein_t::print_atom_types ( ostream &o )
{
  o << "================================" << endl;
  o << "     Atom type information      " << endl;
  o << "--------------------------------" << endl;
  o << "Name \t\t Radius \t\t Std Name  " << endl;
  o << "--------------------------------" << endl;

  for ( uint i = 0 ; i < num_atom_types;i++ )
  {
    o
        << atom_types[i].name << "\t\t"
        << atom_types[i].radius << "\t\t"
        << std_atom_types[atom_types[i].std_type_idx].fullname << endl;
  }

  o << "================================" << endl;
}

void protein_t::print_acid_types ( ostream &o )
{
  o << "=====================" << endl;
  o << "Acid type information" << endl;
  o << "---------------------" << endl;
  o << "Name \t\t Inst Idx   " << endl;
  o << "---------------------" << endl;

  for ( uint i = 0 ; i < num_acid_types;i++ )
  {
    o << acid_types[i].name << "\t\t" << acid_types[i].inst_idx << endl;
  }

  o << "=====================" << endl;
}

void protein_t::print_atoms ( std::ostream &o )
{
  o << "===========================================" << endl;
  o << "              Atom Information             " << endl;
  o << "-------------------------------------------" << endl;
  o << "Name \t\tx\t\ty\t\tz\t\tbonded atoms       " << endl;
  o << "-------------------------------------------" << endl;

  for ( uint i = 0 ; i < num_atoms;i++ )
  {
    o
        << atom_types[atoms[i].type_idx].name << "\t\t"
        << atoms[i].x << "\t\t"
        << atoms[i].y << "\t\t"
        << atoms[i].z << "\t\t";

    for ( uint j = atoms[i].bond_start;j < atoms[i].bond_end;j++ )
    {
      o << bonds[j+1] << "\t\t";
    }

    o << endl;
  }

  o << "===========================================" << endl;

}

void protein_t::print_acids ( std::ostream &o )
{
  o << "===========================================" << endl;
  o << "              Acid Information             " << endl;
  o << "-------------------------------------------" << endl;
  o << "Name \t start \t end                       " << endl;
  o << "-------------------------------------------" << endl;

  for ( uint i = 0 ; i < num_acids;i++ )
  {
    o
        << acid_types[acids[i].type_idx].name << "\t"
        << acids[i].start << "\t"
        << acids[i].end << endl;
  }

  o << "===========================================" << endl;

}

void protein_t::print_chains ( std::ostream &o )
{
  o << "===========================================" << endl;
  o << "              Chain Information            " << endl;
  o << "-------------------------------------------" << endl;
  o << "Name \t start \t end                       " << endl;
  o << "-------------------------------------------" << endl;

  for ( uint i = 0 ; i < num_chains;i++ )
  {
    o
        << chains[i].name << "\t"
        << chains[i].start << "\t"
        << chains[i].end << endl;
  }

  o << "===========================================" << endl;
}

void protein_t::print_stats ( std::ostream &o )
{
  o << "===========================================" << endl;
  o << "              Protein Stats                " << endl;
  o << "-------------------------------------------" << endl;
  o << "Num Atoms       : " << num_atoms  << endl;
  o << "Num Acids       : " << num_acids  << endl;
  o << "Num Chains      : " << num_chains << endl;
  o << "Num Bonds       : " << num_bonds << endl;
  o << "Num Atom Types  : " << num_atom_types << endl;
  o << "Num Acid Types  : " << num_acid_types << endl;
  o << "===========================================" << endl;
}

void protein_t::print_incorrect_bonds ( std::ostream &o )
{
  const double bond_len_treshhold = 2.0;

  double *bond_lengths = new double[num_bonds];

  for ( uint i = 0; i < num_bonds;i++ )
  {
    uint a1 = bonds[2*i+0];
    uint a2 = bonds[2*i+1];

    bond_lengths[i] = sqrt ( pow ( ( atoms[a1].x - atoms[a2].x ), 2 ) +
                             pow ( ( atoms[a1].y - atoms[a2].y ), 2 ) +
                             pow ( ( atoms[a1].z - atoms[a2].z ), 2 ) );
  }

  o << "===========================================" << endl;

  o << "  Bonds exceeding " << bond_len_treshhold << " Angstroms " << endl;
  o << "-------------------------------------------" << endl;
  o << "Len\tind1,ind2\tname1,name2" << endl;
  o << "-------------------------------------------" << endl;


  for ( uint i = 0 ;i < num_bonds;++i )
  {
    if ( bond_lengths[i] >= bond_len_treshhold )
    {
      uint a1  = bonds[2*i+0];
      uint a2  = bonds[2*i+1];

      string &name1 = atom_types[atoms[a1].type_idx].name;
      string &name2 = atom_types[atoms[a2].type_idx].name;

      o << bond_lengths[i] << "\t" << a1 + 1 << "," << a2 + 1 << "\t\t" << name1 << "," << name2 << endl;

    }
  }

  o << "===========================================" << endl;

  delete [] bond_lengths;
}

void protein_t::print_info ( std::ostream &o )
{
  print_atom_types ( o );
  print_acid_types ( o );
  print_chains ( o );
  print_stats ( o );
  //print_acids ( o );
}

bool comp_atomno_acid_end ( uint atomno,acid_t a)
{
  if (  atomno <a.end)
    return true;
  else
    return false;
}

uint protein_t::get_atom_acid_idx(uint atomno)
{
  acid_t *acid = upper_bound(acids,acids+num_acids,atomno,
                             comp_atomno_acid_end);

  return (acid -acids);
}

bool comp_acidno_chain_end ( uint acidno,chain_t c)
{
  if (  acidno <c.end)
    return true;
  else
    return false;
}

uint protein_t::get_acid_chain_idx(uint acidno)
{
  chain_t *chain= upper_bound(chains,chains+num_chains,acidno,
                              comp_acidno_chain_end);

  return chain-chains;
}


bool protein_t::is_backbone_atom ( uint atomno )
{
  string& atom_name = atom_types[atoms[atomno].type_idx].name;

  if ( atom_name == "N" ||
       atom_name == "CA" ||
       atom_name == "C" ||
       atom_name == "O" ||
       atom_name == "OT1" ) // in crd files this is used
    return true;

  return false;
}

bool protein_t::is_ca_atom(uint atomno) const
{
    const string& atom_name = atom_types[atoms[atomno].type_idx].name;

    if (atom_name == "CA" ) // in crd and pdb files, this is used
      return true;

    return false;
}

bool protein_t::is_o_atom(uint atomno) const
{
    const string& atom_name = atom_types[atoms[atomno].type_idx].name;

    if (atom_name == "O" ||
        atom_name == "OT1") // in crd and pdb files, this is used
      return true;

    return false;
}

void protein_t::check_only_four_bb_atoms_per_acid()
{
  for ( uint i = 0 ; i < num_acids;i++ )
  {
    uint num_bb_atoms_in_acid = 0;

    for ( uint j = acids[i].start ; j < acids[i].end;j++ )
    {
      if ( is_backbone_atom ( j ) )
        num_bb_atoms_in_acid++;
    }

    if ( ( num_bb_atoms_in_acid != 4 ) &&
         ( acid_types[acids[i].type_idx].name != "NTER" ) &&
         ( acid_types[acids[i].type_idx].name != "CTER" )
      )
      {
      _LOG ( "Incorrect num of bb atoms in acid" );
      _LOG ( "Acid No       :" << i );
      _LOG ( "Acid Name     :" << acid_types[acids[i].type_idx].name );
      _LOG ( "Num BB atoms  :" << num_bb_atoms_in_acid );
      _LOG ( "Atom start    :" << acids[i].start );
      _LOG ( "Atom end      :" << acids[i].end );
    }
  }
}

void protein_t::collect_bb_atoms()
{
  vector<uint> bb_atoms_idx_vec;

  for ( uint i = 0 ;i < num_atoms;i++ )
  {
    if ( is_backbone_atom ( i ) )
      bb_atoms_idx_vec.push_back ( i );
  }

  num_bb_atoms = bb_atoms_idx_vec.size();
  bb_atoms_idx = new uint[num_bb_atoms];

  copy ( bb_atoms_idx_vec.begin(), bb_atoms_idx_vec.end(), bb_atoms_idx );
}

void protein_t::collect_bb_bonds()
{
  vector<uint> bb_bonds_vec;

  for ( uint i = 0 ;i < num_atoms;i++ )
  {
    if ( is_backbone_atom ( i ) )
    {
      for ( uint j = atoms[i].bond_start;j < atoms[i].bond_end; j += 2 )
      {
        if ( is_backbone_atom ( bonds[j+1] ) )
        {
          bb_bonds_vec.push_back ( i );
          bb_bonds_vec.push_back ( bonds[j+1] );
        }
      }
    }
  }

  num_bb_bonds = bb_bonds_vec.size() / 2;

  bb_bonds     = new uint[num_bb_bonds*2];

  copy ( bb_bonds_vec.begin(), bb_bonds_vec.end(), bb_bonds );
}

void protein_t::get_subset_atom_bonds ( const uint *atom_indxs,
                                        const uint &num_atom_indxs,
                                        uint *&subset_atom_bonds,
                                        uint &num_subset_atom_bonds ) const
{
  set<uint> atomset ( atom_indxs, atom_indxs + num_atom_indxs );

  vector<uint> atomset_bonds;

  for ( set<uint>::iterator it = atomset.begin();it != atomset.end(); ++it )
  {
    uint atomno = *it;

    for ( uint i = atoms[atomno].bond_start ; i < atoms[atomno].bond_end;i += 2 )
    {
      if ( atomset.find ( bonds[i+1] ) != atomset.end() )
      {
        atomset_bonds.push_back ( atomno );
        atomset_bonds.push_back ( bonds[i+1] );
      }
    }
  }

  num_subset_atom_bonds = atomset_bonds.size() / 2;

  subset_atom_bonds = new uint[num_subset_atom_bonds*2];

  copy ( atomset_bonds.begin(), atomset_bonds.end(), subset_atom_bonds );
}


void protein_rd_t::upload_data_items()
{

  // a bo to hold the coord of each atom

  double * atom_coords = new double[m_protein->get_num_atoms() *3];

  for ( uint i = 0 ;i < m_protein->get_num_atoms();i++ )
  {
    atom_coords[3*i+0] = m_protein->get_atoms() [i].x;
    atom_coords[3*i+1] = m_protein->get_atoms() [i].y;
    atom_coords[3*i+2] = m_protein->get_atoms() [i].z;
  }

  m_atom_coord_bo = glutils::buf_obj_t::create_bo
                    (atom_coords,GL_DOUBLE,3,GL_ARRAY_BUFFER,
                     sizeof ( GLdouble ) * 3 * m_protein->get_num_atoms(),0);

  delete []atom_coords;

  // a bo to hold the radii of each atom

  double * atom_radii = new double[m_protein->get_num_atoms() ];

  for ( uint i = 0 ;i < m_protein->get_num_atoms();i++ )
  {
    atom_radii[i]          = m_protein->get_atoms() [i].radius;
  }


  m_atom_radii_bo = glutils::buf_obj_t::create_bo
                    (atom_radii,GL_DOUBLE,1,GL_ARRAY_BUFFER,
                     sizeof ( GLdouble ) * m_protein->get_num_atoms(),0);

  delete []atom_radii;

  // a bo to hold the bonds of each atom
  m_atom_bonds_bo = glutils::buf_obj_t::create_bo
                    (m_protein->get_bonds(),GL_UNSIGNED_INT,2,
                     GL_ELEMENT_ARRAY_BUFFER,
                     sizeof ( GLuint ) * 2 * m_protein->get_num_bonds (),0);

  // vbo for holding the indices of all atoms that are on the backbone

  m_bb_atom_indices_bo = glutils::buf_obj_t::create_bo
                         (m_protein->get_bb_atoms_idx(),GL_UNSIGNED_INT,1,
                          GL_ELEMENT_ARRAY_BUFFER,
                          sizeof ( GLuint ) * m_protein->get_num_bb_atoms(),0);

  // vbo for holding the indices of all atoms that are on the backbone

  m_bb_bond_indices_bo = glutils::buf_obj_t::create_bo
                         (m_protein->get_bb_bonds(),GL_UNSIGNED_INT,2,
                          GL_ELEMENT_ARRAY_BUFFER,
                          sizeof ( GLuint )*2*m_protein->get_num_bb_bonds(),0);

  // a bo to hold the set of c-alpha atoms

  glutils::point_idx_list_t ca_idxs;
  for ( uint i = 0 ;i < m_protein->get_num_atoms();i++ )
    if(m_protein->is_ca_atom(i))
      ca_idxs.push_back(i);

  m_ca_atom_indices_bo = glutils::make_buf_obj(ca_idxs);
}

void protein_rd_t::compute_extent()
{
  double extent[] =
  {
    std::numeric_limits<double>::max(),
    std::numeric_limits<double>::min(),
    std::numeric_limits<double>::max(),
    std::numeric_limits<double>::min(),
    std::numeric_limits<double>::max(),
    std::numeric_limits<double>::min(),
  };

  for ( uint i = 0 ;i < m_protein->get_num_atoms();i++ )
  {

    double radius = m_protein->get_atoms() [i].radius;

    extent[0] = min ( m_protein->get_atoms() [i].x - radius , extent[0] );
    extent[1] = max ( m_protein->get_atoms() [i].x + radius , extent[1] );
    extent[2] = min ( m_protein->get_atoms() [i].y - radius , extent[2] );
    extent[3] = max ( m_protein->get_atoms() [i].y + radius , extent[3] );
    extent[4] = min ( m_protein->get_atoms() [i].z - radius , extent[4] );
    extent[5] = max ( m_protein->get_atoms() [i].z + radius , extent[5] );
  }

  copy ( extent, extent + 6, m_extent );
}

protein_rd_t::protein_rd_t ( boost::shared_ptr<protein_t> p ) : m_protein ( p )
{
  upload_data_items();

  compute_extent();
}

void protein_rd_t::get_extent ( double *extent )
{
  copy ( m_extent, m_extent + 6, extent );
}

protein_rd_t::~protein_rd_t ( )
{
}

glutils::bufobj_ptr_t protein_rd_t::get_coord_bo()
{
  return m_atom_coord_bo;
}

glutils::bufobj_ptr_t protein_rd_t::get_bonds_bo()
{
  return m_atom_bonds_bo;
}

glutils::bufobj_ptr_t protein_rd_t::get_radii_bo()
{
  return m_atom_radii_bo;
}

glutils::bufobj_ptr_t protein_rd_t::get_bb_bonds_bo()
{
  return m_bb_bond_indices_bo;
}

glutils::bufobj_ptr_t protein_rd_t::get_bb_coord_bo()
{
  return m_bb_atom_indices_bo;
}

glutils::bufobj_ptr_t protein_rd_t::get_ca_atoms_bo()
{
    return m_ca_atom_indices_bo;
}

glutils::bufobj_ptr_t protein_rd_t::get_o_atoms_bo()
{
    return m_o_atom_indices_bo;
}

protein_grouping_t::protein_grouping_t(boost::shared_ptr<protein_t> p)
{
  m_protein       = p;
  m_grouping_type = GROUP_ATOMS_COUNT;
  m_num_groups    = 0;
  m_group_colors.clear();

  set_grouping_type((eGroupAtomsBy)0);

  update_atom_color_bo();
}

protein_grouping_t::~protein_grouping_t()
{
}

void protein_grouping_t::set_group_color(glutils::color_t col,int i)
{
  if(i<0 || i>= m_num_groups )
    throw string("unknown group")+__func__;

  m_group_colors[i] =col;
}

std::string protein_grouping_t::get_group_name(int groupno) const
{

  if(groupno <0 || groupno >= (int)m_num_groups)
  {
    throw string("unknown group")+__func__;
    return string();
  }

  stringstream ss;

  switch(m_grouping_type)
  {
  case GROUP_ATOMS_ATOM:
    {
      int atom_type = m_protein->get_atoms()[groupno].type_idx;
      ss<<groupno<<" :: "<<m_protein->get_atom_types()[atom_type].name;
      return ss.str();
    }
  case GROUP_ATOMS_ATOM_TYPE:
    {
      return m_protein->get_atom_types()[groupno].name;
    }

  case GROUP_ATOMS_ACID:
    {
      int acid_type = m_protein->get_acids()[groupno].type_idx;
      ss<<groupno<<" :: "<<m_protein->get_acid_types()[acid_type].name;
      return ss.str();
    }
  case GROUP_ATOMS_ACID_TYPE:
    {
      return m_protein->get_acid_types()[groupno].name;
    }

  case GROUP_ATOMS_CHAIN:
    {
      return m_protein->get_chains()[groupno].name;
    }

  case GROUP_ATOMS_ALL:
    {
      return "all";
    }
  default:
    throw string("unknown group")+__func__;
    return "all";
  }
}

glutils::color_t
    protein_grouping_t::get_group_color(int groupno) const
{
  if(groupno <0 || groupno >= (int)m_num_groups)
  {
    throw string("unknown group")+__func__;
  }

  return m_group_colors[groupno];
}

uint protein_grouping_t::get_updated_num_groups()const
{
  switch(m_grouping_type)
  {
  case GROUP_ATOMS_ATOM:return m_protein->get_num_atoms();
  case GROUP_ATOMS_ATOM_TYPE:return m_protein->get_num_atom_types();
  case GROUP_ATOMS_ACID:return m_protein->get_num_acids();
  case GROUP_ATOMS_ACID_TYPE:return m_protein->get_num_acid_types();
  case GROUP_ATOMS_CHAIN:return m_protein->get_num_chains();
  case GROUP_ATOMS_ALL: return 1;
  default:throw string("unknown group")+__func__;return 0;
  }
}

std::string protein_grouping_t::get_groupby_displaystr
    (protein_grouping_t::eGroupAtomsBy g)
{
  switch(g)
  {
  case GROUP_ATOMS_ATOM:return "By Atom";
  case GROUP_ATOMS_ATOM_TYPE:return "By Atom Type";
  case GROUP_ATOMS_ACID:return "By Acid";
  case GROUP_ATOMS_ACID_TYPE:return "By Acid type";
  case GROUP_ATOMS_CHAIN:return "By Chain";
  case GROUP_ATOMS_ALL: return "All";
  default:throw string("unknown group")+__func__;return 0;
  }
}

void protein_grouping_t::set_grouping_type(protein_grouping_t::eGroupAtomsBy g)
{
  if(g == m_grouping_type)
    return;

  if(g< 0 ||g >= GROUP_ATOMS_COUNT)
  {
    stringstream ss;
    ss<<"invalid group no    "<<endl;
    ss<<"GROUP_ATOMS_COUNT = "<<GROUP_ATOMS_COUNT<<endl;
    ss<<"g                 = "<<g<<endl;
    ss<<__func__<<endl;
    throw ss.str();
  }

  m_grouping_type = g;

  m_num_groups = get_updated_num_groups();

  set_default_color_values();
}

void protein_grouping_t::set_default_color_values()
{

  m_group_colors.resize(m_num_groups);

  double dgroup_ct      = m_num_groups;

  double dgroup_ct_sqrt = sqrt(dgroup_ct);

  int     group_ct_sqrt  = std::ceil(dgroup_ct_sqrt);

  for( int i = 0 ; i < m_num_groups;++i)
  {
    switch(m_grouping_type)
    {
    case GROUP_ATOMS_ATOM:
      {
        uchar col_uc[3];

        int atom_type = m_protein->get_atoms()[i].type_idx;
        int std_type  = m_protein->get_atom_types()[atom_type].std_type_idx;
        const uchar * std_atom_col = protein_t::std_atom_types[std_type].color;

        copy(std_atom_col,std_atom_col+3,col_uc);

        m_group_colors[i][0] = (double)col_uc[0]/(double)UCHAR_MAX;
        m_group_colors[i][1] = (double)col_uc[1]/(double)UCHAR_MAX;
        m_group_colors[i][2] = (double)col_uc[2]/(double)UCHAR_MAX;
        break;
      }
    case GROUP_ATOMS_ATOM_TYPE:
      {
        uchar col_uc[3];

        int std_type  = m_protein->get_atom_types()[i].std_type_idx;
        const uchar * std_atom_col = protein_t::std_atom_types[std_type].color;

        copy(std_atom_col,std_atom_col+3,col_uc);
        m_group_colors[i][0] = (double)col_uc[0]/(double)UCHAR_MAX;
        m_group_colors[i][1] = (double)col_uc[1]/(double)UCHAR_MAX;
        m_group_colors[i][2] = (double)col_uc[2]/(double)UCHAR_MAX;
        break;
      }

    case GROUP_ATOMS_ACID:
    case GROUP_ATOMS_ACID_TYPE:
    case GROUP_ATOMS_CHAIN:
    case GROUP_ATOMS_ALL:
    default:
      {
        m_group_colors[i][0] =
            0.5;
        m_group_colors[i][1] =
            0.25 + ((double)((i/group_ct_sqrt)%group_ct_sqrt))*0.50/group_ct_sqrt;
        m_group_colors[i][2] =
            0.25 + ((double)(i%group_ct_sqrt))*0.50/group_ct_sqrt;
        break;
      }
    }
  }
}

void  protein_grouping_t::update_atom_color_bo()
{
  double * atom_colors = new double [m_protein->get_num_atoms()*3];

  for(uint i = 0 ;i < m_protein->get_num_atoms();++i)
  {
    glutils::color_t col;
    switch(m_grouping_type)
    {
    case GROUP_ATOMS_ATOM:
      col = m_group_colors[i];
      break;
    case GROUP_ATOMS_ATOM_TYPE:
      col = m_group_colors[m_protein->get_atoms()[i].type_idx];
      break;
    case GROUP_ATOMS_ACID:
      col = m_group_colors[m_protein->get_atom_acid_idx(i)];
      break;
    case GROUP_ATOMS_ACID_TYPE:
      {
        uint acid_idx = m_protein->get_atom_acid_idx(i);
        col = m_group_colors[m_protein->get_acids()[acid_idx].type_idx];
        break;
      }

    case GROUP_ATOMS_CHAIN:
      {
        uint acid_idx  = m_protein->get_atom_acid_idx(i);
        col = m_group_colors[m_protein->get_acid_chain_idx(acid_idx)];
        break;
      }
    case GROUP_ATOMS_ALL:
      {
        col = m_group_colors[0];
        break;
      }
    default:
      {
        stringstream ss;
        ss<<"invalid group no    "<<endl;
        ss<<"GROUP_ATOMS_COUNT = "<<GROUP_ATOMS_COUNT<<endl;
        ss<<"g                 = "<<m_grouping_type<<endl;
        ss<<__func__<<endl;
        throw ss.str();
      }
      break;
    }

    atom_colors[3*i+0] = col[0];
    atom_colors[3*i+1] = col[1];
    atom_colors[3*i+2] = col[2];
  }

  m_atom_color_bo = glutils::buf_obj_t::create_bo
                    (atom_colors,GL_DOUBLE,3,GL_ARRAY_BUFFER,
                     3*sizeof(double)*m_protein->get_num_atoms(),0);

  delete []atom_colors;
}

void protein_t::transform_and_save_crd
    ( const char * in_filename,
      const char *out_filename,
      const double* trans_mat )
{
  string line;

  stringstream linestream;

  uint num_atoms = 0;

  if ( string ( in_filename ) == string ( out_filename ) )
  {
    _LOG ( "cannot write to same file" );
    _LOG_VAR ( in_filename );
    _LOG_VAR ( out_filename );
    return;
  }

  fstream infile ( in_filename, ios::in );

  fstream outfile ( out_filename, ios::out );

  {
    if ( infile.eof() )
      throw std::runtime_error ( "unable to read file" );

    getline ( infile, line );

    outfile << line << endl;

    linestream.str ( chompLine ( line ) );

    linestream >> num_atoms;

    if ( num_atoms <= 0 )
      throw std::runtime_error ( "zero atom info in file" );
  }

  uint num_atoms_read = 0;

  try
  {
    while ( !infile.eof() )
    {
      getline ( infile, line );

      stringstream linestream ( chompLine ( line ) );

      if ( linestream.str().size() == 0 )
      {
        outfile << line << endl;
        continue;
      }

      double epsilon, sigma, charge, asp;

      double x, y, z, r;

      string atom_name, acid_name, chain_name;

      uint atomno, acid_no;

      linestream >> atomno;

      if ( atomno - 1 != num_atoms_read )
        throw std::runtime_error ( "incorrectly numbered atom no in atom info " );

      linestream
          >> x
          >> y
          >> z
          >> r
          >> epsilon
          >> sigma
          >> charge
          >> asp
          >> atom_name
          >> acid_name
          >> chain_name
          >> acid_no;

      num_atoms_read++;

      char buf[200];

      double x_new, y_new, z_new;

      x_new = x * trans_mat[0] + y * trans_mat[4] + z * trans_mat[8] + trans_mat[12];

      y_new = x * trans_mat[1] + y * trans_mat[5] + z * trans_mat[9] + trans_mat[13];

      z_new = x * trans_mat[2] + y * trans_mat[6] + z * trans_mat[10] + trans_mat[14];

      atom_name.insert ( atom_name.end(), 4 - atom_name.size() , ' ' );

      acid_name.insert ( acid_name.end(), 4 - acid_name.size() , ' ' );

      chain_name.insert ( chain_name.end(), 4 - chain_name.size() , ' ' );

      sprintf ( buf, " %5d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f     %s       %s     %s       %5d\n",
                atomno, x_new, y_new, z_new, r,
                epsilon, sigma, charge, asp,
                atom_name.c_str(), acid_name.c_str(), chain_name.c_str(), acid_no );

      outfile << buf;

      if ( num_atoms_read == num_atoms )
        break;

    }

    if ( num_atoms_read != num_atoms )
      throw std::runtime_error ( "incorrect no of atoms read from the atom info" );

    while ( !infile.eof() )
    {
      getline ( infile, line );

      outfile << line << endl;
      continue;

    }

  }
  catch ( ... )
  {

  }
}


bool read_crd_file ( const char * filename, protein_t &protein )
{
  try
  {
    string line;

    stringstream linestream;

    uint num_atoms = 0;

    fstream crdfile ( filename, ios::in );

    {
      if ( crdfile.eof() )
        throw std::runtime_error ( "unable to read file" );

      getline ( crdfile, line );

      linestream.str ( chompLine ( line ) );

      linestream >> num_atoms;

      if ( num_atoms <= 0 )
        throw std::runtime_error ( "zero atom info in file" );
    }


    atom_t *atoms       = new atom_t[num_atoms];

    string *atom_names  = new string[num_atoms];

    string *acid_names  = new string[num_atoms];

    string *chain_names = new string[num_atoms];

    double *atom_rads   = new double[num_atoms];

    uint   *acid_nos    = new uint[num_atoms];

    uint   *bonds  = NULL;

    uint    num_bonds = 0;

    vector<uint> *bonds_vec = new vector<uint>;

    uint num_atoms_read = 0;

    try
    {
      while ( !crdfile.eof() )
      {
        getline ( crdfile, line );

        stringstream linestream ( chompLine ( line ) );

        linestream.seekg ( 0, ios::beg );

        if ( linestream.str().size() == 0 )
          continue;

        double epsilon, sigma, charge, asp;

        uint atomno;

        linestream >> atomno;

        if ( atomno - 1 != num_atoms_read )
          throw std::runtime_error ( "incorrectly numbered atom no in atom info " );

        linestream

            >> atoms[num_atoms_read].x
            >> atoms[num_atoms_read].y
            >> atoms[num_atoms_read].z
            >> atom_rads[num_atoms_read]
            >> epsilon
            >> sigma
            >> charge
            >> asp
            >> atom_names[num_atoms_read]
            >> acid_names[num_atoms_read]
            >> chain_names[num_atoms_read]
            >> acid_nos[num_atoms_read];

        atoms[num_atoms_read].radius = atom_rads[num_atoms_read];

        num_atoms_read++;

        if ( num_atoms_read == num_atoms )
          break;

      }

      if ( num_atoms_read != num_atoms )
        throw std::runtime_error ( "incorrect no of atoms read from the atom info" );

      num_atoms_read = 0;

      while ( !crdfile.eof() )
      {
        getline ( crdfile, line );

        stringstream linestream ( chompLine ( line ) );

        if ( linestream.str().size() == 0 )
          continue;

        uint atomno, atom_bond_ct;

        linestream >> atomno;

        if ( atomno - 1 != num_atoms_read )
          throw std::runtime_error ( "incorrectly numbered atom no in bond info" );

        linestream >> atom_bond_ct;

        atoms[num_atoms_read].bond_start = bonds_vec->size();

        atoms[num_atoms_read].bond_end   = bonds_vec->size() + 2 * atom_bond_ct;

        for ( uint i = 0 ; i < atom_bond_ct;i++ )
        {
          uint bonded_atom_no = 0;

          linestream >> bonded_atom_no;

          bonds_vec->push_back ( num_atoms_read );

          bonds_vec->push_back ( bonded_atom_no - 1 );
        }

        num_atoms_read++;

        if ( num_atoms_read == num_atoms )
          break;
      }

      crdfile.close();

      num_bonds = bonds_vec->size() / 2;

      bonds     = new uint[num_bonds*2];

      copy ( bonds_vec->begin(), bonds_vec->end(), bonds );

      if ( num_atoms_read != num_atoms )
        throw std::runtime_error ( "incorrect no of atoms read from in bond info" );

      delete bonds_vec;


    }
    catch ( ... )
    {

      delete[] atoms;
      delete[] atom_names;
      delete[] acid_names;
      delete[] chain_names;
      delete[] atom_rads;
      delete bonds_vec;

      throw;
    }


    uint num_acids  = 0;

    acid_t  * acids  = NULL;

    {
      num_acids = 1;

      for ( uint i = 1 ; i < num_atoms;i++ )
        if ( acid_names[i] != acid_names[i-1] ||
             chain_names[i] != chain_names[i-1] ||
             acid_nos[i] != acid_nos[i-1]  )
          num_acids ++;

      acids  = new acid_t [num_acids];

      uint idx = 0;

      acids[idx].start = 0;

      for ( uint i = 1 ; i < num_atoms;i++ )
        if ( acid_names[i] != acid_names[i-1] ||
             chain_names[i] != chain_names[i-1] ||
             acid_nos[i] != acid_nos[i-1]  )
        {
        acids[idx].end = i;
        ++idx;
        acids[idx].start = i;
      }

      acids[idx].end = num_atoms;
    }

    uint        num_chains = 0;

    chain_t * chains = NULL;

    {
      num_chains = 1;

      for ( uint i = 1 ; i < num_acids;i++ )
        if ( chain_names[acids[i].start] != chain_names[acids[i-1].start] )
          num_chains ++;

      chains = new chain_t[num_chains];

      uint idx = 0;

      chains[idx].start = 0;

      for ( uint i = 1 ; i < num_acids;i++ )
        if ( chain_names[acids[i].start] != chain_names[acids[i-1].start]  )
        {
        chains[idx].name = chain_names[acids[i-1].start];

        chains[idx].end = i;
        ++idx;
        chains[idx].start = i;
      }

      chains[idx].name = chain_names[acids[num_acids-1].start];

      chains[idx].end = num_acids;
    }



    uint      num_atom_types = 0;

    atom_type_t * atom_types = NULL;

    {
      map<string, uint>          atom_name_idx_map;

      for ( uint i = 0 ; i < num_atoms;i++ )
        if ( atom_name_idx_map.find ( atom_names[i] ) == atom_name_idx_map.end() )
          atom_name_idx_map[atom_names[i]] = i;

      num_atom_types = atom_name_idx_map.size();

      atom_types     = new atom_type_t[num_atom_types];

      uint idx = 0;

      for ( map<string, uint>::iterator iter = atom_name_idx_map.begin();
      iter != atom_name_idx_map.end();++iter, ++idx )
      {
        atom_types[idx].name         = iter->first;
        atom_types[idx].radius       = atom_rads[iter->second];
        atoms[iter->second].type_idx = idx;
        atom_types[idx].std_type_idx = get_std_atom_type_idx ( atom_types[idx].name );
      }

      for ( uint i = 0 ; i < num_atoms;i++ )
      {
        atoms[i].type_idx = atoms[atom_name_idx_map[atom_names[i]]].type_idx;
      }

    }


    uint                num_acid_types = 0;
    acid_type_t * acid_types = NULL;

    {
      map<string, uint> acid_name_idx_map;

      for ( uint i = 0 ; i < num_acids;i++ )
        if ( acid_name_idx_map.find ( acid_names[acids[i].start] ) ==
             acid_name_idx_map.end() )
          acid_name_idx_map[acid_names[acids[i].start]] = i;

      num_acid_types = acid_name_idx_map.size();

      acid_types     = new acid_type_t[num_acid_types];

      uint idx = 0;

      for ( map<string, uint>::iterator iter = acid_name_idx_map.begin();
      iter != acid_name_idx_map.end();++iter, ++idx )
      {
        acid_types[idx].name     = iter->first;
        acid_types[idx].inst_idx = iter->second;

        acids[iter->second].type_idx = idx;
      }

      for ( uint i = 0 ; i < num_acids;i++ )
        acids[i].type_idx = acids[acid_name_idx_map[acid_names[acids[i].start]]].type_idx;
    }


    delete[] atom_names;
    delete[] acid_names;
    delete[] chain_names;
    delete[] atom_rads;

    protein.atoms           = atoms;
    protein.num_atoms       = num_atoms;
    protein.acids           = acids;
    protein.num_acids       = num_acids;
    protein.chains          = chains;
    protein.num_chains      = num_chains;
    protein.atom_types      = atom_types;
    protein.num_atom_types  = num_atom_types;
    protein.acid_types      = acid_types;
    protein.num_acid_types  = num_acid_types;
    protein.bonds           = bonds;
    protein.num_bonds       = num_bonds;
    protein.acid_rel_bb_atoms_idx[0] = 0;
    protein.acid_rel_bb_atoms_idx[1] = 2;
    protein.acid_rel_bb_atoms_idx[2] = 4;
    protein.acid_rel_bb_atoms_idx[3] = 5;

    return true;

  }
  catch ( std::runtime_error e )
  {
    _ERROR ( e.what() );
    return false;
  }
}

bool begins_with(const string &str,const string &beg)
{
  string tmp(str.begin(),str.begin() + min(str.size(),beg.size()));

  return tmp == beg;
}

bool read_pdb_file ( const char *filename, protein_t & protein )
{

  static const uint atomline_parts_idx [] =
  {
    0,   //  Record name   "ATOM  "
    6,   //  Integer       serial       Atom  serial number.
    12,   //  Atom          name         Atom name.
    16,   //  Character     altLoc       Alternate location indicator.
    17,   //  Residue name  resName      Residue name.
    21,   //  Character     chainID      Chain identifier.
    22,   //  Integer       resSeq       Residue sequence number.
    26,   //  AChar         iCode        Code for insertion of residues.
    30,   //  Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    38,   //  Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    46,   //  Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    54,   //  Real(6.2)     occupancy    Occupancy.
    60,   //  Real(6.2)     tempFactor   Temperature  factor.
    76,   //  LString(2)    element      Element symbol, right-justified.
    78,   //  LString(2)    charge       Charge  on the atom.
    80
  };

  enum eATOMLINE_PARTS
  {
    ALP_REC_NAME,
    ALP_ATOM_SL_NO,
    ALP_ATOM_NAME,
    ALP_ATOM_ALT_LOC,
    ALP_ACID_NAME,
    ALP_CHAIN_ID,
    ALP_ACID_SEQ_NO,
    ALP_ACID_ICODE,
    ALP_ATOM_X,
    ALP_ATOM_Y,
    ALP_ATOM_Z,
    ALP_OCCUPANCY,
    ALP_TEMP_FACTOR,
    ALP_ELEMENT,
    ALP_CHARGE,
    ALP_COUNT
  };

  static const uint sheet_indexes_in_pdb[]=
  {
      0,//SHEET
      7,//strand no
      11,//sheet id
      14,//no of strands
      17,//init residue name
      21,//init chain id
      22,//residue seq no
      26,//code for insertion of residues
      28,//terminal residue no
      32,//term chain id
      33,//residue seq no
      37,//code for insertion of residues
      38 //strand sense wrt previous
      //rest is related to hydrogen bond...not required
  };

  static const uint helix_indexes_in_pdb[]=
  {
      0,//helix
      7,//helix serial no
      11,//helix id
      15,//initial residue name
      19,//init chain id
      21,//residue seq no
      25,//code for insertion of residues
      27,//term residue name
      31,//chain id
      33,//residue seq no
      37,//code for insertion of residues
      38,//type of helix
      40,//comment
      71,//lenght of helix
      76//eof
  };

  try
  {

    fstream pdbfile ( filename, ios::in );

    if ( pdbfile.eof() )
      throw std::runtime_error ( "unable to read file" );

    vector<string> atom_name_vec;

    vector<string> acid_name_vec;

    vector<uint>   acid_no_vec;

    vector<char>   chain_id_vec;

    vector<atom_t> atom_vec;

    vector<string> sheet_lines;
    vector<string> helix_lines;

    while ( !pdbfile.eof() )
    {
      char atomline[100];

      pdbfile.getline ( atomline, 100,'\n' );

      atomline[99] = '\0';

      if(begins_with(atomline,"SHEET"))
      {
        sheet_lines.push_back(atomline);
        continue;
      }
      else if(begins_with(atomline,"HELIX"))
      {
        helix_lines.push_back(atomline);
        continue;
      }
      else if(!begins_with(atomline,"ATOM"))
      {
        continue;
      }


      string atom_name ( atomline + atomline_parts_idx[ALP_ATOM_NAME],
                         atomline + atomline_parts_idx[ALP_ATOM_NAME+1] );

      string acid_name ( atomline + atomline_parts_idx[ALP_ACID_NAME],
                         atomline + atomline_parts_idx[ALP_ACID_NAME+1] );

      string seqno_str ( atomline + atomline_parts_idx[ALP_ACID_SEQ_NO],
                         atomline + atomline_parts_idx[ALP_ACID_SEQ_NO+1] );

      string xyz_str ( atomline + atomline_parts_idx[ALP_ATOM_X],
                       atomline + atomline_parts_idx[ALP_ATOM_Z+1] );

      atom_name_vec.push_back ( stripWS ( atom_name ) );

      acid_name_vec.push_back ( stripWS ( acid_name ) ) ;

      chain_id_vec.push_back ( atomline[atomline_parts_idx[ALP_CHAIN_ID]] );

      uint seqno;

      stringstream ( stripWS ( seqno_str ) ) >> seqno;

      acid_no_vec.push_back ( seqno );

      stringstream xyzstream ( stripWS ( xyz_str ) );

      atom_t atom;

      xyzstream >> atom.x >> atom.y >> atom.z;

      atom.bond_start = 0;

      atom.bond_end   = 0;

      atom_vec.push_back ( atom );
    }

    pdbfile.close();

    uint num_atoms = atom_vec.size();

    if ( ( atom_name_vec.size() != num_atoms ) ||
         ( acid_name_vec.size() != num_atoms ) ||
         ( acid_no_vec.size()   != num_atoms ) ||
         ( chain_id_vec.size()  != num_atoms )
      )
      throw std::runtime_error ( "list contents inconsistent in size" );

    atom_t * atoms = new atom_t[num_atoms];

    try
    {

      copy ( atom_vec.begin(), atom_vec.end(), atoms );

      atom_vec.clear();

      int   num_acids  = 0;
      int   num_chains = 0;
      acid_t  * acids  = NULL;
      chain_t * chains = NULL;

      map<char,int>         pdbchain_to_chain;
      vector<map<int,int> > pdbresno_to_resno;
      pdbresno_to_resno.resize(100);

      {
        num_acids  = 1;
        num_chains = 1;

        pdbchain_to_chain[chain_id_vec[0]]              = num_chains-1;
        pdbresno_to_resno[num_chains-1][acid_no_vec[0]] = num_acids -1;

        for ( uint i = 1 ; i < num_atoms;i++ )
        {

          if(chain_id_vec[i]  != chain_id_vec[i-1])
          {
            num_chains++;
            pdbchain_to_chain[chain_id_vec[i]] = num_chains-1;
          }

          if ( acid_no_vec[i]  != acid_no_vec[i-1]  ||
               chain_id_vec[i] != chain_id_vec[i-1])
          {
            num_acids++;
            pdbresno_to_resno[num_chains-1][acid_no_vec[i]] = num_acids -1;
          }
        }

        chains = new chain_t[num_chains];
        acids  = new acid_t [num_acids];


        int acid_idx = 0;
        int chain_idx = 0;

        acids[acid_idx].start   = 0;
        chains[chain_idx].start = 0;

        for ( int i = 1 ; i < num_atoms;i++ )
        {
          if ( acid_no_vec[i]  != acid_no_vec[i-1]  ||
               chain_id_vec[i] != chain_id_vec[i-1])
          {
            acids[acid_idx].end = i;
            ++acid_idx;
            acids[acid_idx].start = i;
          }

          if(chain_id_vec[i]  != chain_id_vec[i-1])
          {
            chains[chain_idx].end = acid_idx;
            chain_idx++;
            chains[chain_idx].start = acid_idx;
          }

        }

        chains[chain_idx].end = num_acids;
        acids[acid_idx].end = num_atoms;
      }


      uint      num_atom_types = 0;

      atom_type_t * atom_types = NULL;

      {
        map<string, uint>          atom_name_idx_map;

        for ( uint i = 0 ; i < num_atoms;i++ )
          if ( atom_name_idx_map.find ( atom_name_vec[i] ) == atom_name_idx_map.end() )
            atom_name_idx_map[atom_name_vec[i]] = i;

        num_atom_types = atom_name_idx_map.size();

        atom_types     = new atom_type_t[num_atom_types];

        uint idx = 0;

        for ( map<string, uint>::iterator iter = atom_name_idx_map.begin();
        iter != atom_name_idx_map.end();++iter, ++idx )
        {
          atom_types[idx].name         = iter->first;
          atoms[iter->second].type_idx = idx;
          uint std_type_idx            = get_std_atom_type_idx ( atom_types[idx].name );
          atom_types[idx].std_type_idx = std_type_idx;
          atom_types[idx].radius       = protein_t::std_atom_types[std_type_idx].vWall_radius;
        }

        for ( uint i = 0 ; i < num_atoms;i++ )
        {
          atoms[i].type_idx = atoms[atom_name_idx_map[atom_name_vec[i]]].type_idx;
          atoms[i].radius   = atom_types[atoms[i].type_idx].radius;
        }

      }


      uint                num_acid_types = 0;
      acid_type_t * acid_types = NULL;

      {
        map<string, uint> acid_name_idx_map;

        for ( uint i = 0 ; i < num_acids;i++ )
          if ( acid_name_idx_map.find ( acid_name_vec[acids[i].start] ) ==
               acid_name_idx_map.end() )
            acid_name_idx_map[acid_name_vec[acids[i].start]] = i;

        num_acid_types = acid_name_idx_map.size();

        acid_types     = new acid_type_t[num_acid_types];

        uint idx = 0;

        for ( map<string, uint>::iterator iter = acid_name_idx_map.begin();
        iter != acid_name_idx_map.end();++iter, ++idx )
        {
          acid_types[idx].name     = iter->first;
          acid_types[idx].inst_idx = iter->second;

          acids[iter->second].type_idx = idx;
        }

        for ( uint i = 0 ; i < num_acids;i++ )
          acids[i].type_idx = acids[acid_name_idx_map[acid_name_vec[acids[i].start]]].type_idx;
      }

      {
        protein.num_sheets = sheet_lines.size();

        if (protein.get_num_sheets() != 0)
          protein.sheets = new sheet_t[protein.get_num_sheets()];

        map<string,int> sheetid_to_sheetno;

        for(int i = 0 ; i < protein.get_num_sheets();++i)
        {
          string line = sheet_lines[i];

          string strand_no    = stripWS(string(line.begin() + sheet_indexes_in_pdb[1] ,line.begin() + sheet_indexes_in_pdb[2] ));
          string sheet_id     = stripWS(string(line.begin() + sheet_indexes_in_pdb[2] ,line.begin() + sheet_indexes_in_pdb[3] ));
          string num_strands  = stripWS(string(line.begin() + sheet_indexes_in_pdb[3] ,line.begin() + sheet_indexes_in_pdb[4] ));
          string start_chainno= stripWS(string(line.begin() + sheet_indexes_in_pdb[5] ,line.begin() + sheet_indexes_in_pdb[6] ));
          string start_resno  = stripWS(string(line.begin() + sheet_indexes_in_pdb[6] ,line.begin() + sheet_indexes_in_pdb[7] ));
          string end_chainno  = stripWS(string(line.begin() + sheet_indexes_in_pdb[9] ,line.begin() + sheet_indexes_in_pdb[10]));
          string end_resno    = stripWS(string(line.begin() + sheet_indexes_in_pdb[10],line.begin() + sheet_indexes_in_pdb[11]));

          if(sheetid_to_sheetno.count(sheet_id) == 0)
            sheetid_to_sheetno.insert(make_pair(sheet_id,(int)sheetid_to_sheetno.size()));

          assert(start_chainno.size() == 1);
          assert(end_chainno.size()   == 1);

          sheet_t &strand = protein.sheets[i];

          strand.strand_no     = atoi(strand_no.c_str());
          strand.sheet_no      = sheetid_to_sheetno[sheet_id];
          strand.num_strands   = atoi(num_strands.c_str());

          assert(pdbchain_to_chain.count(start_chainno[0]) == 1);
          assert(pdbchain_to_chain.count(end_chainno[0]) == 1);

          strand.start_chainno = pdbchain_to_chain[start_chainno[0]];
          strand.end_chainno   = pdbchain_to_chain[end_chainno[0]]+1;

          assert(pdbresno_to_resno[strand.start_chainno].count(atoi(start_resno.c_str())) == 1);
          assert(pdbresno_to_resno[strand.end_chainno-1].count(atoi(end_resno.c_str())) == 1);

          strand.start_resno   = pdbresno_to_resno[strand.start_chainno][atoi(start_resno.c_str())];
          strand.end_resno     = pdbresno_to_resno[strand.end_chainno-1][atoi(end_resno.c_str())]+1;

        }

      }

      {
        protein.num_helices = helix_lines.size();

        if (protein.get_num_helices() != 0)
          protein.helices = new helix_t[protein.get_num_helices()];

        for(int i = 0 ; i < protein.get_num_helices();++i)
        {
          string line = helix_lines[i];

          helix_t &helix = protein.helices[i];

          string start_chainno= stripWS(string(line.begin() + helix_indexes_in_pdb[4] ,line.begin() + helix_indexes_in_pdb[5] ));
          string start_resno  = stripWS(string(line.begin() + helix_indexes_in_pdb[5] ,line.begin() + helix_indexes_in_pdb[6] ));
          string end_chainno  = stripWS(string(line.begin() + helix_indexes_in_pdb[8] ,line.begin() + helix_indexes_in_pdb[9]));
          string end_resno    = stripWS(string(line.begin() + helix_indexes_in_pdb[9] ,line.begin() + helix_indexes_in_pdb[10]));

          assert(start_chainno.size() == 1);
          assert(end_chainno.size()   == 1);

          helix.start_chainno = pdbchain_to_chain[start_chainno[0]];          
          helix.end_chainno   = pdbchain_to_chain[end_chainno[0]]+1;

          helix.start_resno   = pdbresno_to_resno[helix.start_chainno][atoi(start_resno.c_str())];
          helix.end_resno     = pdbresno_to_resno[helix.end_chainno-1][atoi(end_resno.c_str())]+1;

        }

      }


      protein.atoms           = atoms;
      protein.num_atoms       = num_atoms;
      protein.acids           = acids;
      protein.num_acids       = num_acids;
      protein.chains          = chains;
      protein.num_chains      = num_chains;
      protein.atom_types      = atom_types;
      protein.num_atom_types  = num_atom_types;
      protein.acid_types      = acid_types;
      protein.num_acid_types  = num_acid_types;
      protein.bonds           = NULL;
      protein.num_bonds       = 0;
      protein.acid_rel_bb_atoms_idx[0] = 0;
      protein.acid_rel_bb_atoms_idx[1] = 1;
      protein.acid_rel_bb_atoms_idx[2] = 2;
      protein.acid_rel_bb_atoms_idx[3] = 3;
    }
    catch ( ... )
    {
      delete[] atoms;
      num_atoms = 0;
      atoms = NULL;
      throw;
    }
  }
  catch ( std::runtime_error e )
  {
    _ERROR ( e.what() );
    return false;
  }

  return true;
}
