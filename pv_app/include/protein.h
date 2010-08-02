/***************************************************************************
 *   Copyright (C) 2009 by nithin   *
 *   nithin@vidyaranya   *
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
#ifndef __PROTEIN_H
#define __PROTEIN_H

#include <string>
#include <iostream>

typedef unsigned int uint;
typedef unsigned char uchar;


struct atom_t
{
  double x, y, z;

  uint   type_idx; // an index into the list of atom types

  double radius;

  uint   bond_start, bond_end; // index into a list of atom indices to which this atom is bonded to
};

struct acid_t
{
  uint type_idx; // an index into the list of amino acid types

  uint start, end;// indices into the list of atoms [)
};

struct chain_t
{
  std::string name;
  uint start, end; // indices into the list of amino acids [)
};

struct atom_type_t
{
  std::string name;
  double radius;
  uint  std_type_idx; //index into a list of standard items uint(-1) if no match found
};

struct acid_type_t
{
  std::string name;
  uint inst_idx; // an index of an instance of the type of amino acid in the list of amino acids
};

struct std_atom_types_t
{
  const char  * name;
  const char  * fullname;
  const uchar * color;
  const double ionic_radius;
  const double vWall_radius;
};

const uint NUM_BB_ATOMS = 4;

class protein_t
{

private:
  atom_t      *atoms;
  uint     num_atoms;

  uint        *bb_atoms_idx;
  uint     num_bb_atoms;

  uint        *bonds;
  uint     num_bonds;

  uint        *bb_bonds;
  uint     num_bb_bonds;

  acid_t      *acids;
  uint     num_acids;

  chain_t     *chains;
  uint     num_chains;

  atom_type_t *atom_types;
  uint     num_atom_types;

  acid_type_t *acid_types;
  uint     num_acid_types;

  uint         acid_rel_bb_atoms_idx[NUM_BB_ATOMS];// indices within an acid seq which contain the NCCO atoms

  void init();

  void destroy();

  void collect_bb_atoms();

  void collect_bb_bonds();

  bool is_backbone_atom ( uint );

  void check_only_four_bb_atoms_per_acid();


public:

  inline uint get_num_atoms() const
  {
    return num_atoms;
  }

  inline uint get_num_atom_types() const
  {
    return num_atom_types;
  }

  inline uint get_num_bonds() const
  {
    return num_bonds;
  }

  inline uint get_num_bb_atoms() const
  {
    return num_bb_atoms;
  }

  inline uint get_num_bb_bonds() const
  {
    return num_bb_bonds;
  }

  inline uint get_num_acids() const
  {
    return num_acids;
  }

  inline uint get_num_acid_types() const
  {
    return num_acid_types;
  }

  inline uint get_num_chains() const
  {
    return num_chains;
  }

  inline const atom_t * get_atoms() const
  {
    return atoms;
  }

  inline const atom_type_t * get_atom_types() const
  {
    return atom_types;
  }

  inline const uint* get_bonds() const
  {
    return bonds;
  }

  inline const uint* get_bb_atoms_idx() const
  {
    return bb_atoms_idx;
  }

  inline const uint* get_bb_bonds() const
  {
    return bb_bonds;
  }

  inline const acid_t * get_acids() const
  {
    return acids;
  }

  inline const acid_type_t * get_acid_types() const
  {
    return acid_types;
  }

  inline const chain_t * get_chains() const
  {
    return chains;
  }

  protein_t ( const char * );

  ~protein_t();

  void print_atoms ( std::ostream &o );

  void print_acids ( std::ostream &o );

  void print_chains ( std::ostream &o );

  void print_atom_types ( std::ostream &o );

  void print_acid_types ( std::ostream &o );

  void print_stats ( std::ostream &o );

  void print_info ( std::ostream &o );

  void print_incorrect_bonds ( std::ostream &o );

  void get_subset_atom_bonds
      ( const uint *, const uint &, uint *&, uint & ) const;

  uint get_atom_acid_idx(uint atomno);

  uint get_acid_chain_idx(uint atomno);

  static std_atom_types_t std_atom_types[];

  friend bool read_crd_file ( const char * , protein_t & );

  friend bool read_pdb_file ( const char * , protein_t & );

  static void transform_and_save_crd(const char *,const char *, const double*);
};

#include <glutils.h>

class protein_rd_t
{
  protein_t * m_protein;

  glutils::bufobj_ptr_t m_atom_coord_bo;

  glutils::bufobj_ptr_t m_atom_radii_bo;

  glutils::bufobj_ptr_t m_atom_bonds_bo;

  glutils::bufobj_ptr_t m_bb_atom_indices_bo;

  glutils::bufobj_ptr_t m_bb_bond_indices_bo;

  void upload_data_items();

  void compute_extent();

  double m_extent[6];

public:
  protein_rd_t ( protein_t * );

  ~protein_rd_t();

  void get_extent ( double * );

  inline protein_t* get_protein()
  {
    return m_protein;
  }

  glutils::bufobj_ptr_t get_coord_bo();

  glutils::bufobj_ptr_t get_radii_bo();

  glutils::bufobj_ptr_t get_bonds_bo();

  glutils::bufobj_ptr_t get_bb_coord_bo();

  glutils::bufobj_ptr_t get_bb_bonds_bo();
};

#include <cpputils.h>
class protein_grouping_t
{
public:

  typedef n_vector_t<double,3> atom_color_t;

  enum eGroupAtomsBy
  {GROUP_ATOMS_ATOM,
   GROUP_ATOMS_ATOM_TYPE,
   GROUP_ATOMS_ACID,
   GROUP_ATOMS_ACID_TYPE,
   GROUP_ATOMS_CHAIN,
   GROUP_ATOMS_ALL,
   GROUP_ATOMS_COUNT
 };

  protein_grouping_t(protein_t * p);
  ~protein_grouping_t();

  std::string        get_group_name(int groupno) const;
  inline uint        get_num_groups() const {return m_num_groups;}

  void               set_grouping_type(eGroupAtomsBy );

  atom_color_t       get_group_color(int groupno) const;
  void               set_group_color(atom_color_t,int);
  void               set_default_color_values();
  double*            create_atom_color_set();

  void               update_atom_color_bo();

  inline glutils::bufobj_ptr_t
      get_atom_color_bo() const
  {
    return m_atom_color_bo;
  }

private:

  uint                  get_updated_num_groups()const ;

  protein_t            *m_protein;
  eGroupAtomsBy         m_grouping_type;
  atom_color_t         *m_group_colors;
  int                   m_num_groups;
  glutils::bufobj_ptr_t m_atom_color_bo;

public:
  static std::string get_groupby_displaystr(eGroupAtomsBy g);

};

#endif//__PROTEIN_H
