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

#ifndef POCKETMODEL_H_INCLUDED
#define POCKETMODEL_H_INCLUDED

#include <vector>
#include <glutils.h>

class protein_rd_t;

const double g_add_atom_radius = 1.4;
const double g_alpha_step_value  = 0.1;


class pocket_model_t
{

  public:
    typedef unsigned int uint;

  private:

    struct pocket_tet_range_t
    {
      uint start;
      uint end;
    };

    struct alpha_pocket_range_t
    {
      uint start;
      uint end;
    };

    struct pockets_t
    {
      alpha_pocket_range_t     *alpha_pocket_ranges;
      uint                  num_alpha_pocket_ranges;

      pocket_tet_range_t       *pocket_tet_ranges;
      uint                  num_pocket_tet_ranges;

      uint *                    tet_idx;
      uint                  num_tet_idx;

      void init();

      void destroy();
    };

  private:
    void read_file ( const char * , const char * );

  public:

    pocket_model_t ( const char *, const char *, protein_rd_t * );

    ~pocket_model_t ();

    int render() const ;

    void setup_render ( const uint &alphanum, const uint &pocno ) ;

    inline void setup_render ( const uint &alphanum )
    {
      setup_render ( alphanum, ( uint ) - 1 );
    }

    inline double get_add_atom_radius()
    {
      return g_add_atom_radius;
    }

    inline double get_alpha_value ( const uint & alphanum )
    {
      return alphanum*g_alpha_step_value;
    }

    inline uint get_num_alpha()
    {
      return m_pockets.num_alpha_pocket_ranges;
    }

    inline uint get_num_pockets ( const uint & alphanum )
    {
      if ( !check_alpha_num ( alphanum ) )
        return 0;

      return ( m_pockets.alpha_pocket_ranges[alphanum].end -
               m_pockets.alpha_pocket_ranges[alphanum].start );
    }

    inline glutils::bufobj_ptr_t get_pocket_atoms()
    {
      return m_atom_indxs;
    }

    inline glutils::bufobj_ptr_t get_pocket_atom_bonds()
    {
      return m_atom_bonds;
    }

    void clear_render();


  private:



    bool check_alpha_num ( const uint & );

    protein_rd_t                 *m_protein_rd;
    glutils::renderable_t        *m_ren;
    pockets_t                     m_pockets;

    glutils::bufobj_ptr_t         m_atom_indxs;
    glutils::bufobj_ptr_t         m_atom_bonds;


};

#endif
