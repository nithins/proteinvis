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
#ifndef ONELEVELMODEL_H_INCLUDED
#define ONELEVELMODEL_H_INCLUDED

#include <glutils.h>

class protein_rd_t;

class GLSLProgram;

namespace glutils
{

  class buf_obj_t;

  class BufobjArrayRenderer;
}

typedef unsigned int uint;
typedef unsigned char uchar;

class onelevel_model_t
{

  public:

    onelevel_model_t ();

    ~onelevel_model_t ();



    inline void render_sf
    ( glutils::buf_obj_t atom_coord,
      glutils::buf_obj_t atom_color,
      glutils::buf_obj_t atom_radii,
      const double & add_radius = 0.0,
      const double & alpha = 0.0 ) const
    {
      render_spheres ( atom_coord, atom_color, atom_radii , add_radius, alpha );
    }

    inline void render_sf
    ( glutils::buf_obj_t atom_coord,
      glutils::buf_obj_t atom_color,
      glutils::buf_obj_t atom_radii,
      glutils::buf_obj_t atom_indxs,
      const double & add_radius = 0.0,
      const double & alpha = 0.0 ) const
    {
      render_spheres ( atom_coord, atom_color, atom_radii, atom_indxs, add_radius, alpha );
    }

    inline void render_sf
    ( glutils::buf_obj_t atom_coord,
      glutils::buf_obj_t atom_color,
      double atom_radius ) const
    {
      render_spheres ( atom_coord, atom_color, atom_radius );
    }

    inline void render_sf
    ( glutils::buf_obj_t atom_coord,
      glutils::buf_obj_t atom_color,
      double atom_radius,
      glutils::buf_obj_t atom_indxs ) const
    {
      render_spheres ( atom_coord, atom_color, atom_radius, atom_indxs );
    }

    inline void render_bs
    ( glutils::buf_obj_t atom_coord,
      glutils::buf_obj_t atom_color,
      double atom_radius,
      glutils::buf_obj_t atom_bonds,
      double bond_radius ) const
    {
      render_spheres ( atom_coord, atom_color, atom_radius );

      glColor3ub ( 0xE0, 0xB0, 0xFF  );

      render_cylinders ( atom_coord, atom_bonds, bond_radius );
    }

    inline void render_bs
    ( glutils::buf_obj_t atom_coord,
      glutils::buf_obj_t atom_color,
      double atom_radius,
      glutils::buf_obj_t atom_indxs,
      glutils::buf_obj_t atom_bonds,
      double bond_radius ) const
    {
      render_spheres ( atom_coord, atom_color, atom_radius, atom_indxs );

      glColor3ub ( 0xE0,0xB0,0xFF  );

      render_cylinders ( atom_coord, atom_bonds, bond_radius );
    }

    void render_sphere
    ( const double &,
      const double &,
      const double &,
      const double &
    );

    void render_cylinder
    ( const double &,
      const double &,
      const double &,
      const double &,
      const double &,
      const double &,
      const double &
    );

  private:

    void render_spheres
    ( glutils::buf_obj_t center_coord,
      glutils::buf_obj_t sphere_color,
      glutils::buf_obj_t sphere_radii,
      const double &add_radius ,
      const double &alpha ) const;

    void render_spheres
    ( glutils::buf_obj_t center_coord,
      glutils::buf_obj_t sphere_color,
      glutils::buf_obj_t sphere_radii,
      glutils::buf_obj_t sphere_indxs ,
      const double &add_radius,
      const double &alpha ) const;

    void render_spheres
    ( glutils::buf_obj_t center_coord ,
      glutils::buf_obj_t sphere_color,
      double sphere_radius ) const;

    void render_spheres
    ( glutils::buf_obj_t center_coord ,
      glutils::buf_obj_t sphere_color,
      double sphere_radius,
      glutils::buf_obj_t sphere_indxs ) const;

    void render_cylinders
    ( glutils::buf_obj_t endpts_coord,
      glutils::buf_obj_t endpts_indxs,
      double cyl_radius
    ) const;


    static GLSLProgram*  s_sphereShader;

    static GLSLProgram*  s_cylinderShader;
};

#endif
