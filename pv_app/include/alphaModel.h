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
#ifndef ALPHAMODEL_H_INCLUDED
#define ALPHAMODEL_H_INCLUDED

namespace glutils
{

  class BufobjArrayRenderer;

  class renderable_t;
}

class protein_rd_t;

class alpha_complex_model_t
{

  private:
    void read_file ( const char * );

  public:

    alpha_complex_model_t ( const char *, protein_rd_t * );

    ~alpha_complex_model_t ();

    int render_tetrahedrons() const;
    int render_triangles() const;
    int render_edges() const;

  private:

    protein_rd_t                 *m_protein_rd;
    glutils::renderable_t        *m_tri_ren;
    glutils::renderable_t        *m_ege_ren;
    glutils::renderable_t        *m_tet_ren;


};

#endif

