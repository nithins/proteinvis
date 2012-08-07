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

#ifndef SECONDARYMODEL_H_INCLUDED
#define SECONDARYMODEL_H_INCLUDED

#include <vector>
#include <list>
#include <glutils.h>
#include <protein.h>
#include <oneLevelModel.h>
#include <QFile>
#include <QDir>
#include <QResource>
#include <logutil.h>
#include <map>


using namespace::std;

class protein_rd_t;

class secondary_model_t
{
  struct chain_rd_t
  {
    glutils::color_t       color;

    glutils::vertex_list_t spline_pts;
    glutils::bufobj_ptr_t  spline_pts_bo; //vb holds all points on the spline

    glutils::vertex_list_t sec_spline_pts;
    glutils::bufobj_ptr_t  sec_spline_pts_bo; //vb holds all points on the spline

  };

  struct helix_rd_t
  {
    glutils::color_t  color;

    int chainno;
    int spt_idx_b;
    int spt_idx_e;

    glutils::vertex_t  axis_b;
    glutils::vertex_t  axis_e;

    glutils::normal_t  axis_dir;
    glutils::normal_t  x_dir;
    glutils::normal_t  y_dir;
    double             radius;
  };

  struct strand_rd_t
  {
    int chainno;
    int spt_idx_b;
    int spt_idx_e;
    glutils::color_t     color;

    glutils::bufobj_ptr_t width_bo;
    vector<double>        width;
  };


private:

  //reference to the protein
  boost::shared_ptr<protein_t> m_protein;


  std::vector<chain_rd_t>  m_chains_rd;
  std::vector<helix_rd_t>  m_helices_rd;
  std::vector<strand_rd_t> m_strands_rd;


  //variables for spin box
  uint num_sheets;
  uint num_helices;

public:
  secondary_model_t(boost::shared_ptr<protein_t> protein_reference);
  ~secondary_model_t();

  //method to render the secondary protein structure
  void RenderSheets();
  void RenderTubes();
  void RenderHelices();
  void RenderImposterHelices();

  static void InitShaders();
  void InitSplines();
  void InitHelices();
  void InitSheets();

};
#endif
