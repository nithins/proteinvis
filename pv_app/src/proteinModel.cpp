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

#include <cstdlib>
#include <cfloat>
#include <fstream>

#include <logutil.h>

#include <proteinModel.h>
#include <oneLevelModel.h>
#include <protein.h>
#include <alphaModel.h>
#include <pocketModel.h>

using namespace std;

const double g_ball_stick_atom_radius = 0.4;

const double g_ball_stick_bond_radius = 0.2;

const double g_small_space_fill_atom_radius = g_ball_stick_atom_radius;

double s_extent[6];
double s_extent_valid = false;

glutils::material_properties_t g_atombond_material_default =
{
  {0,0,0,1},// ambient
  {0,0,0,1},// diffuse
  {1,1,1,1},// specular
  {0,0,0,1},// emission
  50      // shininess
};

glutils::material_properties_t g_surface_material_default =
{
  {0.2,0.2,0.2,1},// ambient
  {0.75,0.75,1.0,1},// diffuse
  {0,0,0,1},      // specular
  {0,0,0,1},      // emission
  1             // shininess
};

protein_model_t::protein_model_t (const string &pf) :
    m_render_model ( RMDL_SPACE_FILL ),
    m_render_mode ( RMDE_FULL ),
    m_surface_render_mode ( SRM_SHOW ),
    m_alpha_complex_render_mode ( ACRM_HIDE ),
    m_pocket_alpha_num ( 0 ),
    m_pocket_num ( 0 ),
    m_pocket_render_mode ( PRM_SHOW_NONE ),
    m_add_atom_radius ( 0.0 ),
    m_alpha_value ( 0.0 ),
    m_atombond_material(g_atombond_material_default),
    m_surface_material(g_surface_material_default)
{
  size_t last_slash_pos = pf.find_last_of ( "/\\" );

  if ( last_slash_pos != string::npos )
    m_protein_name = string ( pf.begin() + last_slash_pos + 1, pf.end() );
  else
    m_protein_name = pf;

  onelevel_model_t::init();

  m_protein.reset(new protein_t ( pf.c_str() ));

  m_protein_rd.reset(new protein_rd_t ( m_protein ));

  m_protein_atoms_grouping.reset(new protein_grouping_t(m_protein));

  if ( s_extent_valid == false )
  {
    m_protein_rd->get_extent ( s_extent );

    s_extent_valid = true;
  }
  else
  {
    double extent[6];
    m_protein_rd->get_extent ( extent );

    s_extent[0] = min ( extent[0], s_extent[0] );
    s_extent[1] = max ( extent[1], s_extent[1] );
    s_extent[2] = min ( extent[2], s_extent[2] );
    s_extent[3] = max ( extent[3], s_extent[3] );
    s_extent[4] = min ( extent[4], s_extent[4] );
    s_extent[5] = max ( extent[5], s_extent[5] );
  }
}

void protein_model_t::load_surface(const std::string &filename)
{
  glutils::vertex_list_t  vlist;
  glutils::tri_idx_list_t tlist;

  glutils::read_off_file ( filename.c_str(), vlist,tlist);

  double max_val = 0.0;

  for ( uint i = 0 ; i < m_protein->get_num_atoms();i++ )
  {
    max_val = max ( max_val, fabs ( m_protein->get_atoms() [i].x ) );
    max_val = max ( max_val, fabs ( m_protein->get_atoms() [i].y ) );
    max_val = max ( max_val, fabs ( m_protein->get_atoms() [i].z ) );
    max_val = max ( max_val, fabs ( m_protein->get_atoms() [i].radius ) );
  }

  uint inum = 0;

  while ( max_val > 10 )
  {
    max_val = max_val / 10;
    inum++;
  }

  double scalarfactor = 1;

  for ( uint i = 0;i < inum; i++ )
  {
    scalarfactor = scalarfactor * 10;
  }

  for ( uint i = 0;i < vlist.size(); i++ )
  {
    vlist[i] *= scalarfactor;
  }

  m_surface_renderer.reset
      (glutils::create_buffered_triangles_ren
       (  glutils::make_buf_obj(vlist), glutils::make_buf_obj(tlist),
          glutils::make_normals_buf_obj(vlist,tlist)));
}

void protein_model_t::load_alpha_shape(const std::string &filename)
{
  m_alpha_complex_model.reset
      (new alpha_complex_model_t( filename, m_protein_rd ));
}

void protein_model_t::load_pockets(const std::string &pocf, const std::string &tetf)
{
  m_pocket_model.reset(new pocket_model_t
                   ( pocf.c_str(), tetf.c_str() , m_protein_rd ));

  update_pocket_render_state();
}

void protein_model_t::render_onelevel() const
{

  glPushAttrib(GL_ENABLE_BIT);

  glEnable ( GL_COLOR_MATERIAL );

  glColorMaterial ( GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE ) ;

  m_atombond_material.render_all(GL_FRONT_AND_BACK);

  switch ( m_render_mode )
  {

  case RMDE_FULL:

    switch ( m_render_model )
    {

        case RMDL_SPACE_FILL:
      onelevel_model_t::render_sf
          ( m_protein_rd->get_coord_bo(),
            m_protein_atoms_grouping->get_atom_color_bo(),
            m_protein_rd->get_radii_bo(),
            m_add_atom_radius,
            m_alpha_value );
      break;


        case RMDL_BALL_STICK:
      onelevel_model_t::render_bs
          ( m_protein_rd->get_coord_bo(),
            m_protein_atoms_grouping->get_atom_color_bo(),
            g_ball_stick_atom_radius,
            m_protein_rd->get_bonds_bo(),
            g_ball_stick_bond_radius );
      break;

        case RMDL_SMALL_SPACE_FILL:
      onelevel_model_t::render_sf
          ( m_protein_rd->get_coord_bo(),
            m_protein_atoms_grouping->get_atom_color_bo(),
            g_small_space_fill_atom_radius );
      break;
    }

    break;

  case RMDE_BACKBONE:

    switch ( m_render_model )
    {

        case RMDL_SPACE_FILL:
      onelevel_model_t::render_sf
          ( m_protein_rd->get_coord_bo(),
            m_protein_atoms_grouping->get_atom_color_bo(),
            m_protein_rd->get_radii_bo(),
            m_protein_rd->get_bb_coord_bo() ,
            m_add_atom_radius,
            m_alpha_value );

      break;


        case RMDL_BALL_STICK:
      onelevel_model_t::render_bs
          ( m_protein_rd->get_coord_bo(),
            m_protein_atoms_grouping->get_atom_color_bo(),
            g_ball_stick_atom_radius,
            m_protein_rd->get_bb_coord_bo() ,
            m_protein_rd->get_bb_bonds_bo(),
            g_ball_stick_bond_radius );
      break;

        case RMDL_SMALL_SPACE_FILL:
      onelevel_model_t::render_sf
          ( m_protein_rd->get_coord_bo(),
            m_protein_atoms_grouping->get_atom_color_bo(),
            g_small_space_fill_atom_radius,
            m_protein_rd->get_bb_coord_bo() );

      break;

    }

    break;

  case RMDE_POCKET_ATOMS:

    if ( m_pocket_model != NULL )
    {

      switch ( m_render_model )
      {

          case RMDL_SPACE_FILL:
        onelevel_model_t::render_sf
            ( m_protein_rd->get_coord_bo(),
              m_protein_atoms_grouping->get_atom_color_bo(),
              m_protein_rd->get_radii_bo(),
              m_pocket_model->get_pocket_atoms(),
              m_add_atom_radius,
              m_alpha_value
              );

        break;

          case RMDL_BALL_STICK:
        onelevel_model_t::render_bs
            ( m_protein_rd->get_coord_bo(),
              m_protein_atoms_grouping->get_atom_color_bo(),
              g_ball_stick_atom_radius,
              m_pocket_model->get_pocket_atoms(),
              m_pocket_model->get_pocket_atom_bonds(),
              g_ball_stick_bond_radius );
        break;

          case RMDL_SMALL_SPACE_FILL:
        onelevel_model_t::render_sf
            ( m_protein_rd->get_coord_bo(),
              m_protein_atoms_grouping->get_atom_color_bo(),
              g_small_space_fill_atom_radius,
              m_pocket_model->get_pocket_atoms() );
        break;

      }
    }

    break;

  }

  glPopAttrib();
}

void protein_model_t::render_surface() const
{
  if ( m_surface_renderer != NULL )
  {
    glPushAttrib ( GL_POLYGON_BIT | GL_ENABLE_BIT );

    glDisable( GL_COLOR_MATERIAL );

    m_surface_material.render_all(GL_FRONT_AND_BACK);

    glEnable ( GL_RESCALE_NORMAL );

    switch ( m_surface_render_mode )
    {

    case SRM_SHOW_WIREFRAME:

      glPolygonMode ( GL_FRONT, GL_LINE );
      glPolygonMode ( GL_BACK, GL_LINE );

      m_surface_renderer->render();
      break;


    case SRM_SHOW:

      m_surface_renderer->render();
      break;
    };

    glPopAttrib();

  }
}

void protein_model_t::update_sf_model_for_pocket()
{
  if ( m_pocket_model == NULL )
    return;

  m_add_atom_radius = m_pocket_model->get_add_atom_radius();

  m_alpha_value     = m_pocket_model->get_alpha_value ( m_pocket_alpha_num );
}

void protein_model_t::update_pocket_render_state()
{
  if ( m_pocket_model == NULL )
    return;

  switch ( m_pocket_render_mode )
  {

  case PRM_SHOW_NONE:

    if ( m_render_mode == RMDE_POCKET_ATOMS )
      m_render_mode = RMDE_NONE;

    m_pocket_model->clear_render();

    break;

  case PRM_SHOW_ONE:
    m_pocket_model->setup_render ( m_pocket_alpha_num, m_pocket_num );

    break;

  case PRM_SHOW_ALL:
    m_pocket_model->setup_render ( m_pocket_alpha_num );

    break;
  }

  if ( m_render_mode == RMDE_POCKET_ATOMS )
  {
    update_sf_model_for_pocket();
  }
}

protein_model_t::~protein_model_t ()
{
}

int protein_model_t::render()
{
  double scalefactor;

  glPushMatrix();

  scalefactor = max ( max ( ( s_extent[1] - s_extent[0] ), ( s_extent[3] - s_extent[2] ) ), ( s_extent[5] - s_extent[4] ) ) / 2.0;

  glTranslated ( -1.0, -1.0, -1.0 );

  glScaled ( 1.0 / scalefactor, 1.0 / scalefactor, 1.0 / scalefactor );

  glTranslated ( -s_extent[0], -s_extent[2], -s_extent[4] );

  render_onelevel();

  if ( m_pocket_model != NULL )
    m_pocket_model->render();

  if ( m_alpha_complex_model != NULL )
  {
    if ( m_alpha_complex_render_mode&ACRM_TETRAHEDRONS )
      m_alpha_complex_model->render_tetrahedrons();

    if ( m_alpha_complex_render_mode&ACRM_TRIANGLES )
      m_alpha_complex_model->render_triangles();

    if ( m_alpha_complex_render_mode&ACRM_EDGES )
      m_alpha_complex_model->render_edges();
  }

  render_surface();

  glPopMatrix();

  return 0;
}

#include <QDir>
#include <QFileDialog>

protein_model_ui_t::protein_model_ui_t(protein_model_ptr_t pptr,QWidget * par)
{
  m_protein_model = pptr;

  setParent(par);

  setupUi(this);

  // its critical that no one emits any signals till m_protein_model::gl_init
  // has been called

  protein_grouping_comboBox->blockSignals(true);

  for(int i = 0;i< protein_grouping_t::GROUP_ATOMS_COUNT; ++i)
  {
    std::string grp_name = protein_grouping_t::get_groupby_displaystr((protein_grouping_t::eGroupAtomsBy)i);

    QString str(grp_name.c_str());

    protein_grouping_comboBox->addItem(str);
  }

  protein_grouping_comboBox->blockSignals(false);
}

void protein_model_ui_t::init_ui()
{
  m_protein_grouping_model =
      new protein_grouping_ui_model_t(m_protein_model->m_protein_atoms_grouping);

  protein_grouping_tableView->setModel(m_protein_grouping_model);

  if ( m_protein_model->m_surface_renderer == NULL )
    showhide_surface_groupBox->setEnabled ( false );
  else
  {
  }

  if ( m_protein_model->m_alpha_complex_model == NULL )
    alpha_complex_groupBox->setEnabled ( false );

  if ( m_protein_model->m_pocket_model == NULL )
  {
    pocatoms_rendermode_radioButton->setEnabled ( false );
    pocket_groupBox->setEnabled ( false );
  }
  else
  {
    pocket_alpha_num_spinBox->setMinimum ( 1 );
    pocket_alpha_num_spinBox->setMaximum
        ( m_protein_model->m_pocket_model->get_num_alpha() );

    pocket_num_spinBox->setMinimum ( 1 );
    pocket_num_spinBox->setMaximum
        ( m_protein_model->m_pocket_model->get_num_pockets
          ( m_protein_model->m_pocket_alpha_num ) );
  }
}

void protein_model_ui_t::pocket_ui_state_updated()
{

  // pull stuff from the ui and update the base class vars

  if ( !pocket_groupBox->isChecked() )
  {
    m_protein_model->m_pocket_render_mode  = protein_model_t::PRM_SHOW_NONE;
  }
  else
  {

    m_protein_model->m_pocket_alpha_num = pocket_alpha_num_spinBox->value() - 1;

    if ( !show_all_pockets_radioButton->isChecked() )
    {
      m_protein_model->m_pocket_num          = pocket_num_spinBox->value() - 1;

      m_protein_model->m_pocket_render_mode  = protein_model_t::PRM_SHOW_ONE;
    }
    else
    {
      m_protein_model->m_pocket_num          = ( uint ) - 1;

      m_protein_model->m_pocket_render_mode  = protein_model_t::PRM_SHOW_ALL;
    }
  }

  // intimate base class

  m_protein_model->update_pocket_render_state();

  // update ui

  pocket_num_spinBox->setMaximum ( m_protein_model->m_pocket_model->get_num_pockets ( m_protein_model->m_pocket_alpha_num ) );

  add_radius_sf_model_doubleSpinBox->setValue ( m_protein_model->m_add_atom_radius );

  alpha_value_sf_model_doubleSpinBox->setValue ( m_protein_model->m_alpha_value );

}

void protein_model_ui_t::alpha_complex_ui_state_updated()
{
  if ( !alpha_complex_groupBox->isChecked() )
  {
    m_protein_model->m_alpha_complex_render_mode = protein_model_t::ACRM_HIDE;
  }
  else
  {
//    if ( show_alpha_complex_tets_checkBox->isChecked() )
//      m_protein_model->m_alpha_complex_render_mode |= protein_model_t::ACRM_TETRAHEDRONS;
//    else
//      m_protein_model->m_alpha_complex_render_mode &= ( ( uint ) - 1 ) xor protein_model_t::ACRM_TETRAHEDRONS;

//    if ( show_alpha_complex_tris_checkBox->isChecked() )
//      m_protein_model->m_alpha_complex_render_mode |= protein_model_t::ACRM_TRIANGLES;
//    else
//      m_protein_model->m_alpha_complex_render_mode &= ( ( uint ) - 1 ) xor protein_model_t::ACRM_TRIANGLES;

//    if ( show_alpha_complex_edges_checkBox->isChecked() )
//      m_protein_model->m_alpha_complex_render_mode |= protein_model_t::ACRM_EDGES;
//    else
//      m_protein_model->m_alpha_complex_render_mode &= ( ( uint ) - 1 ) xor protein_model_t::ACRM_EDGES;
  }
}

void protein_model_ui_t::on_sf_model_radiobutton_clicked ( bool checked )
{
  if ( checked )
  {
    m_protein_model->m_render_model = protein_model_t::RMDL_SPACE_FILL;
  }
}

void protein_model_ui_t::on_ssf_model_radiobutton_clicked ( bool checked )
{
  if ( checked )
  {
    m_protein_model->m_render_model = protein_model_t::RMDL_SMALL_SPACE_FILL;
  }
}

void protein_model_ui_t::on_bs_model_radioButton_clicked ( bool checked )
{
  if ( checked )
  {
    m_protein_model->m_render_model = protein_model_t::RMDL_BALL_STICK;
  }
}

void protein_model_ui_t::on_full_mol_rendermode_radioButton_clicked ( bool checked )
{
  if ( checked )
  {
    m_protein_model->m_render_mode = protein_model_t::RMDE_FULL;
  }
}

void protein_model_ui_t::on_bb_mol_rendermode_radioButton_clicked ( bool checked )
{
  if ( checked )
  {
    m_protein_model->m_render_mode = protein_model_t::RMDE_BACKBONE;
  }
}

void protein_model_ui_t::on_pocatoms_rendermode_radioButton_clicked ( bool checked )
{
  if ( checked )
  {
    // change base variables
    m_protein_model->m_render_mode     = protein_model_t::RMDE_POCKET_ATOMS;

    // intimate base class
    m_protein_model->update_sf_model_for_pocket();

    // update ui
    add_radius_sf_model_doubleSpinBox->setValue ( m_protein_model->m_add_atom_radius );

    alpha_value_sf_model_doubleSpinBox->setValue ( m_protein_model->m_alpha_value );
  }
}

void protein_model_ui_t::on_none_rendermode_radioButton_clicked ( bool checked )
{
  if ( checked )
    m_protein_model->m_render_mode = protein_model_t::RMDE_NONE;
}

void protein_model_ui_t::on_add_radius_sf_model_doubleSpinBox_valueChanged ( double v )
{
  m_protein_model->m_add_atom_radius = v;
}

void protein_model_ui_t::on_alpha_value_sf_model_doubleSpinBox_valueChanged ( double v )
{
  m_protein_model->m_alpha_value = v;
}

void protein_model_ui_t::on_show_surface_radioButton_clicked ( bool checked )
{
  if ( checked )
    m_protein_model->m_surface_render_mode = protein_model_t::SRM_SHOW;
}

void protein_model_ui_t::on_show_surface_wireframe_radioButton_clicked ( bool checked )
{
  if ( checked )
    m_protein_model->m_surface_render_mode = protein_model_t::SRM_SHOW_WIREFRAME;
}

void protein_model_ui_t::on_hide_surface_radioButton_clicked ( bool checked )
{
  if ( checked )
    m_protein_model->m_surface_render_mode = protein_model_t::SRM_HIDE;
}

void protein_model_ui_t::on_pocket_groupBox_clicked ( bool )
{
  pocket_ui_state_updated();
}

void protein_model_ui_t::on_show_all_pockets_radioButton_clicked ( bool )
{
  pocket_ui_state_updated();
}

void protein_model_ui_t::on_pocket_alpha_num_spinBox_valueChanged ( int )
{
  pocket_ui_state_updated();
}

void protein_model_ui_t::on_pocket_num_spinBox_valueChanged ( int )
{
  pocket_ui_state_updated();
}

void protein_model_ui_t::on_alpha_complex_groupBox_clicked ( bool )
{
  alpha_complex_ui_state_updated();
}

void protein_model_ui_t::on_show_alpha_complex_tets_checkBox_clicked ( bool )
{
  alpha_complex_ui_state_updated();
}

void protein_model_ui_t::on_show_alpha_complex_tris_checkBox_clicked ( bool )
{
  alpha_complex_ui_state_updated();
}

void protein_model_ui_t::on_show_alpha_complex_edges_checkBox_clicked ( bool )
{
  alpha_complex_ui_state_updated();
}

void protein_model_ui_t::on_protein_grouping_comboBox_currentIndexChanged(int ind)
{
  m_protein_grouping_model->set_grouping_type(ind);
}

void protein_model_ui_t::on_reload_pushButton_clicked(bool)
{
  m_protein_model->m_protein_atoms_grouping->update_atom_color_bo();
}

void protein_model_ui_t::on_surface_color_colorpicker_colorChanged
    (const QColor &c)
{

}

protein_grouping_ui_model_t::protein_grouping_ui_model_t
    (boost::shared_ptr<protein_grouping_t> p, QObject *parent )
  :QAbstractItemModel(parent)
{
  m_protein_grouping = p;
}

protein_grouping_ui_model_t::~protein_grouping_ui_model_t()
{
}

QVariant protein_grouping_ui_model_t::data(const QModelIndex &index, int role) const
{
  if (!index.isValid() ||index.row() >= (int)m_protein_grouping->get_num_groups())
    return QVariant();

  switch(index.column())
  {
  case 0:
    {
      if (role == Qt::DisplayRole)
        return QString(m_protein_grouping->get_group_name(index.row()).c_str());
      break;
    }
  case 1:
    {
      if(role == Qt::DecorationRole ||
         role == Qt::DisplayRole)
      {
        glutils::color_t col = m_protein_grouping->get_group_color(index.row());
        QColor qcol;
        qcol.setRgbF(col[0],col[1],col[2]);
        return qcol;
      }
      break;
    }
  default:
    break;
  }
  return QVariant();
}

Qt::ItemFlags protein_grouping_ui_model_t::flags(const QModelIndex &index) const
{
  if (!index.isValid())
    return 0;

  switch(index.column())
  {
  case 0:
    return Qt::ItemIsEnabled | Qt::ItemIsSelectable;
  case 1:
    return Qt::ItemIsEnabled | Qt::ItemIsSelectable|Qt::ItemIsEditable;
  default:
    return 0;
  }
}

QVariant protein_grouping_ui_model_t::headerData
    (int section, Qt::Orientation orientation,int role) const

{
  if (orientation == Qt::Horizontal && role == Qt::DisplayRole)
  {
    switch(section)
    {
    case 0:
      return QString("Group");
    case 1:
      return QString("Color");
    }
  }

  return QVariant();
}


QModelIndex protein_grouping_ui_model_t::index
    (int row, int column,
     const QModelIndex &) const
{
  if(row >= (int)m_protein_grouping->get_num_groups())
    return QModelIndex();

  if(column >= 2)
    return QModelIndex();

  return createIndex(row,column,(void*)NULL);

}

QModelIndex protein_grouping_ui_model_t::parent(const QModelIndex &) const
{
  return QModelIndex();
}

int protein_grouping_ui_model_t::rowCount(const QModelIndex & ) const
{
  return m_protein_grouping->get_num_groups();
}

int protein_grouping_ui_model_t::columnCount(const QModelIndex &) const
{
  return 2;
}

bool protein_grouping_ui_model_t::setData
    ( const QModelIndex & index,
      const QVariant & value,
      int role )
{

  if(role == Qt::EditRole
     && index.column() == 1 &&
     index.row() >=0 &&
     index.row() < (int)m_protein_grouping->get_num_groups() &&
     value.canConvert<QColor>())
  {

    QColor color = value.value<QColor>();
    glutils::color_t col(color.redF(),color.greenF(),color.blueF());
    m_protein_grouping->set_group_color(col,index.row());
    emit dataChanged(index,index);
    return true;
  }

  return false;
}

void protein_grouping_ui_model_t::set_grouping_type(int ind)
{
  m_protein_grouping->set_grouping_type((protein_grouping_t::eGroupAtomsBy)ind);

  reset();
}

uint protein_grouping_ui_model_t::get_grouping_type()
{
  return m_protein_grouping->get_grouping_type();
}
