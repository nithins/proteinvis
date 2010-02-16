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

#include <modelcontroller.h>

#include <proteinModel.h>
#include <oneLevelModel.h>
#include <protein.h>
#include <alphaModel.h>
#include <pocketModel.h>

#include <QDir>
#include <QFileDialog>
#include <ui_ProteinModelFrame.h>



using namespace std;

enum eRenderModels
{
  RMDL_SPACE_FILL,
  RMDL_BALL_STICK,
  RMDL_SMALL_SPACE_FILL,
  RMDL_COUNT
};

enum eRenderModes
{
  RMDE_BACKBONE,
  RMDE_FULL,
  RMDE_POCKET_ATOMS,
  RMDE_NONE,
  RMDE_COUNT
};

enum eSurfaceRenderMode
{
  SRM_HIDE,
  SRM_SHOW,
  SRM_SHOW_WIREFRAME,
  SRM_COUNT
};

enum eAlphaComplexRenderMode
{
  ACRM_HIDE         = 0,
  ACRM_TETRAHEDRONS = 1,
  ACRM_TRIANGLES    = 2,
  ACRM_EDGES        = 4,
};

enum ePocketRenderMode
{
  PRM_SHOW_NONE,
  PRM_SHOW_ONE,
  PRM_SHOW_ALL,
};

const double g_ball_stick_atom_radius = 0.4;

const double g_ball_stick_bond_radius = 0.2;

const double g_small_space_fill_atom_radius = g_ball_stick_atom_radius;

const unsigned char g_default_surface_color  [] = {176,148,101};

double s_extent[6];
double s_extent_valid = false;

ProteinModel::ProteinModel ( string pf, string sf, string acf , string pocf, string tetf ) :
    m_onelevel_model ( NULL ),
    m_protein ( NULL ),
    m_surface_renderer ( NULL ),
    m_alpha_complex_model ( NULL ),
    m_pocket_model ( NULL ),
    m_render_model ( RMDL_SPACE_FILL ),
    m_render_mode ( RMDE_FULL ),
    m_surface_render_mode ( SRM_SHOW ),
    m_alpha_complex_render_mode ( ACRM_HIDE ),
    m_pocket_alpha_num ( 0 ),
    m_pocket_num ( 0 ),
    m_pocket_render_mode ( PRM_SHOW_NONE ),
    m_add_atom_radius ( 0.0 ),
    m_alpha_value ( 0.0 ),
    m_protein_filename ( pf ),
    m_surface_filename ( sf ),
    m_alpha_complex_filename ( acf ),
    m_pocket_filename ( pocf ),
    m_tetra_filename ( tetf )
{
  m_controller = IModelController::Create();

  m_protein = new protein_t ( m_protein_filename.c_str() );

  m_protein_rd = new protein_rd_t ( m_protein );

  m_protein_atoms_grouping = new protein_grouping_t(m_protein);

  m_onelevel_model = new onelevel_model_t ();

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

  size_t last_slash_pos = m_protein_filename.find_last_of ( "/\\" );

  if ( last_slash_pos != string::npos )
    m_protein_name = string ( m_protein_filename.begin() + last_slash_pos + 1, m_protein_filename.end() );
  else
    m_protein_name = m_protein_filename;

  setup_surface();

  if ( m_alpha_complex_filename.size() != 0 )
  {
    m_alpha_complex_model = new alpha_complex_model_t
                            ( m_alpha_complex_filename.c_str(), m_protein_rd );
  }

  if ( m_pocket_filename.size() != 0 && m_tetra_filename.size() != 0 )
  {
    m_pocket_model = new pocket_model_t
                     ( m_pocket_filename.c_str(), m_tetra_filename.c_str() , m_protein_rd );

    update_pocket_render_state();
  }

  init_ui();
}

void ProteinModel::setup_surface()
{
  if ( m_surface_filename.size() == 0 )
    return;

  double *vrts;

  uint   *tris;

  uint    vrt_ct, tri_ct;

  copy(g_default_surface_color,g_default_surface_color+3,m_surface_color.elements);

  glutils::read_off_file ( m_surface_filename.c_str(), vrts, vrt_ct, tris, tri_ct );

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

  for ( uint i = 0;i < 3*vrt_ct; i++ )
  {
    vrts[i] *= scalarfactor;
  }

  glutils::bufobj_ptr_t vrt_bo,tri_bo,col_bo;

  vrt_bo  = glutils::buf_obj_t::create_bo
            (vrts, GL_DOUBLE, 3, GL_ARRAY_BUFFER,sizeof( double )*3*vrt_ct, 0 );

  tri_bo = glutils::buf_obj_t::create_bo
           (tris , GL_UNSIGNED_INT, 3, GL_ELEMENT_ARRAY_BUFFER,
            sizeof ( uint ) *3*tri_ct, 0 );

  col_bo = glutils::buf_obj_t::create_bo();

  m_surface_renderer = glutils::create_buffered_tristrip_ren ( vrt_bo, tri_bo, col_bo );

  delete []vrts;
  delete []tris;

}

void ProteinModel::render_onelevel() const
{
  switch ( m_render_mode )
  {

  case RMDE_FULL:

    switch ( m_render_model )
    {

        case RMDL_SPACE_FILL:
      m_onelevel_model->render_sf
          ( m_protein_rd->get_coord_bo(),
            m_protein_atoms_grouping->get_atom_color_bo(),
            m_protein_rd->get_radii_bo(),
            m_add_atom_radius,
            m_alpha_value );
      break;


        case RMDL_BALL_STICK:
      m_onelevel_model->render_bs
          ( m_protein_rd->get_coord_bo(),
            m_protein_atoms_grouping->get_atom_color_bo(),
            g_ball_stick_atom_radius,
            m_protein_rd->get_bonds_bo(),
            g_ball_stick_bond_radius );
      break;

        case RMDL_SMALL_SPACE_FILL:
      m_onelevel_model->render_sf
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
      m_onelevel_model->render_sf
          ( m_protein_rd->get_coord_bo(),
            m_protein_atoms_grouping->get_atom_color_bo(),
            m_protein_rd->get_radii_bo(),
            m_protein_rd->get_bb_coord_bo() ,
            m_add_atom_radius,
            m_alpha_value );

      break;


        case RMDL_BALL_STICK:
      m_onelevel_model->render_bs
          ( m_protein_rd->get_coord_bo(),
            m_protein_atoms_grouping->get_atom_color_bo(),
            g_ball_stick_atom_radius,
            m_protein_rd->get_bb_coord_bo() ,
            m_protein_rd->get_bb_bonds_bo(),
            g_ball_stick_bond_radius );
      break;

        case RMDL_SMALL_SPACE_FILL:
      m_onelevel_model->render_sf
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
        m_onelevel_model->render_sf
            ( m_protein_rd->get_coord_bo(),
              m_protein_atoms_grouping->get_atom_color_bo(),
              m_protein_rd->get_radii_bo(),
              m_pocket_model->get_pocket_atoms(),
              m_add_atom_radius,
              m_alpha_value
              );

        break;

          case RMDL_BALL_STICK:
        m_onelevel_model->render_bs
            ( m_protein_rd->get_coord_bo(),
              m_protein_atoms_grouping->get_atom_color_bo(),
              g_ball_stick_atom_radius,
              m_pocket_model->get_pocket_atoms(),
              m_pocket_model->get_pocket_atom_bonds(),
              g_ball_stick_bond_radius );
        break;

          case RMDL_SMALL_SPACE_FILL:
        m_onelevel_model->render_sf
            ( m_protein_rd->get_coord_bo(),
              m_protein_atoms_grouping->get_atom_color_bo(),
              g_small_space_fill_atom_radius,
              m_pocket_model->get_pocket_atoms() );
        break;

      }
    }

    break;

  }

}

void ProteinModel::update_sf_model_for_pocket()
{
  if ( m_pocket_model == NULL )
    return;

  m_add_atom_radius = m_pocket_model->get_add_atom_radius();

  m_alpha_value     = m_pocket_model->get_alpha_value ( m_pocket_alpha_num );
}

void ProteinModel::update_pocket_render_state()
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

// int ProteinModel::RenderForPick() const
// {
//   m_controller->Render();
//
//   DRAWBB ( -1, 1, -1, 1, -1, 1 );
//   return 24;
// }

ProteinModel::~ProteinModel ()
{
  _LOG("destroying Proteinmodel");

  destroy_ui();

  IModelController::Delete ( m_controller );

  delete m_onelevel_model;

  delete m_protein;

  delete m_protein_rd;

  delete m_protein_atoms_grouping;

  if ( m_surface_renderer != NULL )
    delete m_surface_renderer;

  if ( m_alpha_complex_model != NULL )
    delete m_alpha_complex_model;
}

void ProteinModel::Reset()
{
  m_controller->Reset();
}

int ProteinModel::Render() const
{
  m_controller->Render();

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

  if ( m_surface_renderer != NULL )
  {
    glPushAttrib ( GL_POLYGON_BIT | GL_ENABLE_BIT );

    glEnable ( GL_RESCALE_NORMAL );

    glColor3ubv ( m_surface_color.elements);

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

  glPopMatrix();

  return 0;
}

bool ProteinModel::MousePressedEvent
    ( const int &x, const int &y, const eMouseButton &mb,
      const eKeyFlags &,const eMouseFlags &)
{
  switch(mb)
  {
  case MOUSEBUTTON_RIGHT:
    {
      m_controller->StartTB ( x, y );
      return true;
    }
  case MOUSEBUTTON_LEFT:
    {
      m_controller->StartTrans ( x, y );
      return true;
    }
  default:
    return false;
  }
}

bool ProteinModel::MouseReleasedEvent
    ( const int &x, const int &y, const eMouseButton &mb,
      const eKeyFlags &,const eMouseFlags &)
{
  switch(mb)
  {
  case MOUSEBUTTON_RIGHT:
    {
      m_controller->StopTB ( x, y );
      return true;
    }

  case MOUSEBUTTON_LEFT:
    {
      m_controller->StopTrans( x, y );
      return true;
    }
  default:
    return false;
  }
}

bool ProteinModel::MouseMovedEvent
    ( const int &x, const int &y, const int &, const int &,
      const eKeyFlags &,const eMouseFlags &mf)
{
  if(mf &MOUSEBUTTON_LEFT || mf &MOUSEBUTTON_RIGHT)
  {
    m_controller->Move ( x, y );
    return true;
  }
  return false;
}

bool ProteinModel::WheelEvent
    ( const int &, const int &, const int &d,
      const eKeyFlags &kf,const eMouseFlags &)
{

  if(kf&KEYFLAG_CTRL)
  {

    double gf = 1.0;

    if(d>0)
      gf += 0.1;
    else
      gf -= 0.1;

    m_controller->set_uniform_scale(m_controller->get_uniform_scale()*gf);
    return true;
  }

  return false;
}

void ProteinModel::save_to_file ( std::string filename )
{
  double transmat[16];

  m_controller->getMatrix ( transmat );

  protein_t::transform_and_save_crd ( m_protein_filename.c_str(), filename.c_str(), transmat );
}

void ProteinModel::save_transformation(std::string filename )
{
  fstream outfile(filename.c_str(),ios::out);

  ostream * o = &outfile;

  m_controller->serialize(*o);

  outfile.close();

}
void ProteinModel::load_transformation(std::string filename )
{
  fstream infile(filename.c_str(),ios::in);

  istream * i = &infile;

  m_controller->serialize(*i);

  infile.close();
}

void ProteinModel::init_ui()
{
  setWindowTitle(QString("hello"));

  m_ui = new Ui::ProteinModelFrame();

  m_ui->setupUi ( this );

  m_protein_grouping_model =
      new protein_grouping_ui_model_t(m_protein_atoms_grouping);

  m_ui->protein_grouping_tableView->setModel(m_protein_grouping_model);

  for(int i = 0;i< protein_grouping_t::GROUP_ATOMS_COUNT; ++i)
  {
    m_ui->protein_grouping_comboBox->addItem
        (protein_grouping_t::get_groupby_displaystr((protein_grouping_t::eGroupAtomsBy)i).c_str() );
  }

  if ( m_surface_renderer == NULL )
    m_ui->showhide_surface_groupBox->setEnabled ( false );
  else
  {
    m_ui->surface_color_colorpicker->setCurrentColor
        (QColor(m_surface_color[0],m_surface_color[1],m_surface_color[2]));
  }


  if ( m_alpha_complex_model == NULL )
    m_ui->alpha_complex_groupBox->setEnabled ( false );

  if ( m_pocket_model == NULL )
  {
    m_ui->pocatoms_rendermode_radioButton->setEnabled ( false );
    m_ui->pocket_groupBox->setEnabled ( false );
  }
  else
  {
    m_ui->pocket_alpha_num_spinBox->setMinimum ( 1 );
    m_ui->pocket_alpha_num_spinBox->setMaximum ( m_pocket_model->get_num_alpha() );

    m_ui->pocket_num_spinBox->setMinimum ( 1 );
    m_ui->pocket_num_spinBox->setMaximum ( m_pocket_model->get_num_pockets ( m_pocket_alpha_num ) );
  }
}

void ProteinModel::destroy_ui()
{
  delete m_ui;

  delete m_protein_grouping_model;
}

void ProteinModel::pocket_ui_state_updated()
{

  // pull stuff from the ui and update the base class vars

  if ( !m_ui->pocket_groupBox->isChecked() )
  {
    m_pocket_render_mode  = PRM_SHOW_NONE;
  }
  else
  {

    m_pocket_alpha_num = m_ui->pocket_alpha_num_spinBox->value() - 1;

    if ( !m_ui->show_all_pockets_radioButton->isChecked() )
    {
      m_pocket_num          = m_ui->pocket_num_spinBox->value() - 1;

      m_pocket_render_mode  = PRM_SHOW_ONE;
    }
    else
    {
      m_pocket_num          = ( uint ) - 1;

      m_pocket_render_mode  = PRM_SHOW_ALL;
    }
  }

  // intimate base class

  update_pocket_render_state();

  // update ui

  m_ui->pocket_num_spinBox->setMaximum ( m_pocket_model->get_num_pockets ( m_pocket_alpha_num ) );

  m_ui->add_radius_sf_model_doubleSpinBox->setValue ( m_add_atom_radius );

  m_ui->alpha_value_sf_model_doubleSpinBox->setValue ( m_alpha_value );

}

void ProteinModel::alpha_complex_ui_state_updated()
{
  if ( !m_ui->alpha_complex_groupBox->isChecked() )
  {
    m_alpha_complex_render_mode = ACRM_HIDE;
  }
  else
  {
    if ( m_ui->show_alpha_complex_tets_checkBox->isChecked() )
      m_alpha_complex_render_mode |= ACRM_TETRAHEDRONS;
    else
      m_alpha_complex_render_mode &= ( ( uint ) - 1 ) xor ACRM_TETRAHEDRONS;

    if ( m_ui->show_alpha_complex_tris_checkBox->isChecked() )
      m_alpha_complex_render_mode |= ACRM_TRIANGLES;
    else
      m_alpha_complex_render_mode &= ( ( uint ) - 1 ) xor ACRM_TRIANGLES;

    if ( m_ui->show_alpha_complex_edges_checkBox->isChecked() )
      m_alpha_complex_render_mode |= ACRM_EDGES;
    else
      m_alpha_complex_render_mode &= ( ( uint ) - 1 ) xor ACRM_EDGES;
  }
}

void ProteinModel::on_sf_model_radiobutton_clicked ( bool checked )
{
  if ( checked )
  {
    m_render_model = RMDL_SPACE_FILL;
  }
}

void ProteinModel::on_ssf_model_radiobutton_clicked ( bool checked )
{
  if ( checked )
  {
    m_render_model = RMDL_SMALL_SPACE_FILL;
  }
}

void ProteinModel::on_bs_model_radioButton_clicked ( bool checked )
{
  if ( checked )
  {
    m_render_model = RMDL_BALL_STICK;
  }
}

void ProteinModel::on_full_mol_rendermode_radioButton_clicked ( bool checked )
{
  if ( checked )
  {
    m_render_mode = RMDE_FULL;
  }
}

void ProteinModel::on_bb_mol_rendermode_radioButton_clicked ( bool checked )
{
  if ( checked )
  {
    m_render_mode = RMDE_BACKBONE;
  }
}

void ProteinModel::on_pocatoms_rendermode_radioButton_clicked ( bool checked )
{
  if ( checked )
  {
    // change base variables
    m_render_mode     = RMDE_POCKET_ATOMS;

    // intimate base class
    update_sf_model_for_pocket();

    // update ui
    m_ui->add_radius_sf_model_doubleSpinBox->setValue ( m_add_atom_radius );

    m_ui->alpha_value_sf_model_doubleSpinBox->setValue ( m_alpha_value );
  }
}

void ProteinModel::on_none_rendermode_radioButton_clicked ( bool checked )
{
  if ( checked )
    m_render_mode = RMDE_NONE;
}

void ProteinModel::on_add_radius_sf_model_doubleSpinBox_valueChanged ( double v )
{
  m_add_atom_radius = v;
}

void ProteinModel::on_alpha_value_sf_model_doubleSpinBox_valueChanged ( double v )
{
  m_alpha_value = v;
}

void ProteinModel::on_show_surface_radioButton_clicked ( bool checked )
{
  if ( checked )
    m_surface_render_mode = SRM_SHOW;
}

void ProteinModel::on_show_surface_wireframe_radioButton_clicked ( bool checked )
{
  if ( checked )
    m_surface_render_mode = SRM_SHOW_WIREFRAME;
}

void ProteinModel::on_hide_surface_radioButton_clicked ( bool checked )
{
  if ( checked )
    m_surface_render_mode = SRM_HIDE;
}

void ProteinModel::on_pocket_groupBox_clicked ( bool )
{
  pocket_ui_state_updated();
}

void ProteinModel::on_show_all_pockets_radioButton_clicked ( bool )
{
  pocket_ui_state_updated();
}

void ProteinModel::on_pocket_alpha_num_spinBox_valueChanged ( int )
{
  pocket_ui_state_updated();
}

void ProteinModel::on_pocket_num_spinBox_valueChanged ( int )
{
  pocket_ui_state_updated();
}

void ProteinModel::on_alpha_complex_groupBox_clicked ( bool )
{
  alpha_complex_ui_state_updated();
}

void ProteinModel::on_show_alpha_complex_tets_checkBox_clicked ( bool )
{
  alpha_complex_ui_state_updated();
}

void ProteinModel::on_show_alpha_complex_tris_checkBox_clicked ( bool )
{
  alpha_complex_ui_state_updated();
}

void ProteinModel::on_show_alpha_complex_edges_checkBox_clicked ( bool )
{
  alpha_complex_ui_state_updated();
}

void ProteinModel::on_save_to_file_pushButton_clicked ( bool )
{
  QString filename = QFileDialog::getSaveFileName
                     ( this, tr ( "Save File" ),
                       QDir::homePath(),
                       tr ( "Protein data (*.crd)" ) );

  save_to_file ( filename.toAscii().constData() );
}

void ProteinModel::on_save_trans_pushButton_clicked(bool)
{
  QString filename = QFileDialog::getSaveFileName
                     ( this, tr ( "Save Transformation" ),
                       QDir::homePath(),
                       tr ( "Transformation data (*.trans)" ) );

  save_transformation(filename.toAscii().constData());

}
void ProteinModel::on_load_trans_pushButton_clicked(bool)
{
  QString filename = QFileDialog::getOpenFileName
                     ( this, tr ( "Load Transformation" ),
                       QDir::homePath(),
                       tr ( "Transformation data (*.trans)" ) );

  load_transformation(filename.toAscii().constData());
}

void ProteinModel::on_protein_grouping_comboBox_currentIndexChanged(int ind)
{
  m_protein_grouping_model->set_grouping_type(ind);
}

void ProteinModel::on_reload_pushButton_clicked(bool)
{
  m_protein_atoms_grouping->update_atom_color_bo();
}

void ProteinModel::on_surface_color_colorpicker_colorChanged
    (const QColor &c)
{
  m_surface_color[0] = c.red();
  m_surface_color[1] = c.green();
  m_surface_color[2] = c.blue();
}

protein_grouping_ui_model_t::protein_grouping_ui_model_t(protein_grouping_t * p, QObject *parent )
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
        protein_grouping_t::atom_color_t col =
            m_protein_grouping->get_group_color(index.row());
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
    protein_grouping_t::atom_color_t col
        (color.redF(),color.greenF(),color.blueF());
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
