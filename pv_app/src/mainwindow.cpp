#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/typeof/typeof.hpp>

#include <QFileDialog>
#include <QMessageBox>

#include <logutil.h>

#include <protein.h>
#include <proteinModel.h>
#include <pocketModel.h>

#include <mainwindow.h>

viewer_mainwindow::viewer_mainwindow()
{
  setupUi(this);

  protein_grouping_comboBox->blockSignals(true);

  for(int i = 0;i< protein_grouping_t::GROUP_ATOMS_COUNT; ++i)
  {
    std::string grp_name = protein_grouping_t::get_groupby_displaystr((protein_grouping_t::eGroupAtomsBy)i);

    QString str(grp_name.c_str());

    protein_grouping_comboBox->addItem(str);
  }

  protein_grouping_comboBox->blockSignals(false);

  lightseditor->blockSignals(true);

  lightseditor->set_num_lights(glviewer->get_num_lights());

  for(int i = 0; i < glviewer->m_lights.size();++i)
    lightseditor->set_light(i,glviewer->get_light(i));

  atombond_materialEditor->ambient_colorpicker->setEnabled(false);

  atombond_materialEditor->diffuse_colorpicker->setEnabled(false);

  atombond_materialEditor->emission_colorpicker->setEnabled(false);

  lightseditor->blockSignals(false);

  update_model_ui();
}

viewer_mainwindow::~viewer_mainwindow()
{

}

void viewer_mainwindow::on_render_light_graphic_checkBox_clicked ( bool s)
{
  glviewer->m_show_light_graphics = s;

  glviewer->updateGL();
}

void viewer_mainwindow::on_lightseditor_lightChanged
    (int i,const glutils::light_properties_t & il)
{
  glviewer->set_light(i,il);
}

void viewer_mainwindow::on_actionOpen_Protein_triggered(bool)
{
  QString fname = QFileDialog::getOpenFileName
      (this,tr("Select Protein file"),QDir::currentPath(),"Molecules (*.crd *.pdb)");

  if (fname.size() == 0) return;

  protein_model_data_t pmd;

  pmd.model.reset(new protein_model_t(fname.toStdString()));

  pmd.grouping_model.reset
      (new protein_grouping_ui_model_t(pmd.model->m_protein_atoms_grouping,this));

  m_model_data.push_back(pmd);

  glviewer->add_ren(pmd.model);

  current_protein_comboBox->addItem(pmd.model->name().c_str());
}

void viewer_mainwindow::on_actionClose_Protein_triggered(bool)
{
  int i = current_protein_comboBox->currentIndex();

  protein_model_data_t pmd = m_model_data[i];

  m_model_data.erase(m_model_data.begin()+i);

  glviewer->remove_ren(pmd.model);

  current_protein_comboBox->removeItem(i);
}

void viewer_mainwindow::on_actionLoad_Surface_triggered(bool)
{
  QString fname = QFileDialog::getOpenFileName
      (this,tr("Select Surface file"),QDir::currentPath(),"Surface files (*.off)");

  if (fname.size() == 0) return;

  try
  {
    get_active_model()->load_surface(fname.toStdString());
  }
  catch (std::runtime_error e)
  {
    QString line1(tr("failed to load surface file :: "));
    QString line2(tr("\nMessage :: "));

    line1 += fname;
    line2 += e.what();

    QString line3(tr("\nmake sure you are loading a file in off format"));

    QMessageBox::warning(this,tr("failed to load alpha complex"),line1+line2+line3);
  }

  update_model_ui();
}

void viewer_mainwindow::on_actionLoad_Alpha_Shape_triggered(bool)
{
  QString fname = QFileDialog::getOpenFileName
      (this,tr("Select Alpha Shape file"),QDir::currentPath(),"Surface files (*.alp)");

  if (fname.size() == 0) return;

  try
  {
    get_active_model()->load_alpha_shape(fname.toStdString());
  }
  catch (std::runtime_error e)
  {
    QString line1(tr("failed to load alpha complex file :: "));
    QString line2(tr("\nMessage :: "));

    line1 += fname;
    line2 += e.what();

    QString line3(tr("\nmake sure you are loading the correct alpha shape"\
                     "for the current active protein"));

    QMessageBox::warning(this,tr("failed to load alpha complex"),line1+line2+line3);
  }

  update_model_ui();
}
void viewer_mainwindow::on_actionLoad_Pockets_triggered(bool)
{

  QString tepoc_fname = QFileDialog::getOpenFileName
      (this,tr("Select pockets file"),QDir::currentPath(),"Surface files (*.tepoc.all)");

  if (tepoc_fname.size() == 0) return;

  QString tet_fname = QFileDialog::getOpenFileName
      (this,tr("Select Tet file"),QDir::currentPath(),"Surface files (*.tet)");

  if (tet_fname.size() == 0) return;

  try
  {
    get_active_model()->load_pockets(tepoc_fname.toStdString(),tet_fname.toStdString());
  }
  catch (std::runtime_error e)
  {
    QString line1(tr("failed to pocket files:: "));
    QString line2(tr("\nMessage :: "));

    line1 += tet_fname + "," + tepoc_fname;
    line2 += e.what();

    QString line3(tr("\nmake sure you are loading the correct alpha shape"\
                     "for the current active protein"));

    QMessageBox::warning(this,tr("failed to load alpha complex"),line1+line2+line3);
  }

  update_model_ui();

}

viewer_mainwindow::protein_model_data_t viewer_mainwindow::get_active_model_data()
{
  int i = current_protein_comboBox->currentIndex();

  if (i!= -1 )
    return m_model_data[i];
  else
    return protein_model_data_t();
}

protein_model_ptr_t viewer_mainwindow::get_active_model()
{
  return get_active_model_data().model;
}

void viewer_mainwindow::update_model_ui()
{
  protein_model_ptr_t pm  = get_active_model();

  actionLoad_Alpha_Shape->setEnabled(pm);
  actionLoad_Pockets->setEnabled(pm);
  actionLoad_Surface->setEnabled(pm);
  actionClose_Protein->setEnabled(pm);

  alpha_pocket_tab->setEnabled
      (pm && (pm->m_alpha_complex_model||pm->m_pocket_model));

  update_model_atombond_ui(pm);
  update_model_surface_ui(pm);
  update_model_alpha_ui(pm);
  update_model_pocket_ui(pm);
}

void viewer_mainwindow::update_model_atombond_ui(protein_model_ptr_t pm)
{
  atombond_tab->setEnabled(pm);

  if(!pm ) return;

  switch(pm->m_render_model)
  {
  case protein_model_t::RMDL_SPACE_FILL:
    sf_model_radiobutton->setChecked(true); break;
  case protein_model_t::RMDL_BALL_STICK:
    bs_model_radioButton->setChecked(true); break;
  case protein_model_t::RMDL_SMALL_SPACE_FILL:
    ssf_model_radiobutton->setChecked(true); break;
  };

  switch(pm->m_render_mode)
  {
  case protein_model_t::RMDE_FULL:
    full_mol_rendermode_radioButton->setChecked(true); break;
  case protein_model_t::RMDE_BACKBONE:
    bb_mol_rendermode_radioButton->setChecked(true); break;
  case protein_model_t::RMDE_POCKET_ATOMS:
    pocatoms_rendermode_radioButton->setChecked(true); break;
  case protein_model_t::RMDE_NONE:
    none_rendermode_radioButton->setChecked(true); break;
  };

  pocatoms_rendermode_radioButton->setEnabled(pm->m_pocket_model);

  protein_grouping_tableView->setModel(get_active_model_data().grouping_model.get());

  protein_grouping_comboBox->setCurrentIndex(get_active_model_data().grouping_model->get_grouping_type());

  atombond_materialEditor->setMaterial(pm->m_atombond_material);
}

void viewer_mainwindow::update_model_surface_ui(protein_model_ptr_t pm)
{
  surface_tab->setEnabled(pm && pm->m_surface_renderer);

  if(!(pm && pm->m_surface_renderer)) return;

  switch(pm->m_surface_render_mode)
  {
  case protein_model_t::SRM_SHOW:
    show_surface_radioButton->setEnabled(true);break;
  case protein_model_t::SRM_HIDE:
    hide_surface_radioButton->setEnabled(true);break;
  case protein_model_t::SRM_SHOW_WIREFRAME:
    show_surface_wireframe_radioButton->setEnabled(true);break;
  };

  surface_materialEditor->setMaterial(pm->m_surface_material);
}

void viewer_mainwindow::update_model_alpha_ui(protein_model_ptr_t pm)
{
  alpha_complex_groupBox->setEnabled(pm && pm->m_alpha_complex_model);

  if(!pm ||!pm->m_alpha_complex_model ) return;

  show_alpha_complex_tets_checkBox->setChecked
      (pm->m_alpha_complex_render_mode&protein_model_t::ACRM_TETRAHEDRONS);

  show_alpha_complex_tris_checkBox->setChecked
      (pm->m_alpha_complex_render_mode&protein_model_t::ACRM_TRIANGLES);

  show_alpha_complex_edges_checkBox->setChecked
      (pm->m_alpha_complex_render_mode&protein_model_t::ACRM_EDGES);

}

void viewer_mainwindow::update_model_pocket_ui(protein_model_ptr_t pm)
{
  pocket_groupBox->setEnabled(pm && pm->m_pocket_model);

  if(!pm ||!pm->m_pocket_model ) return;

  show_all_pockets_radioButton->setChecked
      (pm->m_pocket_render_mode == protein_model_t::PRM_SHOW_ALL);

  pocket_alpha_num_spinBox->blockSignals(true);
  pocket_alpha_num_spinBox->setMinimum ( 1 );
  pocket_alpha_num_spinBox->setMaximum( pm->m_pocket_model->get_num_alpha() );
  pocket_alpha_num_spinBox->setValue(pm->m_pocket_alpha_num + 1);
  pocket_alpha_num_spinBox->blockSignals(false);

  pocket_num_spinBox->blockSignals(true);
  pocket_num_spinBox->setMinimum ( 1 );
  pocket_num_spinBox->setMaximum( pm->m_pocket_model->get_num_pockets
                                  ( pm->m_pocket_alpha_num ) );
  pocket_num_spinBox->setValue(pm->m_pocket_num +1);

  pocket_num_spinBox->setEnabled
      (pm->m_pocket_render_mode == protein_model_t::PRM_SHOW_ONE);
  pocket_num_label->setEnabled
      (pm->m_pocket_render_mode == protein_model_t::PRM_SHOW_ONE);
  pocket_num_spinBox->blockSignals(false);


  add_radius_sf_model_doubleSpinBox->blockSignals(true);
  add_radius_sf_model_doubleSpinBox->setValue ( pm->m_add_atom_radius );
  add_radius_sf_model_doubleSpinBox->blockSignals(false);

  alpha_value_sf_model_doubleSpinBox->blockSignals(true);
  alpha_value_sf_model_doubleSpinBox->setValue ( pm->m_alpha_value );
  alpha_value_sf_model_doubleSpinBox->blockSignals(false);
}

void viewer_mainwindow::on_current_protein_comboBox_currentIndexChanged (int )
{
  update_model_ui();
}

void viewer_mainwindow::on_sf_model_radiobutton_clicked ( bool checked )
{
  if ( checked )
  {
    get_active_model()->m_render_model = protein_model_t::RMDL_SPACE_FILL;
    glviewer->updateGL();
  }
}

void viewer_mainwindow::on_ssf_model_radiobutton_clicked ( bool checked )
{
  if ( checked )
  {
    get_active_model()->m_render_model = protein_model_t::RMDL_SMALL_SPACE_FILL;
    glviewer->updateGL();
  }
}

void viewer_mainwindow::on_bs_model_radioButton_clicked ( bool checked )
{
  if ( checked )
  {
    get_active_model()->m_render_model = protein_model_t::RMDL_BALL_STICK;
    glviewer->updateGL();
  }
}

void viewer_mainwindow::on_full_mol_rendermode_radioButton_clicked ( bool checked )
{
  if ( checked )
  {
    get_active_model()->m_render_mode = protein_model_t::RMDE_FULL;
    glviewer->updateGL();
  }
}

void viewer_mainwindow::on_bb_mol_rendermode_radioButton_clicked ( bool checked )
{
  if ( checked )
  {
    get_active_model()->m_render_mode = protein_model_t::RMDE_BACKBONE;
    glviewer->updateGL();
  }
}

void viewer_mainwindow::on_pocatoms_rendermode_radioButton_clicked ( bool checked)
{
  if ( checked )
  {
    get_active_model()->m_render_mode = protein_model_t::RMDE_POCKET_ATOMS;
    glviewer->updateGL();
  }
}

void viewer_mainwindow::on_none_rendermode_radioButton_clicked ( bool checked )
{
  if ( checked )
  {
    get_active_model()->m_render_mode = protein_model_t::RMDE_NONE;
    glviewer->updateGL();
  }
}

void viewer_mainwindow::on_protein_grouping_comboBox_currentIndexChanged(int ind)
{
  get_active_model_data().grouping_model->set_grouping_type(ind);
}

void viewer_mainwindow::on_reload_pushButton_clicked(bool)
{
  get_active_model()->m_protein_atoms_grouping->update_atom_color_bo();
  glviewer->updateGL();
}

void viewer_mainwindow::on_atombond_materialEditor_materialChanged
    (const glutils::material_properties_t &m)
{
  get_active_model()->m_atombond_material = m;
  glviewer->updateGL();
}

void viewer_mainwindow::on_show_surface_radioButton_clicked ( bool checked )
{
  if ( checked )
  {
    get_active_model()->m_surface_render_mode = protein_model_t::SRM_SHOW;
    glviewer->updateGL();
  }
}

void viewer_mainwindow::on_show_surface_wireframe_radioButton_clicked ( bool checked )
{
  if ( checked )
  {
    get_active_model()->m_surface_render_mode = protein_model_t::SRM_SHOW_WIREFRAME;
    glviewer->updateGL();
  }
}

void viewer_mainwindow::on_hide_surface_radioButton_clicked ( bool checked )
{
  if ( checked )
  {
    get_active_model()->m_surface_render_mode = protein_model_t::SRM_HIDE;
    glviewer->updateGL();
  }
}

void viewer_mainwindow::on_surface_materialEditor_materialChanged
    (const glutils::material_properties_t &m)
{
  get_active_model()->m_surface_material = m;
  glviewer->updateGL();
}

void viewer_mainwindow::on_show_alpha_complex_tets_checkBox_clicked ( bool c)
{
  get_active_model()->m_alpha_complex_render_mode ^= protein_model_t::ACRM_TETRAHEDRONS;
  get_active_model()->m_alpha_complex_render_mode &= protein_model_t::ACRM_SHOWALL;

  glviewer->updateGL();
}

void viewer_mainwindow::on_show_alpha_complex_tris_checkBox_clicked ( bool c)
{
  get_active_model()->m_alpha_complex_render_mode ^= protein_model_t::ACRM_TRIANGLES;
  get_active_model()->m_alpha_complex_render_mode &= protein_model_t::ACRM_SHOWALL;

  glviewer->updateGL();
}

void viewer_mainwindow::on_show_alpha_complex_edges_checkBox_clicked ( bool c)
{
  get_active_model()->m_alpha_complex_render_mode ^= protein_model_t::ACRM_EDGES;
  get_active_model()->m_alpha_complex_render_mode &= protein_model_t::ACRM_SHOWALL;

  glviewer->updateGL();
}

void viewer_mainwindow::updated_pocket_ui()
{
  if ( !pocket_groupBox->isChecked() )
  {
    get_active_model()->m_pocket_render_mode  = protein_model_t::PRM_SHOW_NONE;
  }
  else
  {
    get_active_model()->m_pocket_alpha_num = pocket_alpha_num_spinBox->value() - 1;

    if ( !show_all_pockets_radioButton->isChecked() )
    {
      get_active_model()->m_pocket_num          = pocket_num_spinBox->value() - 1;

      get_active_model()->m_pocket_render_mode  = protein_model_t::PRM_SHOW_ONE;
    }
    else
    {
      get_active_model()->m_pocket_num          = ( uint ) - 1;

      get_active_model()->m_pocket_render_mode  = protein_model_t::PRM_SHOW_ALL;
    }
  }

  get_active_model()->update_pocket_render_state();

  glviewer->updateGL();

  update_model_pocket_ui(get_active_model());
}

void viewer_mainwindow::on_pocket_groupBox_clicked ( bool )
{
  updated_pocket_ui();
}

void viewer_mainwindow::on_show_all_pockets_radioButton_clicked ( bool )
{
  updated_pocket_ui();
}

void viewer_mainwindow::on_pocket_alpha_num_spinBox_valueChanged ( int )
{
  updated_pocket_ui();
}

void viewer_mainwindow::on_pocket_num_spinBox_valueChanged ( int )
{
  updated_pocket_ui();
}

void viewer_mainwindow::on_add_radius_sf_model_doubleSpinBox_valueChanged ( double v )
{
  get_active_model()->m_add_atom_radius = v;
  glviewer->updateGL();
}

void viewer_mainwindow::on_alpha_value_sf_model_doubleSpinBox_valueChanged ( double v )
{
  get_active_model()->m_alpha_value = v;
  glviewer->updateGL();
}
