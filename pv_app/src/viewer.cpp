#include <viewer.h>

#include <QKeyEvent>
#include <QMessageBox>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

void glviewer_t::add_ren(glutils::renderable_ptr_t ren)
{
  if(m_init_state == false)
    throw std::runtime_error("gl not yet inited");

  m_rens.push_back(renderable_rd_t(ren));

  updateGL();
}

void glviewer_t::remove_ren(glutils::renderable_ptr_t ren)
{
  using namespace boost::lambda;

  typedef typeof(m_rens) m_ren_t;

  m_ren_t::iterator i =
      std::find_if(m_rens.begin(),m_rens.end(),bind(&m_ren_t::value_type::m_ren,_1) == ren);

  if(i!= m_rens.end())
    m_rens.erase(i);

  updateGL();
}

glutils::material_properties_t g_light_off_material =
{
  {0.2,0.2,0.2,1},// ambient
  {0.5,0.5,0.5,1},// diffuse
  {0,0,0,1},      // specular
  {0,0,0,1},      // emission
  1             // shininess
};

glutils::material_properties_t g_light_on_material =
{
  {0.2,0.2,0.2,1},// ambient
  {0.5,0.5,0.5,1},// diffuse
  {0,0,0,1},      // specular
  {0.75,0.75,0,1},  // emission
  1             // shininess
};



void draw_bulb(double scaleby = 1.0,bool bulbon = false)
{
  glPushMatrix();

//  glRotatef(90,1,0,0);

  glTranslatef(0,0,-scaleby/2);

  GLUquadric * q = gluNewQuadric();

  glColor3f(0.5,0.5,0.5);

  gluCylinder(q,0.1*scaleby,0.125*scaleby,scaleby,10,10);

  glTranslatef(0,0,scaleby + 0.2*scaleby);

  glPushAttrib ( GL_POLYGON_BIT | GL_ENABLE_BIT );

  glutils::material_properties_t material_front,material_back;

  material_front.read_all(GL_FRONT);

  material_back.read_all(GL_BACK);

  glDisable( GL_COLOR_MATERIAL );

  if(bulbon)
    g_light_on_material.render_all(GL_FRONT_AND_BACK);
  else
    g_light_off_material.render_all(GL_FRONT_AND_BACK);

  gluSphere(q,0.2*scaleby,10,10);

  material_front.render_all(GL_FRONT);

  material_back.render_all(GL_BACK);

  glPopAttrib();

  gluDeleteQuadric(q);

  glPopMatrix();
}

void glviewer_t::draw()
{
  float pos[4];

  for(int i = 0; i < m_lights.size();++i)
  {
    glPushMatrix();

    GLenum lightnum = i + GL_LIGHT0;

    glMultMatrixd(m_lights[i]->matrix());

    glGetLightfv(lightnum,GL_POSITION,pos);

    pos[0] = 0;
    pos[1] = 0;
    pos[2] = (pos[3]!= 0)?(0):(1);

    glLightfv(lightnum, GL_POSITION, pos);

    draw_bulb(0.1,glIsEnabled(lightnum));

    glPopMatrix();
  }

  for(uint i = 0 ; i <m_rens.size();++i)
  {
    glPushMatrix();

    glMultMatrixd(m_rens[i].m_frame->matrix());

    m_rens[i].m_ren->render();

    if(axisIsDrawn()) drawAxis();

    glPopMatrix();
  }

}

void glviewer_t::beginSelection(const QPoint& point)
{
  makeCurrent();

  GLint viewport[4];

  glSelectBuffer ( selectBufferSize(), selectBuffer() );

  glRenderMode ( GL_SELECT );

  glMatrixMode ( GL_PROJECTION );

  glLoadIdentity();

  glGetIntegerv ( GL_VIEWPORT, viewport );

  gluPickMatrix ( point.x(), viewport[3] - point.y(),
                  selectRegionWidth(), selectRegionHeight(), viewport );

  camera()->loadProjectionMatrix(false);

  glMatrixMode ( GL_MODELVIEW );

  glLoadIdentity();

  camera()->loadModelViewMatrix();

  glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT );
}

void glviewer_t::drawWithNames()
{
  GLfloat pos[4];

  for(uint i = 0 ; i <m_rens.size();++i)
  {
    glPushMatrix();

    glMultMatrixd(m_rens[i].m_frame->matrix());

    glPushName(i);

    m_rens[i].m_ren->render();

    glPopName();

    glPopMatrix();
  }

  for(uint i = 0 ; i <m_lights.size();++i)
  {
    glPushName(m_rens.size() + i);

    glPushMatrix();

    glMultMatrixd(m_lights[i]->matrix());

    draw_bulb(0.1);

    glPopMatrix();

    glPopName();
  }
}

void glviewer_t::postSelection(const QPoint& point)
{

  if (selectedName() == -1)
    setManipulatedFrame(NULL);
  else
  {
    if(selectedName() < m_rens.size())
      setManipulatedFrame(m_rens[selectedName()].m_frame.get());
    else
    {
      setManipulatedFrame(m_lights[selectedName() - m_rens.size()].get());
      std::cout<<"selected light"<<std::endl;
    }
  }
}

void glviewer_t::init()
{
  if(m_init_state == true)
    throw std::runtime_error("attempted to re init QGLViewer ..rascal\n");

  m_init_state = true;

  glutils::init();

  setSnapshotFormat("PNG");

  setSnapshotQuality(100);

  glEnable ( GL_CULL_FACE );

  glCullFace ( GL_BACK );

  glPolygonMode ( GL_FRONT, GL_FILL );

  glPolygonMode ( GL_BACK, GL_LINE );

  setBackgroundColor(Qt::white);

  restoreStateFromFile();

  setHandlerKeyboardModifiers(QGLViewer::CAMERA, Qt::AltModifier);

  setHandlerKeyboardModifiers(QGLViewer::FRAME,  Qt::NoModifier);

  setHandlerKeyboardModifiers(QGLViewer::CAMERA, Qt::ControlModifier);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  m_lights.resize(2);

  // Light0 is a classical directionnal light
  glEnable(GL_LIGHT0);
  const GLfloat light_ambient0[4]  = {0.2f, 0.2f, 0.2f, 1.0f};
  const GLfloat light_diffuse0[4]  = {0.8f, 0.8f, 0.8f, 1.0f};
  const GLfloat light_specular0[4] = {0.8f, 0.8f, 0.8f, 1.0f};

  glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient0);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);
  glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse0);

  m_lights[0].reset(new qglviewer::ManipulatedFrame());
  m_lights[0]->setPosition(-0.5, 0.5, 0);

  // Light1 is a point light
  glEnable(GL_LIGHT1);
  const GLfloat light_ambient[4]  = {0.2f, 0.2f, 0.2f, 1.0};
  const GLfloat light_diffuse[4]  = {0.4f, 0.4f, 0.4f, 1.0};
  const GLfloat light_specular[4] = {0.8f, 0.8f, 0.8f, 1.0};

  glLightf( GL_LIGHT1, GL_CONSTANT_ATTENUATION, 0.5);
  glLightf( GL_LIGHT1, GL_LINEAR_ATTENUATION, 1.0);
  glLightf( GL_LIGHT1, GL_QUADRATIC_ATTENUATION, 1.5);
  glLightfv(GL_LIGHT1, GL_AMBIENT,  light_ambient);
  glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
  glLightfv(GL_LIGHT1, GL_DIFFUSE,  light_diffuse);

  m_lights[1].reset(new qglviewer::ManipulatedFrame());
  m_lights[1]->setPosition(0.5, 0.5, 0);
  m_lights[1]->setOrientation(qglviewer::Quaternion(qglviewer::Vec(0,0,1), -m_lights[1]->position()));

  QGLViewer::init();
}

glviewer_t::glviewer_t(QWidget * par):
    m_is_recording(false),
    m_bf_cull(true),
    m_wireframe(false),
    m_init_state(false)
{
  setParent(par);
}

glviewer_t::~glviewer_t()
{

}

void glviewer_t::keyPressEvent(QKeyEvent *e)
{
  const Qt::KeyboardModifiers modifiers = e->modifiers();

  if ((e->key()==Qt::Key_C) && (modifiers==Qt::ControlModifier))
  {
    m_is_recording = !m_is_recording;

    if(m_is_recording)
      connect(this, SIGNAL(drawFinished(bool)),this, SLOT(saveSnapshot(bool)));
    else
      disconnect(this, SIGNAL(drawFinished(bool)),this, SLOT(saveSnapshot(bool)));
  }
  else if ((e->key()==Qt::Key_B) && (modifiers==Qt::ControlModifier))
  {
    m_bf_cull = !m_bf_cull;

    if(m_bf_cull)
      glEnable ( GL_CULL_FACE );
    else
      glDisable ( GL_CULL_FACE );
  }
  else if ((e->key()==Qt::Key_W) && (modifiers==Qt::ControlModifier))
  {
    m_wireframe = !m_wireframe;

    if(m_wireframe)
      glPolygonMode ( GL_FRONT_AND_BACK, GL_LINE );
    else
    {
      glPolygonMode ( GL_FRONT, GL_FILL );
      glPolygonMode ( GL_BACK, GL_LINE );
    }
  }
  else
  {
    QGLViewer::keyPressEvent(e);
  }

  updateGL();
}

QString glviewer_t::helpString() const
{
  QString text("<h2>Protein Viewer</h2>");
  return text;
}

#include <protein.h>
#include <proteinModel.h>
#include <pocketModel.h>
#include <QFileDialog>

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

  connect(glviewer,SIGNAL(viewerInitialized()),this,SLOT(init_ui()));

  update_model_ui();
}

void viewer_mainwindow::init_ui()
{
  update_light_ui_items(light_number_spinBox->value());
}

viewer_mainwindow::~viewer_mainwindow()
{

}

enum eLightType
{
  LIGHTTYPE_DIRECTIONAL,
  LIGHTTYPE_POSITIONAL,
  LIGHTTYPE_NONE
};

void viewer_mainwindow::update_light_ui_items(int lightno)
{
  if(glIsEnabled(GL_LIGHT0+lightno) == false)
  {
    light_type_comboBox->setCurrentIndex(LIGHTTYPE_NONE);
    return;
  }

  GLfloat pos[4];

  glGetLightfv(GL_LIGHT0+lightno,GL_POSITION,pos);

  if(pos[3] == 0)
    light_type_comboBox->setCurrentIndex(LIGHTTYPE_DIRECTIONAL);
  else
    light_type_comboBox->setCurrentIndex(LIGHTTYPE_POSITIONAL);

  GLfloat c[4];

  ambient_light_button->blockSignals(true);
  diffuse_light_button->blockSignals(true);
  specular_light_button->blockSignals(true);

  glGetLightfv(GL_LIGHT0+lightno,GL_AMBIENT,c);
  ambient_light_button->setCurrentColor(QColor::fromRgbF(c[0],c[1],c[2],c[3]));

  glGetLightfv(GL_LIGHT0+lightno,GL_DIFFUSE,c);
  diffuse_light_button->setCurrentColor(QColor::fromRgbF(c[0],c[1],c[2],c[3]));

  glGetLightfv(GL_LIGHT0+lightno,GL_SPECULAR,c);
  specular_light_button->setCurrentColor(QColor::fromRgbF(c[0],c[1],c[2],c[3]));

  ambient_light_button->blockSignals(false);
  diffuse_light_button->blockSignals(false);
  specular_light_button->blockSignals(false);

}

void viewer_mainwindow::on_light_number_spinBox_valueChanged ( int i)
{
  update_light_ui_items(i);
}

void viewer_mainwindow::on_ambient_light_button_colorChanged(QColor c)
{
  GLfloat col[] = {c.redF(),c.greenF(),c.blueF(),c.alphaF()};

  glLightfv(GL_LIGHT0+light_number_spinBox->value(),GL_AMBIENT,col);

  glviewer->updateGL();
}
void viewer_mainwindow::on_diffuse_light_button_colorChanged(QColor c)
{
  GLfloat col[] = {c.redF(),c.greenF(),c.blueF(),c.alphaF()};

  glLightfv(GL_LIGHT0+light_number_spinBox->value(),GL_DIFFUSE,col);

  glviewer->updateGL();
}
void viewer_mainwindow::on_specular_light_button_colorChanged(QColor c)
{
  GLfloat col[] = {c.redF(),c.greenF(),c.blueF(),c.alphaF()};

  glLightfv(GL_LIGHT0+light_number_spinBox->value(),GL_SPECULAR,col);

  glviewer->updateGL();
}

void viewer_mainwindow::on_light_type_comboBox_currentIndexChanged(int light_type)
{
  GLenum light_num = GL_LIGHT0+light_number_spinBox->value();

  GLfloat pos[4];

  if(light_type == LIGHTTYPE_NONE)
  {
    glDisable(light_num);

    glviewer->updateGL();

    return;
  }

  glEnable(light_num);

  glGetLightfv(light_num,GL_POSITION,pos);

  if(light_type == LIGHTTYPE_POSITIONAL && pos[3] == 0.0)
    pos[3] = 1.0;

  if(light_type == LIGHTTYPE_DIRECTIONAL && pos[3] != 0.0)
    pos[3] = 0.0;

  glLightfv(light_num,GL_POSITION,pos);

  glviewer->updateGL();
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

  pocket_alpha_num_spinBox->setValue(pm->m_pocket_alpha_num);

  show_all_pockets_radioButton->setChecked
      (pm->m_pocket_render_mode == protein_model_t::PRM_SHOW_ALL);

  pocket_alpha_num_spinBox->setMinimum ( 1 );
  pocket_alpha_num_spinBox->setMaximum( pm->m_pocket_model->get_num_alpha() );

  pocket_num_spinBox->setMinimum ( 1 );
  pocket_num_spinBox->setMaximum( pm->m_pocket_model->get_num_pockets
                                  ( pm->m_pocket_alpha_num ) );

  add_radius_sf_model_doubleSpinBox->setValue ( pm->m_add_atom_radius );

  alpha_value_sf_model_doubleSpinBox->setValue ( pm->m_alpha_value );
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
