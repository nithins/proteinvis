#include <viewer.h>

#include <QKeyEvent>
#include <QMessageBox>

void glviewer_t::add_ren(glutils::renderable_ptr_t ren)
{
  m_rens.push_back(renderable_rd_t(ren));
}

void glviewer_t::draw()
{
  float pos[4];

  for(int i = 0; i < m_lights.size();++i)
  {
    GLenum lightnum = i + GL_LIGHT0;

    glGetLightfv(lightnum,GL_POSITION,pos);

    m_lights[i]->getPosition(pos[0], pos[1], pos[2]);

    glLightfv(lightnum, GL_POSITION, pos);
  }


  for(uint i = 0 ; i <m_rens.size();++i)
  {
    glPushMatrix();

    glMultMatrixd(m_rens[i].m_frame->matrix());

    m_rens[i].m_ren->render();

    if(axisIsDrawn()) drawAxis();

    glPopMatrix();
  }

  for(int i = 0; i < m_lights.size();++i)
    drawLight(i+GL_LIGHT0);
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

    drawLight(GL_LIGHT0+i);

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
      std::cout<<"selected light\n";
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

  for(uint i = 0 ; i <m_rens.size();++i)
  {
    m_rens[i].m_ren->gl_init();
  }

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
  const GLfloat light_ambient[4]  = {0.8f, 0.2f, 0.2f, 1.0};
  const GLfloat light_diffuse[4]  = {1.0, 0.4f, 0.4f, 1.0};
  const GLfloat light_specular[4] = {0.8, 0.8, 0.8, 1.0};

  glLightf( GL_LIGHT1, GL_CONSTANT_ATTENUATION, 0.5);
  glLightf( GL_LIGHT1, GL_LINEAR_ATTENUATION, 1.0);
  glLightf( GL_LIGHT1, GL_QUADRATIC_ATTENUATION, 1.5);
  glLightfv(GL_LIGHT1, GL_AMBIENT,  light_ambient);
  glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
  glLightfv(GL_LIGHT1, GL_DIFFUSE,  light_diffuse);

  m_lights[1].reset(new qglviewer::ManipulatedFrame());
  m_lights[1]->setPosition(0.5, 0.5, 0);
  m_lights[1]->setOrientation(qglviewer::Quaternion(qglviewer::Vec(0,0,1), -m_lights[1]->position()));

  glEnable ( GL_COLOR_MATERIAL );

  glColorMaterial ( GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);

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

#include <proteinModel.h>

viewer_mainwindow::viewer_mainwindow()
{
  setupUi(this);

  connect(glviewer,SIGNAL(viewerInitialized()),this,SLOT(init_ui()));
}

void viewer_mainwindow::init_ui()
{
  light_number_spinBox->setMaximum(8);

  light_number_spinBox->setMinimum(0);

  update_light_ui_items(light_number_spinBox->value());

  for(uint i = 0 ; i < m_models.size();++i)
  {
    m_models[i]->init_ui();
  }
}

viewer_mainwindow::~viewer_mainwindow()
{
}

void viewer_mainwindow::add_model(protein_model_ptr_t model)
{
  glviewer->add_ren(model);

  protein_model_ui_t * sa = new protein_model_ui_t(model,this);

  configuration_tabWidget->addTab(sa,model->name().c_str());

  m_models.push_back(sa);
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

  glGetLightfv(GL_LIGHT0+lightno,GL_AMBIENT,c);
  ambient_light_button->setCurrentColor(QColor::fromRgbF(c[0],c[1],c[2],c[3]));

  glGetLightfv(GL_LIGHT0+lightno,GL_DIFFUSE,c);
  diffuse_light_button->setCurrentColor(QColor::fromRgbF(c[0],c[1],c[2],c[3]));

  glGetLightfv(GL_LIGHT0+lightno,GL_SPECULAR,c);
  specular_light_button->setCurrentColor(QColor::fromRgbF(c[0],c[1],c[2],c[3]));
}

void viewer_mainwindow::on_light_number_spinBox_valueChanged ( int i)
{
  update_light_ui_items(i);
}

void viewer_mainwindow::on_ambient_light_button_colorChanged(QColor c)
{
  GLfloat col[] = {c.redF(),c.greenF(),c.blueF(),c.alphaF()};

  glLightfv(GL_LIGHT0+light_number_spinBox->value(),GL_AMBIENT,col);
}
void viewer_mainwindow::on_diffuse_light_button_colorChanged(QColor c)
{
  GLfloat col[] = {c.redF(),c.greenF(),c.blueF(),c.alphaF()};

  glLightfv(GL_LIGHT0+light_number_spinBox->value(),GL_DIFFUSE,col);

}
void viewer_mainwindow::on_specular_light_button_colorChanged(QColor c)
{
  GLfloat col[] = {c.redF(),c.greenF(),c.blueF(),c.alphaF()};

  glLightfv(GL_LIGHT0+light_number_spinBox->value(),GL_SPECULAR,col);
}

void viewer_mainwindow::on_light_type_comboBox_currentIndexChanged(int light_type)
{

  GLenum light_num = GL_LIGHT0+light_number_spinBox->value();
  GLfloat pos[4];

  if(light_type == LIGHTTYPE_NONE)
  {
    glDisable(light_num);
    return;
  }

  glEnable(light_num);

  glGetLightfv(light_num,GL_POSITION,pos);

  if(light_type == LIGHTTYPE_POSITIONAL && pos[3] == 0.0)
    pos[3] = 1.0;

  if(light_type == LIGHTTYPE_DIRECTIONAL && pos[3] != 0.0)
    pos[3] = 0.0;

  glLightfv(light_num,GL_POSITION,pos);
}


