#include <viewer.h>

#include <QKeyEvent>
#include <QMessageBox>

void glviewer_t::add_ren(glutils::renderable_ptr_t ren)
{
  m_rens.push_back(renderable_item_t(ren));

  if(m_init_state == true) ren->gl_init();
}

void glviewer_t::draw()
{
  float pos[4] = {1.0, 0.5, 1.0, 0.0};

  glLightfv(GL_LIGHT0, GL_POSITION, pos);

  pos[3] = 1.0;

  light1->getPosition(pos[0], pos[1], pos[2]);

  glLightfv(GL_LIGHT1, GL_POSITION, pos);

  glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION,
            light1->inverseTransformOf(qglviewer::Vec(0,0,1)));

  light2->getPosition(pos[0], pos[1], pos[2]);

  glLightfv(GL_LIGHT2, GL_POSITION, pos);

  for(uint i = 0 ; i <m_rens.size();++i)
  {
    glPushMatrix();

    glMultMatrixd(m_rens[i].m_frame->matrix());

    m_rens[i].m_ren->render();

    if(axisIsDrawn()) drawAxis();

    glPopMatrix();
  }

  drawLight(GL_LIGHT0);

  if (light1->grabsMouse())
    drawLight(GL_LIGHT1, 1.2f);
  else
    drawLight(GL_LIGHT1);

  if (light2->grabsMouse())
    drawLight(GL_LIGHT2, 1.2f);
  else
    drawLight(GL_LIGHT2);

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
}

void glviewer_t::postSelection(const QPoint& point)
{
  if (selectedName() == -1)
    setManipulatedFrame(NULL);
  else
    setManipulatedFrame(m_rens[selectedName()].m_frame);

}

void glviewer_t::init()
{
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

  m_init_state = true;

  setHandlerKeyboardModifiers(QGLViewer::CAMERA, Qt::AltModifier);

  setHandlerKeyboardModifiers(QGLViewer::FRAME,  Qt::NoModifier);

  setHandlerKeyboardModifiers(QGLViewer::CAMERA, Qt::ControlModifier);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // Light0 is the default ambient light
  glEnable(GL_LIGHT0);

  // Light1 is a spot light
  glEnable(GL_LIGHT1);
  const GLfloat light_ambient[4]  = {0.8f, 0.2f, 0.2f, 1.0};
  const GLfloat light_diffuse[4]  = {1.0, 0.4f, 0.4f, 1.0};
  const GLfloat light_specular[4] = {1.0, 0.0, 0.0, 1.0};

  glLightf( GL_LIGHT1, GL_SPOT_EXPONENT,  3.0);
  glLightf( GL_LIGHT1, GL_SPOT_CUTOFF,    20.0);
  glLightf( GL_LIGHT1, GL_CONSTANT_ATTENUATION, 0.5);
  glLightf( GL_LIGHT1, GL_LINEAR_ATTENUATION, 1.0);
  glLightf( GL_LIGHT1, GL_QUADRATIC_ATTENUATION, 1.5);
  glLightfv(GL_LIGHT1, GL_AMBIENT,  light_ambient);
  glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
  glLightfv(GL_LIGHT1, GL_DIFFUSE,  light_diffuse);

  // Light2 is a classical directionnal light
  glEnable(GL_LIGHT2);
  const GLfloat light_ambient2[4]  = {0.2f, 0.2f, 2.0, 1.0};
  const GLfloat light_diffuse2[4]  = {0.8f, 0.8f, 1.0, 1.0};
  const GLfloat light_specular2[4] = {0.0, 0.0, 1.0, 1.0};

  glLightfv(GL_LIGHT2, GL_AMBIENT,  light_ambient2);
  glLightfv(GL_LIGHT2, GL_SPECULAR, light_specular2);
  glLightfv(GL_LIGHT2, GL_DIFFUSE,  light_diffuse2);

  light1 = new qglviewer::ManipulatedFrame();
  light2 = new qglviewer::ManipulatedFrame();
  setMouseTracking(true);

  light1->setPosition(0.5, 0.5, 0);
  // Align z axis with -position direction : look at scene center
  light1->setOrientation(qglviewer::Quaternion(qglviewer::Vec(0,0,1), -light1->position()));

  light2->setPosition(-0.5, 0.5, 0);
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

viewer_mainwindow::viewer_mainwindow()
{
  setupUi(this);
}

viewer_mainwindow::~viewer_mainwindow()
{
  for(int i = configuration_toolBox->count()-1; i >=0 ; --i)
    configuration_toolBox->removeItem(i);
}

void viewer_mainwindow::add_ren(glutils::renderable_ptr_t ren )
{
  glviewer->add_ren(ren);
}

void viewer_mainwindow::add_frame(QFrame *f,const std::string &title)
{
  configuration_toolBox->addItem(f,title.c_str());
}

