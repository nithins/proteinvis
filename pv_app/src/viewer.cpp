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
