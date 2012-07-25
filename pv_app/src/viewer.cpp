#include <viewer.h>
#include <proteinModel.h>

#include <QKeyEvent>
#include <QMessageBox>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/foreach.hpp>

glutils::material_properties_t g_light_off_material =
{
  {0.2,0.2,0.2,1}, // ambient
  {0.5,0.5,0.5,1}, // diffuse
  {0,0,0,1},       // specular
  {0,0,0,1},       // emission
  1                // shininess
};

glutils::material_properties_t g_light_on_material =
{
  {0.2,0.2,0.2,1}, // ambient
  {0.5,0.5,0.5,1}, // diffuse
  {0,0,0,1},       // specular
  {0.75,0.75,0,1}, // emission
  1                // shininess
};

// two directional lights by default
// as a convention directional lights are directed along 0,0,1 model coords
// and positional lights are located at 0,0,0 model coords
// another frame is associated with the light to move it around
// and orient it in world coords
glutils::light_properties_t g_lights_default[] =
{
  {
    {0.2f, 0.2f, 0.2f, 1.0f}, // ambient
    {0.5f, 0.5f, 0.5f, 1.0f}, // diffuse
    {1.0f, 1.0f, 1.0f, 1.0f}, // specular
    {0.0f, 0.0f, -1.0f, 0.0f}, // position
    1,0,0,                    // attenution c + l + q
    true                      // enabled
  },

  {
    {0.2f, 0.2f, 0.2f, 1.0f}, // ambient
    {0.5f, 0.5f, 0.5f, 1.0f}, // diffuse
    {1.0f, 1.0f, 1.0f, 1.0f}, // specular
    {0.0f, 0.0f, -1.0f, 0.0f}, // position
    1,0,0,                    // attenution c + l + q
    true                      // enabled
  }
};

glviewer_t::glviewer_t(QWidget * par):
    m_is_recording(false),
    m_bf_cull(true),
    m_wireframe(false),
    m_init_state(false),
    m_max_extent(0.0),
    m_show_light_graphics(false)
{
  setParent(par);

  m_lights.push_back(light_data_t(g_lights_default[0]));
  m_lights.push_back(light_data_t(g_lights_default[1]));

  m_lights[0].frame->setPosition(0.5, 0.5, 0);
  m_lights[1].frame->setPosition(-0.5, 0.5, 0);
  m_lights[1].frame->setOrientation
      (qglviewer::Quaternion(qglviewer::Vec(0,0,1),-m_lights[1].frame->position()));

}

void glviewer_t::add_ren(glutils::renderable_ptr_t ren)
{
  if(m_init_state == false)
    throw std::runtime_error("gl not yet inited");

  m_rens.push_back(renderable_rd_t(ren));

  update_extents();

  setManipulatedFrame((*m_rens.rbegin()).frame.get());

  updateGL();
}

void glviewer_t::remove_ren(glutils::renderable_ptr_t ren)
{
  using namespace boost::lambda;

  setManipulatedFrame(NULL);

  BOOST_AUTO(i,std::find_if(m_rens.begin(),m_rens.end(),bind(&renderable_rd_t::ren,_1) == ren));

  if(i!= m_rens.end()) m_rens.erase(i);

  update_extents();

  updateGL();
}

void glviewer_t::update_extents()
{
  m_max_extent = 0;

  BOOST_FOREACH(typeof(*m_rens.begin()) &rd,m_rens)
  {
    if(rd.ren->get_extent(rd.extent[0].data()) == false)
      throw std::runtime_error("get_extent failed");

    rd.center = (rd.extent.upper_corner() + rd.extent.lower_corner())/2;

    BOOST_AUTO(rd_span,rd.extent.span());

    double rd_max_extent = *std::max_element(rd_span.begin(),rd_span.end());

    if(rd_max_extent > m_max_extent)
    {
      m_max_extent = rd_max_extent;
    }
  }
}

int glviewer_t::get_num_lights() const
{
  return m_lights.size();
}

const glviewer_t::light_properties_t& glviewer_t::get_light(int n) const
{
  assert(n<m_lights.size()&&"n > numlights");

  return m_lights[n].props;
}

void glviewer_t::set_light(int n,const light_properties_t & nl)
{
  assert(n<m_lights.size()&&"n > numlights");

  light_properties_t &l=m_lights[n].props;

  l = nl;

  if(l.position[3] == 0)
    l.position = glutils::vertex4f_t(0,0,1,0);
  else
    l.position = glutils::vertex4f_t(0,0,0,1);

  // lights are just blacked out not switched off
  if(l.enabled == false)
  {
    l.ambient  = glutils::color4f_t(0,0,0,1);
    l.diffuse  = glutils::color4f_t(0,0,0,1);
    l.specular = glutils::color4f_t(0,0,0,1);
  }

  updateGL();
}

void draw_bulb(double scaleby = 1.0,bool bulbon = false)
{
  glPushMatrix();

  glTranslatef(0,0,-scaleby/2);

  GLUquadric * q = gluNewQuadric();

  glColor3f(0.5,0.5,0.5);

  gluCylinder(q,0.125*scaleby,0.1*scaleby,scaleby,10,10);

  glTranslatef(0,0,-0.2*scaleby );

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
  for(int i = 0 ; i < m_lights.size();++i)
  {
    const light_data_t & l = m_lights[i];

    glPushMatrix();

    glMultMatrixd(l.frame->matrix());

    l.props.render(i);

    if(m_show_light_graphics)
      draw_bulb(0.1,l.props.enabled);

    glPopMatrix();
  }

  for(uint i = 0 ; i <m_rens.size();++i)
  {
    glPushMatrix();

    glMultMatrixd(m_rens[i].frame->matrix());

    if(axisIsDrawn()) drawAxis();

    glScalef(1.0/m_max_extent,1.0/m_max_extent,1.0/m_max_extent);

    BOOST_AUTO(const &c,m_rens[i].center);

    glTranslated(-c[0],-c[1],-c[2]);

    m_rens[i].ren->render();

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
  for(uint i = 0 ; i <m_rens.size();++i)
  {
    glPushMatrix();

    glMultMatrixd(m_rens[i].frame->matrix());

    glScalef(1.0/m_max_extent,1.0/m_max_extent,1.0/m_max_extent);

    BOOST_AUTO(const &c,m_rens[i].center);

    glTranslated(-c[0],-c[1],-c[2]);

    glPushName(i);

    m_rens[i].ren->render();

    glPopName();

    glPopMatrix();
  }

  if(m_show_light_graphics == true)
  {
    for(uint i = 0 ; i <m_lights.size();++i)
    {
      glPushName(m_rens.size() + i);

      glPushMatrix();

      glMultMatrixd(m_lights[i].frame->matrix());

      draw_bulb(0.1);

      glPopMatrix();

      glPopName();
    }
  }
}

void glviewer_t::postSelection(const QPoint& point)
{
  if (selectedName() == -1)
    setManipulatedFrame(NULL);
  else
  {
    if(selectedName() < m_rens.size())
      setManipulatedFrame(m_rens[selectedName()].frame.get());
    else
      setManipulatedFrame(m_lights[selectedName() - m_rens.size()].frame.get());
  }
}

void glviewer_t::init()
{
  if(m_init_state == true)
    throw std::runtime_error("attempted to re init QGLViewer ..rascal\n");

  m_init_state = true;

  glutils::init();

  protein_model_t::initgl();

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

  setAnimationPeriod(0);

  glLoadIdentity();
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
