#ifndef VIEWER_H_INCLUDED
#define VIEWER_H_INCLUDED

#include <QGLViewer/qglviewer.h>

#include <glutils.h>

#include <aabb.h>

class glviewer_t : public QGLViewer
{
public:

  typedef boost::shared_ptr<qglviewer::ManipulatedFrame>
      manipulated_frame_ptr_t;

  typedef glutils::light_properties_t     light_properties_t;

public:

  bool m_is_recording;
  bool m_bf_cull;
  bool m_wireframe;
  bool m_show_light_graphics;

  bool m_init_state;

  struct renderable_rd_t
  {
    glutils::renderable_ptr_t        ren;
    manipulated_frame_ptr_t          frame;
    aabb::aabb_t<double,3>           extent;
    aabb::aabb_t<double,3>::point_t  center;//derived from m_extent


    renderable_rd_t(glutils::renderable_ptr_t r):
        ren(r),
        frame(new qglviewer::ManipulatedFrame)
    {
      frame->setWheelSensitivity(0.0); // no zooming
    }
  };

  struct light_data_t
  {
    glutils::light_properties_t props;
    manipulated_frame_ptr_t     frame;

    light_data_t(const glutils::light_properties_t &p):
        props(p),
        frame(new qglviewer::ManipulatedFrame)
    {
      frame->setWheelSensitivity(0.0);// no zooming
    }
  };

  std::vector<renderable_rd_t>           m_rens;
  std::vector<light_data_t>              m_lights;

  // derived from m_extent's of the m_rens
  double                                 m_max_extent;

public:

  glviewer_t(QWidget *par);

  void add_ren(glutils::renderable_ptr_t model);

  void remove_ren(glutils::renderable_ptr_t model);

  void update_extents();

  int get_num_lights() const ;

  const light_properties_t& get_light(int n) const;

  void set_light(int n,const light_properties_t & l);

protected:

  virtual void draw();

  virtual void beginSelection(const QPoint& point);
  virtual void drawWithNames();
//  virtual void endSelection(const QPoint& point);
  virtual void postSelection(const QPoint& point);

  virtual void init();
  virtual QString helpString() const;
  virtual void keyPressEvent(QKeyEvent *e);
};
#endif
