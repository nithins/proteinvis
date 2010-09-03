#ifndef VIEWER_H_INCLUDED
#define VIEWER_H_INCLUDED

#include <vector>

#include <QGLViewer/qglviewer.h>

#include <glutils.h>

class glviewer_t : public QGLViewer
{
public:

  typedef boost::shared_ptr<qglviewer::ManipulatedFrame>
      manipulated_frame_ptr_t;

  bool m_is_recording;
  bool m_bf_cull;
  bool m_wireframe;

  bool m_init_state;

  struct renderable_rd_t
  {
    glutils::renderable_ptr_t m_ren;
    manipulated_frame_ptr_t   m_frame;

    renderable_rd_t(glutils::renderable_ptr_t ren):
        m_ren(ren),
        m_frame(new qglviewer::ManipulatedFrame)
    {}
  };

  std::vector<renderable_rd_t>           m_rens;
  std::vector<manipulated_frame_ptr_t>   m_lights;

public:

  glviewer_t(QWidget *par);
  ~glviewer_t();

  void add_ren(glutils::renderable_ptr_t model);
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

#include <ui_viewer.h>

class protein_model_ui_t;

class protein_model_t;

typedef boost::shared_ptr<protein_model_t> protein_model_ptr_t;

class viewer_mainwindow:
    public QMainWindow,
    public Ui::viewer_MainWindow
{
  Q_OBJECT;

  std::vector<protein_model_ui_t *> m_models;

public:
  viewer_mainwindow();
  ~viewer_mainwindow();

  void add_model(protein_model_ptr_t model);

  void update_light_ui_items(int lightno);

private slots:
  void init_ui();
  void on_light_number_spinBox_valueChanged ( int );

  void on_ambient_light_button_colorChanged(QColor);
  void on_diffuse_light_button_colorChanged(QColor);
  void on_specular_light_button_colorChanged(QColor);

  void on_light_type_comboBox_currentIndexChanged(int);

};



#endif
