#ifndef VIEWER_H_INCLUDED
#define VIEWER_H_INCLUDED

#include <vector>

#include <QGLViewer/qglviewer.h>

#include <glutils.h>

class glviewer_t : public QGLViewer
{
private:

  bool m_is_recording;
  bool m_bf_cull;
  bool m_wireframe;

  bool m_init_state;

  struct renderable_item_t
  {
    glutils::renderable_ptr_t     m_ren;
    qglviewer::ManipulatedFrame * m_frame;

    renderable_item_t(glutils::renderable_ptr_t ren):
        m_ren(ren),
        m_frame(new qglviewer::ManipulatedFrame)
    {}
  };

  std::vector<renderable_item_t> m_rens;

  qglviewer::ManipulatedFrame* light1;
  qglviewer::ManipulatedFrame* light2;

public:

  glviewer_t(QWidget *par);
  ~glviewer_t();

  void add_ren(glutils::renderable_ptr_t ren);
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

class viewer_mainwindow:
    public QDialog,
    public Ui::viewer_mainwindow_Dialog
{
  Q_OBJECT;

public:
  viewer_mainwindow();
  ~viewer_mainwindow();

  void add_ren(glutils::renderable_ptr_t ren );
  void add_frame(QFrame *,const std::string &title);
};



#endif
