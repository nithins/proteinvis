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

  void remove_ren(glutils::renderable_ptr_t model);
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

class protein_grouping_ui_model_t;

typedef boost::shared_ptr<protein_grouping_ui_model_t>
    protein_grouping_ui_model_ptr_t;

class viewer_mainwindow:
    public QMainWindow,
    public Ui::viewer_MainWindow
{
  Q_OBJECT;

  struct protein_model_data_t
  {
    protein_model_ptr_t             model;
    protein_grouping_ui_model_ptr_t grouping_model;
  };

  std::vector<protein_model_data_t> m_model_data;

public:
  viewer_mainwindow();
  ~viewer_mainwindow();


// light ui

private slots:
  void init_ui();
  void on_light_number_spinBox_valueChanged ( int );
  void on_ambient_light_button_colorChanged(QColor);
  void on_diffuse_light_button_colorChanged(QColor);
  void on_specular_light_button_colorChanged(QColor);
  void on_light_type_comboBox_currentIndexChanged(int);

private:

// update_* push changes to ui
  void update_light_ui_items(int lightno);

  void update_model_ui();
  void update_model_atombond_ui(protein_model_ptr_t );
  void update_model_surface_ui(protein_model_ptr_t );
  void update_model_alpha_ui(protein_model_ptr_t );
  void update_model_pocket_ui(protein_model_ptr_t );

// updated_* pull changes from ui
  void updated_pocket_ui();


  protein_model_ptr_t  get_active_model();
  protein_model_data_t get_active_model_data();

private slots:
// protein_model ui
  void on_actionOpen_Protein_triggered(bool);
  void on_actionClose_Protein_triggered(bool);
  void on_actionLoad_Surface_triggered(bool);
  void on_actionLoad_Alpha_Shape_triggered(bool);
  void on_actionLoad_Pockets_triggered(bool);
  void on_current_protein_comboBox_currentIndexChanged (int );


// atombond page ui
  void on_full_mol_rendermode_radioButton_clicked ( bool );
  void on_bb_mol_rendermode_radioButton_clicked ( bool );
  void on_none_rendermode_radioButton_clicked ( bool );

  void on_sf_model_radiobutton_clicked ( bool );
  void on_ssf_model_radiobutton_clicked ( bool );
  void on_bs_model_radioButton_clicked ( bool );

  void on_protein_grouping_comboBox_currentIndexChanged(int);
  void on_reload_pushButton_clicked(bool);
  void on_atombond_materialEditor_materialChanged
      (const glutils::material_properties_t &m);

// surface page ui
  void on_show_surface_radioButton_clicked ( bool );
  void on_show_surface_wireframe_radioButton_clicked ( bool );
  void on_hide_surface_radioButton_clicked ( bool );
  void on_surface_materialEditor_materialChanged
      (const glutils::material_properties_t &m);

// alpha complex ui
  void on_show_alpha_complex_tets_checkBox_clicked ( bool );
  void on_show_alpha_complex_tris_checkBox_clicked ( bool );
  void on_show_alpha_complex_edges_checkBox_clicked ( bool );

// pockets ui
  void on_pocket_groupBox_clicked ( bool );
  void on_show_all_pockets_radioButton_clicked ( bool );
  void on_pocket_alpha_num_spinBox_valueChanged ( int );
  void on_pocket_num_spinBox_valueChanged ( int );

// alpha and pockets ui
  void on_add_radius_sf_model_doubleSpinBox_valueChanged ( double );
  void on_alpha_value_sf_model_doubleSpinBox_valueChanged ( double );
};



#endif
