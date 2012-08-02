#ifndef MAINWINDOW_H_INCLUDED
#define MAINWINDOW_H_INCLUDED

#include <ui_mainwindow.h>
#include <cpputils.h>

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


private:

// update_* push changes to ui
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
// light ui
  void on_render_light_graphic_checkBox_clicked ( bool s);
  void on_lightseditor_lightChanged(int i,const glutils::light_properties_t &);

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
  void on_calpha_rendermode_radioButton_clicked ( bool );
  void on_pocatoms_rendermode_radioButton_clicked ( bool );
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

  void on_sec_tubes_checkbox_clicked(bool checked);
  void on_sec_helices_checkbox_clicked(bool checked);
  void on_sec_sheets_checkbox_clicked(bool checked);

  void on_spinBox_valueChanged(int arg1);
  void on_spinBox_2_valueChanged(int arg1);
};

#endif
