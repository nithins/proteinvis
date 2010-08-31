#ifndef  __PROTIEN_MODEL_H
#define  __PROTIEN_MODEL_H

#include <iostream>
#include <string>

#include <QFrame>

#include <glutils.h>
#include <cpputils.h>

#include <model.h>
#include <input.h>

class onelevel_model_t;

class protein_t;

class protein_rd_t;

class protein_grouping_t;

class alpha_complex_model_t;

class pocket_model_t;

class IModelController;

namespace Ui
{
  class ProteinModelFrame;
}
class protein_grouping_ui_model_t;

typedef glutils::renderable_t array_renderer_t;

class ProteinModel:
    public QFrame,
    public IModel,
    public IHandleInput,
    public glutils::renderable_t

{
public:
  typedef n_vector_t<unsigned char,3> uc_color_t;

protected:

  onelevel_model_t      *m_onelevel_model;
  protein_t             *m_protein;
  protein_rd_t          *m_protein_rd;
  protein_grouping_t    *m_protein_atoms_grouping;
  array_renderer_t      *m_surface_renderer;
  alpha_complex_model_t *m_alpha_complex_model;
  pocket_model_t        *m_pocket_model;


  // other state stuff
  uint                   m_render_model;
  uint                   m_render_mode;
  uint                   m_surface_render_mode;
  uint                   m_alpha_complex_render_mode;
  uint                   m_pocket_alpha_num;
  uint                   m_pocket_num;
  uint                   m_pocket_render_mode;
  double                 m_add_atom_radius;
  double                 m_alpha_value;
  uc_color_t             m_surface_color;

  std::string            m_protein_filename;
  std::string            m_surface_filename;
  std::string            m_alpha_complex_filename;
  std::string            m_pocket_filename;
  std::string            m_tetra_filename;
  std::string            m_protein_name;

  IModelController      *m_controller;

  void render_onelevel() const;

  void update_pocket_render_state();
  void update_sf_model_for_pocket();

  void save_to_file(std::string filename );
  void save_transformation(std::string filename );
  void load_transformation(std::string filename );


public:
  ProteinModel ( std::string , std::string , std::string , std::string, std::string );
  virtual ~ProteinModel ();

  int render()
  {
    Render();
  }

  void gl_init();

  virtual int  Render() const;
  //   virtual int  RenderForPick() const;
  virtual void Reset();

  virtual std::string Name() const
  {
    return m_protein_name.c_str();
  }

  bool MousePressedEvent
      ( const int &x, const int &y, const eMouseButton &mb,
        const eKeyFlags &,const eMouseFlags &);

  bool MouseReleasedEvent
      ( const int &x, const int &y, const eMouseButton &mb,
        const eKeyFlags &,const eMouseFlags &);

  bool MouseMovedEvent
      ( const int &x, const int &y, const int &, const int &,
        const eKeyFlags &,const eMouseFlags &mf);

  virtual bool WheelEvent
    ( const int &, const int &, const int &d,
      const eKeyFlags &kf,const eMouseFlags &);

private:
  //setup stuff
  void setup_surface();

  // ui stuff

private slots:

  void on_full_mol_rendermode_radioButton_clicked ( bool );
  void on_bb_mol_rendermode_radioButton_clicked ( bool );
  void on_none_rendermode_radioButton_clicked ( bool );
  void on_pocatoms_rendermode_radioButton_clicked ( bool );

  void on_sf_model_radiobutton_clicked ( bool );
  void on_ssf_model_radiobutton_clicked ( bool );
  void on_bs_model_radioButton_clicked ( bool );
  void on_add_radius_sf_model_doubleSpinBox_valueChanged ( double );
  void on_alpha_value_sf_model_doubleSpinBox_valueChanged ( double );

  void on_show_surface_radioButton_clicked ( bool );
  void on_show_surface_wireframe_radioButton_clicked ( bool );
  void on_hide_surface_radioButton_clicked ( bool );

  void on_pocket_groupBox_clicked ( bool );
  void on_show_all_pockets_radioButton_clicked ( bool );
  void on_pocket_alpha_num_spinBox_valueChanged ( int );
  void on_pocket_num_spinBox_valueChanged ( int );

  void on_alpha_complex_groupBox_clicked ( bool );
  void on_show_alpha_complex_tets_checkBox_clicked ( bool );
  void on_show_alpha_complex_tris_checkBox_clicked ( bool );
  void on_show_alpha_complex_edges_checkBox_clicked ( bool );

  void on_save_to_file_pushButton_clicked ( bool );

  void on_protein_grouping_comboBox_currentIndexChanged(int i);
  void on_reload_pushButton_clicked(bool);
  void on_surface_color_colorpicker_colorChanged(const QColor &);

  void on_save_trans_pushButton_clicked(bool);
  void on_load_trans_pushButton_clicked(bool);


private:
  void init_ui();
  void destroy_ui();

  void pocket_ui_state_updated();
  void alpha_complex_ui_state_updated();

  protein_grouping_ui_model_t * m_protein_grouping_model;

  Ui::ProteinModelFrame * m_ui;

  Q_OBJECT;

public:
  virtual QFrame * getQFrame()
  {
    return this;
  }

};


#include <QAbstractListModel>

class protein_grouping_ui_model_t:public QAbstractItemModel
{
  Q_OBJECT

public:

  protein_grouping_ui_model_t(protein_grouping_t * p, QObject *parent = 0);

  ~protein_grouping_ui_model_t();

  QVariant data(const QModelIndex &index, int role) const;

  Qt::ItemFlags flags(const QModelIndex &index) const;

  QVariant headerData
      (int section, Qt::Orientation orientation,
       int role = Qt::DisplayRole) const;

  QModelIndex index
      (int row, int column,
       const QModelIndex &) const;

  QModelIndex parent(const QModelIndex &) const;

  int rowCount(const QModelIndex & ) const;

  int columnCount(const QModelIndex &) const;

  virtual bool setData
      ( const QModelIndex & index,
        const QVariant & value,
        int role );

  void set_grouping_type(int ind);

private:
  protein_grouping_t *m_protein_grouping;
};
#endif//__PROTIEN_MODEL_H
