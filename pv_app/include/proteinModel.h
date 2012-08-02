#ifndef  __PROTIEN_MODEL_H
#define  __PROTIEN_MODEL_H

#include <iostream>
#include <string>

#include <glutils.h>
#include <cpputils.h>
#include <secondaryModel.h>

class protein_t;

class protein_rd_t;

class protein_grouping_t;

class alpha_complex_model_t;

class pocket_model_t;

class protein_grouping_ui_model_t;

class protein_model_t:public glutils::renderable_t
{
public:
  friend class protein_model_ui_t;
  friend class viewer_mainwindow;

  enum eRenderModels
  {
    RMDL_SPACE_FILL,
    RMDL_BALL_STICK,
    RMDL_SMALL_SPACE_FILL,
    RMDL_COUNT
  };

  enum eRenderModes
  {
    RMDE_BACKBONE,
    RMDE_FULL,
    RMDE_CALPHA,
    RMDE_POCKET_ATOMS,
    RMDE_NONE,
    RMDE_COUNT
  };

  enum eSurfaceRenderMode
  {
    SRM_HIDE,
    SRM_SHOW,
    SRM_SHOW_WIREFRAME,
    SRM_COUNT
  };

  enum eAlphaComplexRenderMode
  {
    ACRM_HIDE         = 0,
    ACRM_TETRAHEDRONS = 1,
    ACRM_TRIANGLES    = 2,
    ACRM_EDGES        = 4,
    ACRM_SHOWALL      = 7,
  };

  enum ePocketRenderMode
  {
    PRM_SHOW_NONE,
    PRM_SHOW_ONE,
    PRM_SHOW_ALL,
  };


  enum eSecondaryRenderModes
  {
    SEC_NONE=0,
    SEC_TUBES=1,
    SEC_SHEETS=2,
    SEC_HELICES=4,
    SEC_ALL=7,
  };

protected:

  boost::shared_ptr<protein_t>             m_protein;
  boost::shared_ptr<protein_rd_t>          m_protein_rd;
  boost::shared_ptr<protein_grouping_t>    m_protein_atoms_grouping;
  boost::shared_ptr<glutils::renderable_t> m_surface_renderer;
  boost::shared_ptr<alpha_complex_model_t> m_alpha_complex_model;
  boost::shared_ptr<pocket_model_t>        m_pocket_model;
  boost::shared_ptr<secondary_model_t>     m_secondary_model;

  // other state stuff
  eRenderModels                   m_render_model;
  eRenderModes                    m_render_mode;
  eSurfaceRenderMode              m_surface_render_mode;
  uint                            m_alpha_complex_render_mode;
  ePocketRenderMode               m_pocket_render_mode;
  eSecondaryRenderModes           m_sec_renderMode;

  uint                            m_pocket_alpha_num;
  uint                            m_pocket_num;

  double                 m_add_atom_radius;
  double                 m_alpha_value;

  glutils::material_properties_t  m_atombond_material;
  glutils::material_properties_t  m_surface_material;

  std::string            m_protein_name;

  void render_onelevel() const;
  void render_surface() const;
  void render_secondary() const;

  void update_pocket_render_state();
  void update_sf_model_for_pocket();

public:
  protein_model_t ( const std::string & pf);
  virtual ~protein_model_t ();

  int render();
  static void initgl();

  std::string name() const
  {
    return m_protein_name.c_str();
  }

  void updateSecondaryHelix(int no);
  void updateSecondarySheet(int no);
  bool get_extent ( double * );

private:
  void load_surface(const std::string &filename);

  void load_alpha_shape(const std::string &filename);

  void load_pockets(const std::string &pocf,const std::string &tetf);
};

typedef boost::shared_ptr<protein_model_t> protein_model_ptr_t;

#include <QAbstractListModel>

class protein_grouping_ui_model_t:public QAbstractItemModel
{
  Q_OBJECT

public:

  protein_grouping_ui_model_t(boost::shared_ptr<protein_grouping_t> p, QObject *parent = 0);

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

  uint get_grouping_type();

private:
  boost::shared_ptr<protein_grouping_t> m_protein_grouping;
};

#endif//__PROTIEN_MODEL_H

