#ifndef QTGLLIGHTEDITOR_INCLUDED
#define QTGLLIGHTEDITOR_INCLUDED

#include <ui_qtgllighteditor.h>

#include <QGroupBox>

#include <glutils.h>

#include <vector>

class QtGlLightEditor:public QGroupBox,public Ui::GlLightEditor
{
  Q_OBJECT;

public://typedefs

  typedef glutils::light_properties_t     light_properties_t;

  typedef std::vector<light_properties_t> light_properties_list_t;

public: // ctor/dtors

  QtGlLightEditor(QWidget *par);

public:// setter/getters

  void set_num_lights(int numlights);

  void set_light(int lightnum , const light_properties_t & lp);

Q_SIGNALS:
    void lightChanged(int i,const glutils::light_properties_t &);

private Q_SLOTS:

  void on_ambient_colorpicker_colorChanged(const QColor &);
  void on_diffuse_colorpicker_colorChanged(const QColor &);
  void on_specular_colorpicker_colorChanged(const QColor &);
  void on_type_comboBox_currentIndexChanged(int);
  void on_number_spinBox_valueChanged(int i);


private:
  light_properties_list_t m_lights;

private:
  light_properties_t & get_active_light();

  int get_active_light_idx();

  void push_to_ui(const light_properties_t &l);

};

#endif
