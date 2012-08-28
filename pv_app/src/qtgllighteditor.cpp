#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/typeof/typeof.hpp>

#include <cpputils.h>

#include <qtgllighteditor.h>

using namespace boost;
using namespace boost::lambda;

typedef QtGlLightEditor::light_properties_t      light_properties_t;
typedef QtGlLightEditor::light_properties_list_t light_properties_list_t;

enum eLightType
{
  LIGHTTYPE_DIRECTIONAL,
  LIGHTTYPE_POSITIONAL,
  LIGHTTYPE_NONE
};

inline eLightType get_light_type(const light_properties_t &l)
{
  if(l.enabled == false)
    return LIGHTTYPE_NONE;

  if(l.position[3] == 0)
    return LIGHTTYPE_DIRECTIONAL;
  else
    return LIGHTTYPE_POSITIONAL;
}

inline glutils::color4f_t from_qcolor(const QColor qc)
{
  return la::make_vec<float>(qc.redF(),qc.greenF(),qc.blueF(),qc.alphaF());
}

inline QColor to_qcolor(const glutils::color4f_t c)
{
  return QColor::fromRgbF(c[0],c[1],c[2],c[3]);
}

// ctor/dtors
QtGlLightEditor::QtGlLightEditor(QWidget *par)
{
  setParent(par);
  setupUi(this);
}

// setter/getters
void QtGlLightEditor::set_num_lights(int numlights)
{
  int os = m_lights.size();

  m_lights.resize(numlights);

  number_spinBox->setMaximum(m_lights.size()-1);

  number_spinBox->setMinimum(std::min((size_t)0,m_lights.size()-1));

  if(os >= numlights) return;

  BOOST_AUTO(disable_light_ftor,bind(&light_properties_t::enabled,_1)=false);

  std::for_each(m_lights.begin()+os,m_lights.end(),disable_light_ftor);
}

void QtGlLightEditor::set_light(int lightnum , const light_properties_t & lp)
{
  m_lights[lightnum] = lp;

  if(lightnum == get_active_light_idx())
    push_to_ui(get_active_light());
}

// slots
void QtGlLightEditor::on_ambient_colorpicker_colorChanged(const QColor &c)
{
  get_active_light().ambient = from_qcolor(c);
  emit lightChanged(get_active_light_idx(),get_active_light());
}

void QtGlLightEditor::on_diffuse_colorpicker_colorChanged(const QColor &c)
{
  get_active_light().diffuse = from_qcolor(c);
  emit lightChanged(get_active_light_idx(),get_active_light());
}

void QtGlLightEditor::on_specular_colorpicker_colorChanged(const QColor &c)
{
  get_active_light().specular = from_qcolor(c);
  emit lightChanged(get_active_light_idx(),get_active_light());
}

void QtGlLightEditor::on_type_comboBox_currentIndexChanged(int i)
{
  glutils::light_properties_t &l = get_active_light();

  l.enabled     = (i != LIGHTTYPE_NONE);
  l.position[3] = (i == LIGHTTYPE_DIRECTIONAL)?(0):(1);

  colors_gridLayout->setEnabled(l.enabled);
  emit lightChanged(get_active_light_idx(),get_active_light());
}

void QtGlLightEditor::on_number_spinBox_valueChanged(int)
{
  push_to_ui(get_active_light());
}

// private utility
int QtGlLightEditor::get_active_light_idx()
{
  return number_spinBox->value();
}

light_properties_t & QtGlLightEditor::get_active_light()
{
  return m_lights[get_active_light_idx()];
}

void QtGlLightEditor::push_to_ui(const light_properties_t &l)
{
  blockSignals(true);

  type_comboBox->setCurrentIndex(get_light_type(l));
  ambient_colorpicker->setCurrentColor(to_qcolor(l.ambient));
  diffuse_colorpicker->setCurrentColor(to_qcolor(l.diffuse));
  specular_colorpicker->setCurrentColor(to_qcolor(l.specular));
  colors_gridLayout->setEnabled(l.enabled);

  blockSignals(false);
}
