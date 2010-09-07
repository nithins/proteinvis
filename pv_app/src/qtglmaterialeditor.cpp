#include <qtglmaterialeditor.h>

QtGlMaterialEditor::QtGlMaterialEditor(QWidget *par)
{
  setParent(par);
  setupUi(this);
}

void QtGlMaterialEditor::on_ambient_colorpicker_colorChanged(const QColor &)
{
  emit_material_changed();
}

void QtGlMaterialEditor::on_diffuse_colorpicker_colorChanged(const QColor &)
{
  emit_material_changed();
}

void QtGlMaterialEditor::on_specular_colorpicker_colorChanged(const QColor &)
{
  emit_material_changed();
}

void QtGlMaterialEditor::on_emission_colorpicker_colorChanged(const QColor &)
{
  emit_material_changed();
}

void QtGlMaterialEditor::on_shininess_spinBox_valueChanged ( int )
{
  emit_material_changed();
}

glutils::color4f_t from_qcolor(const QColor qc)
{
  glutils::color4f_t c;

  c[0] = qc.redF();
  c[1] = qc.greenF();
  c[2] = qc.blueF();
  c[3] = qc.alphaF();

  return c;
}

QColor to_qcolor(const glutils::color4f_t c)
{
  return QColor::fromRgbF(c[0],c[1],c[2],c[3]);
}

void QtGlMaterialEditor::setMaterial(const glutils::material_properties_t & m)
{
  blockSignals(true);

  ambient_colorpicker->setCurrentColor(to_qcolor(m.ambient));
  diffuse_colorpicker->setCurrentColor(to_qcolor(m.diffuse));
  specular_colorpicker->setCurrentColor(to_qcolor(m.specular));
  emission_colorpicker->setCurrentColor(to_qcolor(m.emission));
  shininess_spinBox->setValue(m.shininess);

  blockSignals(false);
}

glutils::material_properties_t QtGlMaterialEditor::getMaterial()
{
  glutils::material_properties_t m;

  m.ambient   = from_qcolor(ambient_colorpicker->currentColor());
  m.diffuse   = from_qcolor(diffuse_colorpicker->currentColor());
  m.specular  = from_qcolor(specular_colorpicker->currentColor());
  m.emission  = from_qcolor(emission_colorpicker->currentColor());
  m.shininess = shininess_spinBox->value();

  return m;
}

void QtGlMaterialEditor::emit_material_changed()
{
  emit materialChanged(getMaterial());
}
