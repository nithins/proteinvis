#ifndef QTGLMATERIALEDITOR_INCLUDED
#define QTGLMATERIALEDITOR_INCLUDED

#include <ui_qtglmaterialeditor.h>

#include <QGroupBox>

#include <glutils.h>

class QtGlMaterialEditor:public QGroupBox,public Ui::GlMaterialEditor
{
  Q_OBJECT;

public:

  QtGlMaterialEditor(QWidget *par);

  void setMaterial(const glutils::material_properties_t & m);

  glutils::material_properties_t getMaterial();

Q_SIGNALS:
    void materialChanged(const glutils::material_properties_t &);

private Q_SLOTS:

  void on_ambient_colorpicker_colorChanged(const QColor &);
  void on_diffuse_colorpicker_colorChanged(const QColor &);
  void on_specular_colorpicker_colorChanged(const QColor &);
  void on_emission_colorpicker_colorChanged(const QColor &);
  void on_shininess_spinBox_valueChanged ( int );

private:
  void emit_material_changed();
};

#endif
