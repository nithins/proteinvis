/***************************************************************************
 *   Copyright (C) 2012 by Pranav D Bagur,
 *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <fstream>
#include <string>
#include <sstream>
#include <set>
#include <GL/glew.h>
#include <GLSLProgram.h>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>

#include <logutil.h>
#include <cpputils.h>

#include <protein.h>
#include <secondaryModel.h>
#include <malloc.h>
#include <math.h>
#include <pv_config.h>

using  namespace glutils;

// auto generated in config.h .. cmake does it
const int g_segs_btw_ctrlPts = SECONDARY_NUM_SPLINESEGS;

GLSLProgram *s_sheetShader = NULL;
GLSLProgram *s_sheetTipsShader = NULL;
GLSLProgram *s_tubeShader = NULL;
GLSLProgram *s_tubeCapShader = NULL;
GLSLProgram *s_helixShader = NULL;
GLSLProgram *s_helixCapShader = NULL;
GLSLProgram *s_helixImposterShader = NULL;

secondary_model_t::secondary_model_t(boost::shared_ptr<protein_t> protein)
{
  this->m_protein=protein;

  m_chains_rd.resize(m_protein->get_num_chains());

  InitSplines();
  InitSheets();
  InitHelices();
  InitLoops();
}

secondary_model_t::~secondary_model_t()
{
  //destructor
}

vertex_t Interpolate(vertex_t p,vertex_t q,vertex_t r,vertex_t s,double t)
{
  typedef n_vector_t<double,4> vertex4_t;

  vertex4_t t_row(t*t*t,t*t,t,1);

  return (
      p * dot_product(t_row, vertex4_t(-1, 3,-3, 1)) +
      q * dot_product(t_row, vertex4_t( 3,-6, 0, 4)) +
      r * dot_product(t_row, vertex4_t(-3, 3, 3, 1)) +
      s * dot_product(t_row, vertex4_t( 1, 0, 0, 0))
        )/6;
}

void DetailPtGen(const vertex_t &p,const vertex_t &q,
                 const vertex_t &r,const vertex_t &s,
                 vertex_list_t &spts)
{
  for (int i=0;i<g_segs_btw_ctrlPts;i++)
    spts.push_back(Interpolate(p,q,r,s,double(i)/double(g_segs_btw_ctrlPts)));
}

void BSplines(vertex_t *cpts,const int & num_cpts,vertex_list_t &spts)
{
  DetailPtGen(2*cpts[0] - cpts[1],cpts[0],cpts[1],cpts[2],spts);

  for (int i=1;i<num_cpts-2;i++)
    DetailPtGen(cpts[i-1],cpts[i+0],cpts[i+1],cpts[i+2],spts);

  DetailPtGen(cpts[num_cpts-3],cpts[num_cpts-2],cpts[num_cpts-1],
              2*cpts[num_cpts-1]-cpts[num_cpts-2],spts);

  spts.push_back(Interpolate(cpts[num_cpts-3],cpts[num_cpts-2],cpts[num_cpts-1],
                             2*cpts[num_cpts-1]-cpts[num_cpts-2],1.0));
}

inline vertex_t atom_to_vertex(const atom_t & a)
{
  return vertex_t(a.x,a.y,a.z);
}

void secondary_model_t::InitShaders()
{
  if(s_sheetShader != NULL)
    return;

  //initialize sheet shaders
  string sheet_log,sheet_tips_log;

  QFile sheet_vert ( ":/shaders/sheet_vert.glsl" );
  QFile sheet_geom ( ":/shaders/sheet_geom.glsl" );
  QFile sheet_frag ( ":/shaders/sheet_frag.glsl" );

  sheet_vert.open ( QIODevice::ReadOnly );
  sheet_geom.open ( QIODevice::ReadOnly );
  sheet_frag.open ( QIODevice::ReadOnly );

  QString sheet_vert_str =sheet_vert.readAll();
  QString sheet_geom_str =sheet_geom.readAll();
  QString sheet_frag_str =sheet_frag.readAll();

  sheet_vert.close();
  sheet_geom.close();
  sheet_frag.close();

  s_sheetShader = GLSLProgram::createFromSourceStrings
      (
        sheet_vert_str.toStdString(),
        sheet_geom_str.toStdString(),
        sheet_frag_str.toStdString(),
        GL_LINE_STRIP_ADJACENCY,
        GL_TRIANGLE_STRIP
        );

  s_sheetShader->GetProgramLog ( sheet_log );

  _LOG_VAR ( sheet_log );

  QString sheet_tips_vert_str = sheet_vert_str;
  QString sheet_tips_geom_str = sheet_geom_str;

  sheet_tips_geom_str.replace("//#define ENABLE_TIPS","#define ENABLE_TIPS");

  s_sheetTipsShader = GLSLProgram::createFromSourceStrings
      (
        sheet_tips_vert_str.toStdString(),
        sheet_tips_geom_str.toStdString(),
        sheet_frag_str.toStdString(),
        GL_LINE_STRIP_ADJACENCY,
        GL_TRIANGLE_STRIP
        );

  s_sheetTipsShader->GetProgramLog ( sheet_tips_log );

  _LOG_VAR (sheet_tips_log);

  //initialize helix shader
  string helix_log,helix_cap_log;

  QFile helix_vert ( ":/shaders/helix_vert.glsl" );
  QFile helix_geom ( ":/shaders/helix_geom.glsl" );
  QFile helix_frag ( ":/shaders/helix_frag.glsl" );

  helix_vert.open ( QIODevice::ReadOnly );
  helix_geom.open ( QIODevice::ReadOnly );
  helix_frag.open ( QIODevice::ReadOnly );

  QString helix_vert_str =helix_vert.readAll();
  QString helix_geom_str =helix_geom.readAll();
  QString helix_frag_str =helix_frag.readAll();

  s_helixShader = GLSLProgram::createFromSourceStrings
      (
        helix_vert_str.toStdString(),
        helix_geom_str.toStdString(),
        helix_frag_str.toStdString(),
        GL_LINE_STRIP_ADJACENCY,
        GL_TRIANGLE_STRIP
        );

  helix_geom_str.replace("//#define NEED_CAPS","#define NEED_CAPS");

  s_helixCapShader = GLSLProgram::createFromSourceStrings
          (
            helix_vert_str.toStdString(),
            helix_geom_str.toStdString(),
            helix_frag_str.toStdString(),
            GL_LINE_STRIP_ADJACENCY,
            GL_TRIANGLE_STRIP
            );



  s_helixShader->GetProgramLog ( helix_log );
  s_helixCapShader->GetProgramLog ( helix_log );

  helix_vert.close();
  helix_geom.close();
  helix_frag.close();

  if( helix_log.find("error") != string::npos)
    throw std::runtime_error("failed compiling helix shader\n"+helix_log);

  if( helix_cap_log.find("error") != string::npos)
    throw std::runtime_error("failed compiling helix cap shader\n"+helix_cap_log);

  //initialize tubes shader
  string cyl_adj_log;

  QFile cyl_adj_vert ( ":/shaders/cylinder_adj_vert.glsl" );
  QFile cyl_adj_geom ( ":/shaders/cylinder_adj_geom.glsl" );
  QFile cyl_adj_frag ( ":/shaders/cylinder_adj_frag.glsl" );

  cyl_adj_vert.open ( QIODevice::ReadOnly );
  cyl_adj_geom.open ( QIODevice::ReadOnly );
  cyl_adj_frag.open ( QIODevice::ReadOnly );

  QString cyl_vert_src = cyl_adj_vert.readAll();
  QString cyl_geom_src = cyl_adj_geom.readAll();
  QString cyl_frag_src = cyl_adj_frag.readAll();

  cyl_adj_vert.close();
  cyl_adj_geom.close();
  cyl_adj_frag.close();

  s_tubeShader = GLSLProgram::createFromSourceStrings
      (
        cyl_vert_src.toStdString(),
        cyl_geom_src.toStdString(),
        cyl_frag_src.toStdString(),
        GL_LINES_ADJACENCY,
        GL_TRIANGLE_STRIP
        );


  s_tubeShader->GetProgramLog ( cyl_adj_log );

  if( cyl_adj_log.find("error") != string::npos)
    throw std::runtime_error("failed compiling cylinder shader\n"+cyl_adj_log);

  //initialize tubes shader
  string cyl_cap_log;

  cyl_geom_src.replace("//#define ENABLE_CAPS","#define ENABLE_CAPS");
  cyl_frag_src.replace("//#define ENABLE_CAPS","#define ENABLE_CAPS");

  s_tubeCapShader = GLSLProgram::createFromSourceStrings
      (
        cyl_vert_src.toStdString(),
        cyl_geom_src.toStdString(),
        cyl_frag_src.toStdString(),
        GL_TRIANGLES,
        GL_TRIANGLE_STRIP
        );

  s_tubeCapShader->GetProgramLog ( cyl_cap_log );

  if( cyl_cap_log.find("error") != string::npos)
    throw std::runtime_error("failed compiling cylinder caps shader\n"+cyl_cap_log);

  //initialize tubes shader
  string helix_imposter_log;

  QFile helix_imposter_vert ( ":/shaders/helix_ideal_vert.glsl");
  QFile helix_imposter_geom ( ":/shaders/helix_ideal_geom.glsl");
  QFile helix_imposter_frag ( ":/shaders/helix_ideal_frag.glsl");


  helix_imposter_vert.open ( QIODevice::ReadOnly );
  helix_imposter_geom.open ( QIODevice::ReadOnly );
  helix_imposter_frag.open ( QIODevice::ReadOnly );

  s_helixImposterShader = GLSLProgram::createFromSourceStrings
      (
        string ( helix_imposter_vert.readAll().constData() ),
        string ( helix_imposter_geom.readAll().constData() ),
        string ( helix_imposter_frag.readAll().constData() ),
        GL_TRIANGLES,
        GL_TRIANGLE_STRIP
        );

  helix_imposter_vert.close();
  helix_imposter_geom.close();
  helix_imposter_frag.close();

  s_helixImposterShader->GetProgramLog ( helix_imposter_log );

  _LOG_VAR ( helix_imposter_log );
}

void secondary_model_t::InitSplines()
{
  color_list_t chain_colors;

  for(int i=0;i<m_protein->get_num_chains();i++)
    chain_colors.push_back(color_t(double(rand()%128)/128.0f,
                                   double(rand()%128)/128.0f,
                                   double(rand()%128)/128.0f));


  for (int i = 0 ;i  < m_protein->get_num_chains(); ++i)
  {
    vertex_list_t caatom_pos;
    vertex_list_t  oatom_pos;

    chain_t  chain  = m_protein->get_chains()[i];
    acid_t  acid_s  = m_protein->get_acids()[chain.start];
    acid_t  acid_e  = m_protein->get_acids()[chain.end-1];

    int atom_idx_s  = acid_s.start;
    int atom_idx_e  = acid_e.end;

    for(int j = atom_idx_s ; j < atom_idx_e ; ++j)
    {
      atom_t atom = m_protein->get_atoms()[j];

      if(m_protein->is_ca_atom(j))
        caatom_pos.push_back(atom_to_vertex(atom));

      if(m_protein->is_o_atom(j))
        oatom_pos.push_back(atom_to_vertex(atom));
    }

    assert(oatom_pos.size() == caatom_pos.size());

    normal_t prev_n(0,0,0);

    for(int j = 0 ; j < oatom_pos.size(); ++j)
    {
      normal_t n = euclid_normalize(oatom_pos[j]-caatom_pos[j]);

      if(dot_product(prev_n,n) <0 )
        n *=-1;

      oatom_pos[j] = caatom_pos[j]+n;
      prev_n = n;
    }

    BSplines(caatom_pos.data(),caatom_pos.size(),m_chains_rd[i].spline_pts);
    BSplines(oatom_pos.data() ,oatom_pos.size() ,m_chains_rd[i].sec_spline_pts);

    m_chains_rd[i].spline_pts_bo     = make_buf_obj(m_chains_rd[i].spline_pts);
    m_chains_rd[i].sec_spline_pts_bo = make_buf_obj(m_chains_rd[i].sec_spline_pts);
    m_chains_rd[i].color             = chain_colors[i];
  }
}

void secondary_model_t::InitLoops()
{

  for (int i = 0 ;i  < m_chains_rd.size(); ++i)
  {
    chain_rd_t &chain_rd = m_chains_rd[i];

    quad_idx_list_t &loop_idxs     = chain_rd.loop_idxs;
    tri_idx_list_t  &loop_cap_idxs = chain_rd.loop_cap_idxs;
    vertex_list_t   &spts          = chain_rd.spline_pts;

    vector<bool>    is_loop_pt(spts.size(),true);

    for(int j = 0 ; j < m_helices_rd.size(); ++j )
    {
      if(m_helices_rd[j].chainno !=i)
        continue;

      int b = m_helices_rd[j].spt_idx_b;
      int e = m_helices_rd[j].spt_idx_e;

      fill(is_loop_pt.begin()+b+g_segs_btw_ctrlPts,
           is_loop_pt.begin()+e-1.5*g_segs_btw_ctrlPts,false);
    }

    for(int j = 0 ; j < m_strands_rd.size(); ++j )
    {
      if(m_strands_rd[j].chainno !=i)
        continue;

      int b = m_strands_rd[j].spt_idx_b;
      int e = m_strands_rd[j].spt_idx_e;

      fill(is_loop_pt.begin()+b+1,
           is_loop_pt.begin()+e-1,false);
    }

    if(is_loop_pt[0] && is_loop_pt[1])
    {
      loop_cap_idxs.push_back(tri_idx_t(0,1,2));
    }

    if(is_loop_pt[spts.size()-2] && is_loop_pt[spts.size()-1])
    {
      loop_cap_idxs.push_back(tri_idx_t(spts.size()-1,spts.size()-2,spts.size()-3));
    }


    for(int j = 1 ; j < spts.size()-2; ++j)
    {
      if(!(is_loop_pt[j] && is_loop_pt[j+1]))
        continue;

      if(is_loop_pt[j-1] && is_loop_pt[j+2])
      {
        loop_idxs.push_back(quad_idx_t(j-1,j,j+1,j+2));
      }
      else if( is_loop_pt[j-1])
      {
        loop_cap_idxs.push_back(tri_idx_t(j+1,j,j-1));
      }
      else if( is_loop_pt[j+2])
      {
        loop_cap_idxs.push_back(tri_idx_t(j,j+1,j+2));
      }
    }

    chain_rd.loop_idxs_bo     = make_buf_obj(loop_idxs);
    chain_rd.loop_caps_idxs_bo= make_buf_obj(loop_cap_idxs);
  }
}


void secondary_model_t::InitSheets()
{
  //get the number of sheets
  int num_strands = m_protein->get_num_sheets();

  if (num_strands == 0)
    return;


  set<int> sheetno_set;
  for(int i=0;i<num_strands;i++)
    sheetno_set.insert(m_protein->get_sheets()[i].sheet_no);

  int num_sheets  = sheetno_set.size();

  color_list_t sheet_colors;

  for(int i=0;i<num_sheets;i++)
    sheet_colors.push_back(color_t(double(rand()%128)/128.0f,
                                   double(rand()%128)/128.0f,
                                   double(rand()%128)/128.0f));

  m_strands_rd.resize(num_strands);

  for(int i=0;i<num_strands;i++)
  {
    const sheet_t & strand = m_protein->get_sheets()[i];
    strand_rd_t &strand_rd = m_strands_rd[i];

    if(strand.start_chainno+1 != strand.end_chainno)
    {
      cout<<"init chain is not same as term"<<endl;
      cout<<"code for sheet fragements from multiple chains needed"<<endl;
      throw std::logic_error("missing code");
    }

    int chainno = strand.start_chainno;
    const chain_rd_t &chain_rd = m_chains_rd[chainno];

    int res_b = strand.start_resno;
    int res_e = strand.end_resno;

    int chain_res_b = m_protein->get_chains()[chainno].start;
    int chain_res_e = m_protein->get_chains()[chainno].end;

    assert(chain_res_b <= res_b && res_b <  chain_res_e);
    assert(chain_res_b <= res_e && res_e <= chain_res_e);

    res_b = res_b - chain_res_b;
    res_e = res_e - chain_res_b;

    int spt_b = res_b*g_segs_btw_ctrlPts;
    int spt_e = res_e*g_segs_btw_ctrlPts;

    assert(0<=spt_b && spt_b <  chain_rd.spline_pts.size());
//    assert(0< spt_e && spt_e <= m_chains_rd[chainno].spline_pts.size());
    spt_e  = min<int>(spt_e,chain_rd.spline_pts.size());

    strand_rd.chainno     = chainno;
    strand_rd.spt_idx_b   = spt_b;
    strand_rd.spt_idx_e   = spt_e;
    strand_rd.color       = sheet_colors[strand.sheet_no];

    vector<double> &width = strand_rd.width;
    width.resize(spt_e-spt_b);

    int tpts = g_segs_btw_ctrlPts*3/4;

    for(int j = tpts; j< (spt_e-spt_b)-tpts;++j )
      width[j] = 1.2;

    for(int j = 0; j< tpts;++j )
      *(width.begin()+j) = 0.4+double(j)/double(tpts)*0.8;

    for(int j = 0; j< tpts;++j )
      *(width.end()-(tpts-j)) = 2.4-double(j)/double(tpts-1)*2;

    strand_rd.width_bo= buf_obj_t::create_bo
        (width.data(),GL_DOUBLE,1,GL_ARRAY_BUFFER,width.size()*sizeof(double),0);

  }
}

inline vertex_t line_plane_ixn(vertex_t ld, vertex_t lp,vertex_t pn, vertex_t pp)
{
  double t = (dot_product(pn,pp) - dot_product(pn,lp))/(dot_product(pn,ld));
  return lp + t*ld;
}

inline vertex_t closest_line_pt(vertex_t ld,vertex_t lp,vertex_t pt)
{
  return lp+ld*dot_product(pt-lp,ld)/dot_product(ld,ld);
}

inline vertex_t closest_plane_pt(vertex_t pn,vertex_t pp,vertex_t pt)
{
  return line_plane_ixn(pn,pt,pn,pp);
}

void secondary_model_t::InitHelices()
{
  int num_helices = m_protein->get_num_helices();

  if(num_helices ==0)
    return;

  m_helices_rd.resize(num_helices);

  color_list_t chain_colors;

  for(int i = 0 ; i < m_protein->get_num_chains(); ++i)
    chain_colors.push_back(color_t(double(rand()%128)/128.0f,
                                   double(rand()%128)/128.0f,
                                   double(rand()%128)/128.0f));


  for(int i=0;i<m_protein->get_num_helices();i++)
  {
    const helix_t & helix = m_protein->get_helices()[i];

    if(helix.start_chainno +1!= helix.end_chainno)
    {
      cout<<"init chain is not same as term"<<endl;
      cout<<"code for helix fragemnts from multiple chains needed"<<endl;
      throw std::logic_error("missing code");
    }

    const int &chainno  = helix.start_chainno;
    vertex_list_t &spts = m_chains_rd[chainno].spline_pts;

    int res_b = helix.start_resno;
    int res_e = helix.end_resno;

    int chain_res_b = m_protein->get_chains()[chainno].start;
    int chain_res_e = m_protein->get_chains()[chainno].end;

    assert(chain_res_b <= res_b && res_b<  chain_res_e);
    assert(chain_res_b <= res_e && res_e<= chain_res_e);

    res_b = res_b - chain_res_b;
    res_e = res_e - chain_res_b;

    int spt_b = res_b*g_segs_btw_ctrlPts;//  + g_segs_btw_ctrlPts/2;
    int spt_e = res_e*g_segs_btw_ctrlPts;//  - g_segs_btw_ctrlPts;

    int b_offset = g_segs_btw_ctrlPts/2;
    int e_offset = g_segs_btw_ctrlPts;

    assert(0<=spt_b && spt_b <  m_chains_rd[chainno].spline_pts.size());
    assert(0< spt_e && spt_e <= m_chains_rd[chainno].spline_pts.size());

    m_helices_rd[i].spt_idx_b      = spt_b;
    m_helices_rd[i].spt_idx_e      = spt_e;
    m_helices_rd[i].chainno        = chainno;
    m_helices_rd[i].color          = chain_colors[chainno];

    // to determine the helix axis dir a triple of pts is used to compute
    // the binormal vector at point #2 of the triple.
    // a second binormal at a point diametrically opposite on the helix
    // is added so that the component not aligning with the axis is cancelled
    // we know that the helix turns for every 3.6 residues
    // so after 1.8 turns we get the pt diametrically opposite
    // the average of all such dirs is used

    normal_t axis_dir(0,0,0);

    int spt_jmp  = ceil(double(g_segs_btw_ctrlPts)/6.0);
    int spt_opp  = round(double(g_segs_btw_ctrlPts)*1.8);

    for(int j = spt_b + b_offset + spt_jmp; j < spt_e-e_offset-spt_jmp-spt_opp; ++j)
    {
      vertex_t p = spts[j-spt_jmp];
      vertex_t q = spts[j];
      vertex_t r = spts[j+spt_jmp];

      vertex_t u = spts[j+spt_opp-spt_jmp];
      vertex_t v = spts[j+spt_opp];
      vertex_t w = spts[j+spt_opp+spt_jmp];

      normal_t  n1 = euclid_normalize(cross_product(r-q,p-q));
      normal_t  n2 = euclid_normalize(cross_product(w-v,u-v));

      axis_dir += euclid_normalize(n1+n2);
    }

    axis_dir = euclid_normalize(axis_dir);

    // a point on the axis is determined by projecting all points on the
    // helix to a plane with normal helix_dir and then taking their mean

    vertex_t axis_pt(0,0,0);
    for(int j = spt_b + b_offset; j < spt_e-e_offset; ++j)
      axis_pt += closest_plane_pt(axis_dir,vertex_t(0,0,0),spts[j]);
    axis_pt /= spt_e-spt_b-e_offset-b_offset;

    // the radius is the average distance from the mean to all the projections
    double radius = 0;
    for(int j = spt_b+b_offset; j < spt_e-e_offset; ++j)
      radius += euclid_norm(closest_plane_pt(axis_dir,vertex_t(0,0,0),spts[j])-axis_pt);
    radius /= spt_e-spt_b-e_offset-b_offset;

    // the first and last points of the spline section is projected to the axis
    m_helices_rd[i].axis_b   = closest_line_pt(axis_dir,axis_pt,spts[spt_b]);
    m_helices_rd[i].axis_e   = closest_line_pt(axis_dir,axis_pt,spts[spt_e-1]);
    m_helices_rd[i].radius   = radius;
    m_helices_rd[i].axis_dir = axis_dir;
    m_helices_rd[i].x_dir    = radius*euclid_normalize(spts[spt_b] - m_helices_rd[i].axis_b);
    m_helices_rd[i].y_dir    = cross_product(axis_dir,m_helices_rd[i].x_dir);
  }


}

void secondary_model_t::RenderSheets()
{
  glPushAttrib ( GL_ENABLE_BIT );

  GLuint Width_ATTR = s_sheetShader->getAttributeLocation("Width");

  s_sheetShader->use();
  for(int i = 0; i < m_strands_rd.size(); ++i)
  {
    strand_rd_t &strand = m_strands_rd[i];
    chain_rd_t  &chain  = m_chains_rd[strand.chainno];

    int offset  = strand.spt_idx_b;
    int count   = strand.spt_idx_e-strand.spt_idx_b;
    color_t col = strand.color;

    glColor3d(col[0],col[1],col[2]);

    chain.spline_pts_bo->bind_to_vertex_pointer(offset);
    chain.sec_spline_pts_bo->bind_to_normal_pointer(offset);
    strand.width_bo->bind_to_vertex_attrib_pointer(Width_ATTR);

    glDrawArrays(GL_LINE_STRIP_ADJACENCY,0,count);

    strand.width_bo->unbind_from_vertex_attrib_pointer(Width_ATTR);
    chain.sec_spline_pts_bo->unbind_from_normal_pointer();
    chain.spline_pts_bo->unbind_from_vertex_pointer();
  }

  s_sheetShader->disable();
  s_sheetTipsShader->use();

  Width_ATTR = s_sheetTipsShader->getAttributeLocation("Width");

  for(int i = 0; i < m_strands_rd.size(); ++i)
  {
    strand_rd_t &strand = m_strands_rd[i];
    chain_rd_t  &chain  = m_chains_rd[strand.chainno];

    int beg   = strand.spt_idx_b;
    int end   = strand.spt_idx_e;

    vertex_list_t  & spts= chain.spline_pts;
    vertex_list_t  &sspts= chain.sec_spline_pts;
    vector<double> &width= strand.width;

    color_t col = strand.color;

    glColor3d(col[0],col[1],col[2]);

    glBegin(GL_LINES_ADJACENCY);

    for(int j = 0 ; j <4 ;++j)
    {
      glVertexAttrib1f(Width_ATTR,width[j]);
      glNormal3f(sspts[beg+j][0],sspts[beg+j][1],sspts[beg+j][2]);
      glVertex3f(spts [beg+j][0],spts [beg+j][1],spts [beg+j][2]);
    }

    for(int j = 0 ; j <4 ;++j)
    {
      int karr[] = {end-3,end-2,end-1,end-1};
      int k = karr[j];

      glVertexAttrib1f(Width_ATTR,width[k-beg]);
      glNormal3f(sspts[k][0],sspts[k][1],sspts[k][2]);
      glVertex3f(spts [k][0],spts [k][1],spts [k][2]);
    }

    for(int j = 0 ; j <4 ;++j)
    {
      int karr[] = {end-2,end-1,end-1,end-1};
      int k = karr[j];

      glVertexAttrib1f(Width_ATTR,width[k-beg]);
      glNormal3f(sspts[k][0],sspts[k][1],sspts[k][2]);
      glVertex3f(spts [k][0],spts [k][1],spts [k][2]);
    }


    glEnd();
  }


  s_sheetTipsShader->disable();
  glPopAttrib();
}



void secondary_model_t::RenderHelices()
{
  glPushAttrib ( GL_ENABLE_BIT );

  s_helixShader->use();

  for(int i = 0; i < m_helices_rd.size(); ++i)
  {
    helix_rd_t &helix_rd = m_helices_rd[i];

    color_t col = helix_rd.color;
    glColor3d(col[0],col[1],col[2]);


    bufobj_ptr_t bo = m_chains_rd[helix_rd.chainno].spline_pts_bo;

    int offset = helix_rd.spt_idx_b+1;
    int count  = helix_rd.spt_idx_e-helix_rd.spt_idx_b-2;
    normal_t n = helix_rd.axis_dir;

    s_helixShader->sendUniform("g_helixUp",(float)n[0],(float)n[1],(float)n[2]);
    bo->bind_to_vertex_pointer();
    glDrawArrays(GL_LINE_STRIP_ADJACENCY,offset,count);
    bo->unbind_from_vertex_pointer();
  }
  s_helixShader->disable();

  s_helixCapShader->use();

  for(int i = 0; i < m_helices_rd.size(); ++i)
  {
    helix_rd_t &helix_rd = m_helices_rd[i];

    color_t col = helix_rd.color;
    glColor3d(col[0],col[1],col[2]);


    vertex_list_t &spts= m_chains_rd[helix_rd.chainno].spline_pts;

    int beg = helix_rd.spt_idx_b;
    int end  = helix_rd.spt_idx_e;
    normal_t n = helix_rd.axis_dir;

    s_helixCapShader->sendUniform("g_helixUp",(float)n[0],(float)n[1],(float)n[2]);

    glBegin(GL_LINES_ADJACENCY);
    glVertex3f(spts[beg  ][0],spts[beg  ][1],spts[beg  ][2]);
    glVertex3f(spts[beg+1][0],spts[beg+1][1],spts[beg+1][2]);
    glVertex3f(spts[beg+2][0],spts[beg+2][1],spts[beg+2][2]);
    glVertex3f(spts[beg+3][0],spts[beg+3][1],spts[beg+3][2]);

    glVertex3f(spts[end-4][0],spts[end-4][1],spts[end-4][2]);
    glVertex3f(spts[end-3][0],spts[end-3][1],spts[end-3][2]);
    glVertex3f(spts[end-2][0],spts[end-2][1],spts[end-2][2]);
    glVertex3f(spts[end-1][0],spts[end-1][1],spts[end-1][2]);
    glEnd();

  }
  s_helixCapShader->disable();


  glPopAttrib();

}


void draw_plane(const vertex_t &pn,const vertex_t &pt,const double & S)
{
  vertex_t n = euclid_normalize(pn);

  vertex_t s1(0,0,1);

  if(dot_product(s1,n) >0.99 )
    s1 = vertex_t(0,1,0);

  vertex_t s2 = cross_product(n,s1);
  s1 = cross_product(s2,n);

  vertex_t p;

  glBegin(GL_QUADS);

  p = pt + S*(-1*s1 -s2); glVertex3f(p[0],p[1],p[2]);
  p = pt + S*(-1*s1 +s2); glVertex3f(p[0],p[1],p[2]);
  p = pt + S*(+1*s1 +s2); glVertex3f(p[0],p[1],p[2]);
  p = pt + S*(+1*s1 -s2); glVertex3f(p[0],p[1],p[2]);

  glEnd();
}

void secondary_model_t::RenderTubes()
{
  glPushAttrib ( GL_ENABLE_BIT );
  glDisable(GL_LIGHTING);
  s_tubeShader->use();

  for(int i = 0 ;i < m_protein->get_num_chains(); ++i)
  {
    const chain_rd_t &chain_rd = m_chains_rd[i];

//    color_t col(0.2,0.7,0.2);
//    glColor3d( col[0],col[1],col[2]);

    glColor3d( chain_rd.color[0],chain_rd.color[1],chain_rd.color[2]);

    const bufobj_ptr_t &vbo = chain_rd.spline_pts_bo;
    const bufobj_ptr_t &ibo = chain_rd.loop_idxs_bo;

    vbo->bind_to_vertex_pointer();
    glBindBuffer ( ibo->target(), ibo->id() );
    glDrawElements ( GL_LINES_ADJACENCY, ibo->get_num_items()*4, ibo->src_type(), 0 );
    glBindBuffer ( ibo->target(), 0);
    vbo->unbind_from_vertex_pointer();
  }

  s_tubeShader->disable();

  s_tubeCapShader->use();

  for(int i = 0 ;i < m_protein->get_num_chains(); ++i)
  {
    const chain_rd_t &chain_rd = m_chains_rd[i];

//    color_t col(0.2,0.7,0.2);
//    glColor3d( col[0],col[1],col[2]);
    glColor3d( chain_rd.color[0],chain_rd.color[1],chain_rd.color[2]);

    const bufobj_ptr_t &vbo = chain_rd.spline_pts_bo;
    const bufobj_ptr_t &ibo = chain_rd.loop_caps_idxs_bo;

    vbo->bind_to_vertex_pointer();
    glBindBuffer ( ibo->target(), ibo->id() );
    glDrawElements ( GL_TRIANGLES, ibo->get_num_items()*3, ibo->src_type(), 0);
    glBindBuffer ( ibo->target(), 0);
    vbo->unbind_from_vertex_pointer();
  }

  s_tubeCapShader->disable();

  glPopAttrib();
}

void secondary_model_t::RenderImposterHelices()
{
  glPushAttrib ( GL_ENABLE_BIT );

  s_helixImposterShader->use();

  for(int i = 0; i < m_helices_rd.size(); ++i)
  {
    helix_rd_t &helix_rd = m_helices_rd[i];

    color_t col = helix_rd.color;
    glColor3d(col[0],col[1],col[2]);

    vertex_t p = helix_rd.x_dir;
    vertex_t q = helix_rd.axis_b;
    vertex_t r = helix_rd.axis_e;

    glBegin(GL_TRIANGLES);

    glVertex3d(p[0],p[1],p[2]);
    glVertex3d(q[0],q[1],q[2]);
    glVertex3d(r[0],r[1],r[2]);

    glEnd();

  }
  s_helixImposterShader->disable();
  glPopAttrib();
}

