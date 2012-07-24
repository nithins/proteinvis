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
#include <boost/bind.hpp>

#include <logutil.h>
#include <cpputils.h>

#include <protein.h>
#include<secondaryModel.h>
#include <malloc.h>
#include <math.h>

using  namespace glutils;

secondary_model_t::secondary_model_t(boost::shared_ptr<protein_t> protein)
{
  this->m_protein=protein;



  atoms=m_protein->get_atoms();
  num_atoms=m_protein->get_num_atoms();

  InitShaders();
  InitSplines();
  InitSheets(0,true);
  InitHelices(0,true);
  InitTubes();

}

secondary_model_t::~secondary_model_t()
{
  //destructor
}

void secondary_model_t::GetCaAtoms(const uint **ca_atoms_idx_ref)
{
  const uint *ca_atoms_idx=*ca_atoms_idx_ref;
  atom_t *arr=new atom_t[num_ca_atoms];

  for(int i=0;i<num_ca_atoms;i++)
  {
    arr[i]=atoms[ca_atoms_idx[i]];
  }

  ca_atoms=arr;
}

void secondary_model_t::GetOAtoms(const uint **o_atoms_idx_ref)
{
  const uint *o_atoms_idx=*o_atoms_idx_ref;
  atom_t *arr=new atom_t[num_ca_atoms];
  std::string prevChainId=atoms[o_atoms_idx[0]].chainId;
  string currentChainId;
  for(int i=0;i<num_ca_atoms;i++)
  {
    arr[i]=atoms[o_atoms_idx[i]];
    std::stringstream out;
    out << arr[i].residueSeqNo;
    currentChainId=arr[i].chainId;
    if(currentChainId!=prevChainId)
    {
      //indexMap[(i-1)*g_segs_btw_ctrlPts]=(i+1)*g_segs_btw_ctrlPts;
      tubeMap[i-1]=i;
      prevChainId=currentChainId;
      //continue;
    }
    cAlphaMapping[currentChainId+out.str()]=i;

    prevChainId=currentChainId;
  }

  o_atoms=arr;
}

D3DXVECTOR3 secondary_model_t::Interpolate(D3DXVECTOR3 *point1,D3DXVECTOR3 *point2,D3DXVECTOR3 *point3,D3DXVECTOR3 *point4,float amount)
{
  D3DXVECTOR4 multip(0,0,0,0);
  D3DXVECTOR3 result(0,0,0);
  D3DXMATRIX bezierBasis(-1, 3, -3, 1, 3, -6, 3, 0, -3, 0, 3, 0, 1, 4, 1, 0);
  D3DXVECTOR4 tMat(amount*amount*amount,amount*amount,amount,1);
  mulMatrixVec(&multip,&tMat,&bezierBasis);


  D3DXVECTOR4 arg(point1->x,point2->x,point3->x,point4->x);
  result.x=(float)D3DXVec4Dot(&arg,&multip)/6;

  D3DXVECTOR4 arg2(point1->y,point2->y,point3->y,point4->y);
  result.y=(float)D3DXVec4Dot(&arg2,&multip)/6;

  D3DXVECTOR4 arg3(point1->z,point2->z,point3->z,point4->z);
  result.z=(float)D3DXVec4Dot(&arg3,&multip)/6;
  return result;
}

int secondary_model_t::DetailPtGen(D3DXVECTOR3 *point1,D3DXVECTOR3 *point2,D3DXVECTOR3 *point3,D3DXVECTOR3 *point4,D3DXVECTOR3 **retResult)
{
  list<D3DXVECTOR3> detailPts;


  for (int i=0;i<(int)g_segs_btw_ctrlPts;i++)
  {
    detailPts.push_back(Interpolate(point1,point2,point3,point4,(float)i/(float)g_segs_btw_ctrlPts));
  }

  int noOfPts=detailPts.size();

  D3DXVECTOR3 *arr=(D3DXVECTOR3 *)malloc(sizeof(D3DXVECTOR3)*noOfPts);
  copy(detailPts.begin(),detailPts.end(),arr);
  detailPts.clear();
  *retResult=arr;
  return noOfPts;
}

int secondary_model_t::BSplines(atom_t *atomPts,D3DXVECTOR3 **detailedRef,int atomPtsLength)
{
  list<D3DXVECTOR3> completeDetailedPtList;
  D3DXVECTOR3 *currentSegRef;
  D3DXVECTOR3 *controlPts=(D3DXVECTOR3* )malloc(sizeof(D3DXVECTOR3)*atomPtsLength);
  int noOfPts=0;
  int totalNoOfControlPts=0;
  int totalMissed=0;
  for (int i=0;i<atomPtsLength-3;i++)
  {


    if(tubeMap[i+3]!=0)
    {
      i=tubeMap[i+3]-1;
      indexMap[completeDetailedPtList.size()]=completeDetailedPtList.size()+1;
      totalMissed+=g_segs_btw_ctrlPts;
      continue;
    }

    missedMap[i]=completeDetailedPtList.size();
    controlPts[i]=D3DXVECTOR3((float)atomPts[i].x,(float)atomPts[i].y,(float)atomPts[i].z);
    controlPts[i+1]=D3DXVECTOR3((float)atomPts[i+1].x,(float)atomPts[i+1].y,(float)atomPts[i+1].z);
    controlPts[i+2]=D3DXVECTOR3((float)atomPts[i+2].x,(float)atomPts[i+2].y,(float)atomPts[i+2].z);
    controlPts[i+3]=D3DXVECTOR3((float)atomPts[i+3].x,(float)atomPts[i+3].y,(float)atomPts[i+3].z);


    noOfPts=DetailPtGen(&controlPts[i+0],&controlPts[i+1],&controlPts[i+2],&controlPts[i+3],&currentSegRef);
    totalNoOfControlPts+=noOfPts;

    for(int j=0;j<noOfPts;j++)
    {
      if(currentSegRef==NULL)
        break;
      completeDetailedPtList.push_back(currentSegRef[j]);
    }
    free(currentSegRef);
  }

  free(controlPts);
  D3DXVECTOR3 *arr=(D3DXVECTOR3* )malloc(sizeof(D3DXVECTOR3)*totalNoOfControlPts);
  copy(completeDetailedPtList.begin(),completeDetailedPtList.end(),arr);
  //*detailedRef=controlPts;


  *detailedRef=arr;


  //return atomPtsLength;
  return totalNoOfControlPts;
}






void secondary_model_t::InitShaders()
{

  //initialize sheet shader
  string sheet_log;

  QFile sheet_vert ( ":/shaders/sheet.vert" );
  QFile sheet_geom ( ":/shaders/sheet.geom" );
  QFile sheet_frag ( ":/shaders/sheet.frag" );

  sheet_vert.open ( QIODevice::ReadOnly );
  sheet_geom.open ( QIODevice::ReadOnly );
  sheet_frag.open ( QIODevice::ReadOnly );

  s_sheetShader = GLSLProgram::createFromSourceStrings
      (
        string ( sheet_vert.readAll().constData() ),
        string ( sheet_geom.readAll().constData() ),
        string ( sheet_frag.readAll().constData() ),GL_TRIANGLES_ADJACENCY_EXT,GL_QUADS
        );

  sheet_vert.close();
  sheet_geom.close();
  sheet_frag.close();

  s_sheetShader->GetProgramLog ( sheet_log );

  _LOG_VAR ( sheet_log );

  //initialize helix shader
  string helix_log;

  QFile helix_vert ( ":/shaders/helix.vert" );
  QFile helix_geom ( ":/shaders/helix.geom" );
  QFile helix_frag ( ":/shaders/helix.frag" );

  helix_vert.open ( QIODevice::ReadOnly );
  helix_geom.open ( QIODevice::ReadOnly );
  helix_frag.open ( QIODevice::ReadOnly );

  s_helixShader = GLSLProgram::createFromSourceStrings
      (
        string ( helix_vert.readAll().constData() ),
        string ( helix_geom.readAll().constData() ),
        string ( helix_frag.readAll().constData() ),GL_LINES_ADJACENCY_EXT,GL_QUADS
        );

  helix_vert.close();
  helix_geom.close();
  helix_frag.close();

  s_helixShader->GetProgramLog ( helix_log );

  _LOG_VAR ( helix_log );



  //initialize tubes shader
  string cyl_log;

  QFile cyl_vert ( "/home/nithin/projects/proteinvis/pv_app/resources/cylinder_adj_vert.glsl" );
  QFile cyl_geom ( "/home/nithin/projects/proteinvis/pv_app/resources/cylinder_adj_geom.glsl" );
  QFile cyl_frag ( "/home/nithin/projects/proteinvis/pv_app/resources/cylinder_adj_frag.glsl" );

  cyl_vert.open ( QIODevice::ReadOnly );
  cyl_geom.open ( QIODevice::ReadOnly );
  cyl_frag.open ( QIODevice::ReadOnly );

  s_tubeShader = GLSLProgram::createFromSourceStrings
      (
        string ( cyl_vert.readAll().constData() ),
        string ( cyl_geom.readAll().constData() ),
        string ( cyl_frag.readAll().constData() ),
        GL_LINES_ADJACENCY,
        GL_TRIANGLE_STRIP
        );

  cyl_vert.close();
  cyl_geom.close();
  cyl_frag.close();

  s_tubeShader->GetProgramLog ( cyl_log );

  _LOG_VAR(cyl_log);
  _LOG_VAR(cyl_log.find("error"))

  if( cyl_log.find("error") != string::npos)
    throw std::runtime_error("failed compiling cylinder shader\n"+cyl_log);

  //initialize tubes shader
  string helix_imposter_log;


  QFile helix_imposter_vert ( ":/shaders/cylinder.vert" );
  QFile helix_imposter_geom ( ":/shaders/glslConvertCylinder.geom" );
  QFile helix_imposter_frag ( ":/shaders/glslHelix3.frag" );

  helix_imposter_vert.open ( QIODevice::ReadOnly );
  helix_imposter_geom.open ( QIODevice::ReadOnly );
  helix_imposter_frag.open ( QIODevice::ReadOnly );

  s_helixImposterShader = GLSLProgram::createFromSourceStrings
      (
        string ( helix_imposter_vert.readAll().constData() ),
        string ( helix_imposter_geom.readAll().constData() ),
        string ( helix_imposter_frag.readAll().constData() ),
        GL_LINES,
        GL_QUADS
        );

  helix_imposter_vert.close();
  helix_imposter_geom.close();
  helix_imposter_frag.close();

  s_helixImposterShader->GetProgramLog ( helix_imposter_log );

  _LOG_VAR ( helix_imposter_log );


}

void secondary_model_t::InitSplines()
{
  const uint *ca_atoms_idx=m_protein->get_ca_atoms_idx();
  num_ca_atoms=m_protein->get_num_ca_atoms();

  const uint *o_atoms_idx=m_protein->get_o_atoms_idx();
  uint num_o_atoms=m_protein->get_num_o_atoms();

  if(num_ca_atoms!=num_o_atoms)
    throw std::runtime_error("invalid input from the parser");

  GetCaAtoms(&ca_atoms_idx);
  GetOAtoms(&o_atoms_idx);
  num_spline_bo_pts=BSplines(ca_atoms,&splineOneControlPts,num_ca_atoms);



  D3DXVECTOR3 prevOxygenNormalized(0,0,0);

  atom_t *shiftedPts=(atom_t *)malloc(sizeof(atom_t)*num_ca_atoms);
  for(int i=0;i<num_ca_atoms;i++)
  {
    shiftedPts[i].bond_end=ca_atoms[i].bond_end;
    shiftedPts[i].bond_start=ca_atoms[i].bond_start;
    shiftedPts[i].radius=ca_atoms[i].radius;
    shiftedPts[i].type_idx=ca_atoms[i].type_idx;

    //shift the position in the direction of Oxygen atoms
    D3DXVECTOR3 oxygenDir(o_atoms[i].x-ca_atoms[i].x,o_atoms[i].y-ca_atoms[i].y,o_atoms[i].z-ca_atoms[i].z);
    D3DXVECTOR3 normalizedOxygenDir(0,0,0);
    Normalize(&oxygenDir,&normalizedOxygenDir);

    float dotProd=D3DXVec3Dot(&prevOxygenNormalized,&normalizedOxygenDir);
    if(dotProd<=0)
    {
      normalizedOxygenDir.x *=-1;
      normalizedOxygenDir.y *=-1;
      normalizedOxygenDir.z *=-1;
    }
    prevOxygenNormalized=normalizedOxygenDir;
    shiftedPts[i].x=ca_atoms[i].x+normalizedOxygenDir.x;
    shiftedPts[i].y=ca_atoms[i].y+normalizedOxygenDir.y;
    shiftedPts[i].z=ca_atoms[i].z+normalizedOxygenDir.z;
  }

  int secondSplineCount=BSplines(shiftedPts,&splineTwoControlPts,num_ca_atoms);

  if(secondSplineCount!=num_spline_bo_pts)
    throw std::runtime_error("second spline count doesn't match the first one...");







  //temp for debugging
  double *coords=new double[num_spline_bo_pts*6];
  int ct=0;
  for(int i=0;i<num_spline_bo_pts;i++)
  {
    coords[ct++]=(double)splineOneControlPts[i].x;
    coords[ct++]=(double)splineOneControlPts[i].y;
    coords[ct++]=(double)splineOneControlPts[i].z;

    coords[ct++]=(double)splineTwoControlPts[i].x;
    coords[ct++]=(double)splineTwoControlPts[i].y;
    coords[ct++]=(double)splineTwoControlPts[i].z;
  }
  m_spline_dir_pts_bo=glutils::buf_obj_t::create_bo(coords,GL_DOUBLE,3,GL_ARRAY_BUFFER,sizeof( GLdouble )*num_spline_bo_pts*6,0);
}




void secondary_model_t::InitSheets(int no,bool isInit)
{
  sheet_indices_t **sheets=m_protein->get_sheets();
  uint num_strands=m_protein->get_num_sheets();


  //get the number of sheets
  map<string,int> sheetIdsMap;
  int count=1;
  for(int i=0;i<num_strands;i++)
  {
    if(sheetIdsMap[sheets[i]->sheetId.c_str()]==0)
      sheetIdsMap[sheets[i]->sheetId.c_str()]=count++;
  }

  count--;
  num_sheets=count;




  //form arrays of cols for each sheet
  int *reds=new int[num_sheets];
  int *greens=new int[num_sheets];
  int *blues=new int[num_sheets];
  for(int i=0;i<count;i++)
  {
    int red=rand()%128;
    int green=rand()%128;
    int blue=rand()%128;
    reds[i]=red;
    greens[i]=green;
    blues[i]=blue;
  }





  list<double> sheetVertsList;
  list<float> idsList;
  list<float> colorList;
  num_sheet_bo_pts=0;
  int sheetId;

  for(int i=0;i<num_strands;i++)
  {
    sheetId=sheetIdsMap[sheets[i]->sheetId.c_str()];
    if(sheetId!=no && no!=0)
      continue;

    std::stringstream out,out2;
    out << sheets[i]->initResidueSeqNo;
    int begin=cAlphaMapping[sheets[i]->initChainId+out.str()];
    out2 << sheets[i]->termResidueSeqNo;
    int end=cAlphaMapping[sheets[i]->termChainId+out2.str()];

    begin=missedMap[begin];
    end=missedMap[end];

    if(isInit)
      indexMap[begin]=end;

    int count=0;
    for(int j=begin;j<=end-3;j++)
    {
      num_sheet_bo_pts+=6;
      //p0
      sheetVertsList.push_back(splineOneControlPts[j].x);
      sheetVertsList.push_back(splineOneControlPts[j].y);
      sheetVertsList.push_back(splineOneControlPts[j].z);

      //p1
      sheetVertsList.push_back(splineOneControlPts[j+1].x);
      sheetVertsList.push_back(splineOneControlPts[j+1].y);
      sheetVertsList.push_back(splineOneControlPts[j+1].z);

      //p2
      sheetVertsList.push_back(splineTwoControlPts[j+1].x-splineOneControlPts[j+1].x);
      sheetVertsList.push_back(splineTwoControlPts[j+1].y-splineOneControlPts[j+1].y);
      sheetVertsList.push_back(splineTwoControlPts[j+1].z-splineOneControlPts[j+1].z);

      //p3
      sheetVertsList.push_back(splineOneControlPts[j+2].x);
      sheetVertsList.push_back(splineOneControlPts[j+2].y);
      sheetVertsList.push_back(splineOneControlPts[j+2].z);

      //p4
      sheetVertsList.push_back(splineTwoControlPts[j+2].x-splineOneControlPts[j+2].x);
      sheetVertsList.push_back(splineTwoControlPts[j+2].y-splineOneControlPts[j+2].y);
      sheetVertsList.push_back(splineTwoControlPts[j+2].z-splineOneControlPts[j+2].z);


      //p5
      sheetVertsList.push_back(splineOneControlPts[j+3].x);
      sheetVertsList.push_back(splineOneControlPts[j+3].y);
      sheetVertsList.push_back(splineOneControlPts[j+3].z);


      for(int k=0;k<6;k++)
      {
        int no=(int)(sheets[i]->sheetId.c_str()[0])-(int)'A';
        colorList.push_back(reds[no]);
        colorList.push_back(greens[no]);
        colorList.push_back(blues[no]);
      }


      if(j==begin)
      {
        for(int k=0;k<6;k++)
          idsList.push_back(0);

        continue;
      }

      if((end-3-j)>=(g_segs_btw_ctrlPts/2))
      {
        for(int k=0;k<6;k++)
          idsList.push_back(-1);
      }
      else
      {
        count++;
        for(int k=0;k<6;k++)
          idsList.push_back(count);

      }

    }


  }

  double *sheetVerts=new double[num_sheet_bo_pts*3];
  copy(sheetVertsList.begin(),sheetVertsList.end(),sheetVerts);
  m_sheet_Pts_bo=glutils::buf_obj_t::create_bo(sheetVerts,GL_DOUBLE,3,GL_ARRAY_BUFFER,sizeof( GLdouble )*num_sheet_bo_pts*3,0);


  float *ids=new float[idsList.size()];
  copy(idsList.begin(),idsList.end(),ids);
  m_sheet_ids_bo=glutils::buf_obj_t::create_bo(ids,GL_FLOAT,1,GL_ARRAY_BUFFER,sizeof( GLfloat )*idsList.size(),0);



  float *colArray=new float[colorList.size()];
  copy(colorList.begin(),colorList.end(),colArray);
  m_sheet_cols_bo=glutils::buf_obj_t::create_bo(colArray,GL_FLOAT,3,GL_ARRAY_BUFFER,sizeof( GLfloat )*colorList.size(),0);



  sheetVertsList.clear();
  idsList.clear();
  colorList.clear();
}

void secondary_model_t::InitHelices(int no,bool isInit)
{
  helix_indices_t **helices=m_protein->get_helices();
  num_helices=m_protein->get_num_helices();



  list<double> helixVertsList;
  list<double> helixImposterVertsList;
  num_helix_bo_pts=0;
  for(int i=0;i<num_helices;i++)
  {

    if(helices[i]->helixSerialNo!=no && no!=0)
      continue;


    std::stringstream out,out2;
    out << helices[i]->initResidueSeqNo;
    int begin=cAlphaMapping[helices[i]->initChainId+out.str()];
    out2 << helices[i]->termResidueSeqNo;
    int end=cAlphaMapping[helices[i]->termChainId+out2.str()];



    int offset=helices[i]->length*g_segs_btw_ctrlPts/10;
    begin=missedMap[begin];
    end=missedMap[end];
    /* begin*=g_segs_btw_ctrlPts;
         end*=g_segs_btw_ctrlPts;

         begin-=beginOff;
         end-=endOff;*/

    //get the helix axis
    D3DXVECTOR3 p1=splineOneControlPts[(end-1-offset)];
    D3DXVECTOR3 p0=splineOneControlPts[(end-2-offset)];
    D3DXVECTOR3 p2=splineOneControlPts[(end-offset)];

    D3DXVECTOR3 v1=D3DSubtract(&p0,&p1);
    D3DXVECTOR3 v2=D3DSubtract(&p2,&p1);

    D3DXVECTOR3 n1(0,0,0);
    D3DXVECTOR3 temp=D3DAdd(&v1,&v2);
    Normalize(&temp,&n1);

    D3DXVECTOR3 centre=D3DAdd(&p1,&n1);


    p1=splineOneControlPts[(begin+1+offset)];
    p0=splineOneControlPts[(begin+2+offset)];
    p2=splineOneControlPts[(begin+offset)];

    v1=D3DSubtract(&p0,&p1);
    v2=D3DSubtract(&p2,&p1);

    D3DXVECTOR3 n2(0,0,0);
    temp=D3DAdd(&v1,&v2);
    Normalize(&temp,&n2);

    D3DXVECTOR3 centre2=D3DAdd(&p1,&n2);


    //D3DXVECTOR3 helixUp=D3DXCross(&n1,&n2);
    D3DXVECTOR3 helixUp=D3DSubtract(&centre2,&centre);



    if(isInit)
      indexMap[begin]=end;


    helixImposterVertsList.push_back(splineOneControlPts[begin].x);
    helixImposterVertsList.push_back(splineOneControlPts[begin].y);
    helixImposterVertsList.push_back(splineOneControlPts[begin].z);

    helixImposterVertsList.push_back(splineOneControlPts[end].x);
    helixImposterVertsList.push_back(splineOneControlPts[end].y);
    helixImposterVertsList.push_back(splineOneControlPts[end].z);

    for(int j=begin;j<=end-2;j++)
    {

      num_helix_bo_pts+=4;
      //p0
      helixVertsList.push_back(splineOneControlPts[j].x);
      helixVertsList.push_back(splineOneControlPts[j].y);
      helixVertsList.push_back(splineOneControlPts[j].z);

      //p1
      helixVertsList.push_back(splineOneControlPts[j+1].x);
      helixVertsList.push_back(splineOneControlPts[j+1].y);
      helixVertsList.push_back(splineOneControlPts[j+1].z);

      //p2
      helixVertsList.push_back(splineOneControlPts[j+2].x);
      helixVertsList.push_back(splineOneControlPts[j+2].y);
      helixVertsList.push_back(splineOneControlPts[j+2].z);

      //axis
      helixVertsList.push_back(helixUp.x);
      helixVertsList.push_back(helixUp.y);
      helixVertsList.push_back(helixUp.z);

    }

  }


  double *helixVerts=new double[num_helix_bo_pts*3];
  copy(helixVertsList.begin(),helixVertsList.end(),helixVerts);
  m_helix_Pts_bo=glutils::buf_obj_t::create_bo(helixVerts,GL_DOUBLE,3,GL_ARRAY_BUFFER,sizeof( GLdouble )*num_helix_bo_pts*3,0);

  num_helix_imposter_bo_pts=helixImposterVertsList.size()/3;
  double *helixImposterVerts=new double[helixImposterVertsList.size()];
  copy(helixImposterVertsList.begin(),helixImposterVertsList.end(),helixImposterVerts);
  m_helix_imposter_bo=glutils::buf_obj_t::create_bo(helixImposterVerts,GL_DOUBLE,3,GL_ARRAY_BUFFER,sizeof( GLdouble )*helixImposterVertsList.size(),0);

  helixVertsList.clear();

}


void secondary_model_t::InitTubes()
{
  vertex_list_t pts;

  for(int i=0;i<num_spline_bo_pts;i++)
  {
    pts.push_back(vertex_t(splineOneControlPts[i].x,
                           splineOneControlPts[i].y,
                           splineOneControlPts[i].z));

  }
  m_spline_pts_bo = make_buf_obj(pts);
}

void secondary_model_t::Render()
{
  glPushAttrib(GL_ENABLE_BIT);
  glDisable(GL_LIGHTING);
  glColor3ub ( 120, 0, 0 );
  m_spline_pts_bo->bind_to_vertex_pointer();
  glDrawArrays(GL_LINES, 0, num_spline_bo_pts-1);
  m_spline_pts_bo->unbind_from_vertex_pointer();
  glPopAttrib();


  /*glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_LIGHTING);
    glColor3ub ( 0, 120, 0 );
    m_spline_dir_pts_bo->bind_to_vertex_pointer();
    glDrawArrays(GL_LINES, 0, num_spline_bo_pts*2);
    m_spline_dir_pts_bo->unbind_from_vertex_pointer();
    glPopAttrib();*/

  RenderTubes();
  RenderSheets();
  RenderHelices();
}


void secondary_model_t::RenderSheets()
{
  if(num_sheet_bo_pts<0)
    return;


  glPushAttrib ( GL_ENABLE_BIT );
  s_sheetShader->use();
  GLuint indices_attrib = s_sheetShader->getAttributeLocation ( "indices" );
  m_sheet_ids_bo->bind_to_vertex_attrib_pointer(indices_attrib);
  GLuint pos_attrib = s_sheetShader->getAttributeLocation ( "position" );
  m_sheet_Pts_bo->bind_to_vertex_attrib_pointer(pos_attrib);
  GLuint cols_attrib = s_sheetShader->getAttributeLocation ( "colorsInp" );
  m_sheet_cols_bo->bind_to_vertex_attrib_pointer(cols_attrib);
  glDrawArrays(GL_TRIANGLES_ADJACENCY_EXT,0,num_sheet_bo_pts);
  m_sheet_Pts_bo->unbind_from_vertex_attrib_pointer(pos_attrib);
  m_sheet_ids_bo->unbind_from_vertex_attrib_pointer(indices_attrib);
  m_sheet_cols_bo->unbind_from_vertex_attrib_pointer(cols_attrib);
  s_sheetShader->disable();
  glPopAttrib();

}



void secondary_model_t::RenderHelices()
{
  if(num_helix_bo_pts<0)
    return;

  glPushAttrib ( GL_ENABLE_BIT );
  s_helixShader->use();
  m_helix_Pts_bo->bind_to_vertex_pointer();
  glDrawArrays(GL_LINES_ADJACENCY_EXT,0,num_helix_bo_pts);
  m_helix_Pts_bo->unbind_from_vertex_pointer();
  s_helixShader->disable();
  glPopAttrib();
}


void secondary_model_t::RenderTubes()
{
  glPushAttrib ( GL_ENABLE_BIT );
  glDisable(GL_LIGHTING);
  s_tubeShader->use();
  glColor3ub ( 120, 60, 18 );
  m_spline_pts_bo->bind_to_vertex_pointer();
  glDrawArrays(GL_LINE_STRIP_ADJACENCY,0,num_spline_bo_pts);
  m_spline_pts_bo->unbind_from_vertex_pointer();
  s_tubeShader->disable();
  glPopAttrib();
}

void secondary_model_t::RenderImposterHelices()
{
  //glDisable(GL_CULL_FACE);
  glPushAttrib ( GL_ENABLE_BIT );
  s_helixImposterShader->use();
  s_helixImposterShader->sendUniform ( "ug_cylinder_radius", ( float ) 1 );
  //s_tubeShader->sendUniform ( "offset",offset );
  m_helix_imposter_bo->bind_to_vertex_pointer();
  glDrawArrays(GL_LINES,0,num_helix_imposter_bo_pts);
  m_helix_imposter_bo->unbind_from_vertex_pointer();
  s_helixImposterShader->disable();
  glPopAttrib();
}

