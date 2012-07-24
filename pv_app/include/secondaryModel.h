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

#ifndef SECONDARYMODEL_H_INCLUDED
#define SECONDARYMODEL_H_INCLUDED

#include <vector>
#include <list>
#include <glutils.h>
#include <DxUtils.h>
#include <protein.h>
#include <oneLevelModel.h>
#include <GLSLProgram.h>
#include <QFile>
#include <QDir>
#include <QResource>
#include <logutil.h>
#include <map>

using namespace::std;

class protein_rd_t;
class GLSLProgram;

const int g_segs_btw_ctrlPts = 30;

class secondary_model_t
{
private:
    //reference to the protein
    boost::shared_ptr<protein_t> m_protein;

    //reference to the list of atoms in the protein
    const atom_t      *atoms;
    uint     num_atoms;

    //all c-alpha atoms in the protein
    atom_t *ca_atoms;
    atom_t *o_atoms;
    uint num_ca_atoms;

    //spline control points
    D3DXVECTOR3 *splineOneControlPts;
    D3DXVECTOR3 *splineTwoControlPts;


    //vbos to hold the structures
    glutils::bufobj_ptr_t m_spline_pts_bo; //vb holds all points on the spline
    glutils::bufobj_ptr_t m_sheet_Pts_bo;//vb indices of all the sheets
    glutils::bufobj_ptr_t m_helix_Pts_bo;//vb indices of all the helices


    glutils::bufobj_ptr_t m_sheet_cols_bo;//colors for sheets
    glutils::bufobj_ptr_t m_sheet_ids_bo;//to know if it is an arrowhead


    int num_spline_bo_pts;
    int num_sheet_bo_pts;
    int num_helix_bo_pts;

    //shader parameters
    GLSLProgram *s_sheetShader;
    GLSLProgram *s_tubeShader;
    GLSLProgram *s_helixShader;
    GLSLProgram *s_helixImposterShader;

    //other data structures
    map<string,uint> cAlphaMapping;//mapping to get the correct point on spline using the pdb data
    map<int,int> indexMap;//used to get the tube indices
    map<int,int> tubeMap;//used to get the tube indices
    map<int,int> missedMap;

    //variables for spin box
    uint num_sheets;
    uint num_helices;

public:
    secondary_model_t(boost::shared_ptr<protein_t> protein_reference);
    ~secondary_model_t();

    //method to render the secondary protein structure
    void RenderSheets();
    void RenderTubes();
    void RenderHelices();
    void RenderImposterHelices();
    //methods to fill buffers with only some elements as required by spinboxes
    void InitHelices(int no,bool isInit);
    void InitSheets(int no,bool isInit);

    void Render();

private:
    //methods to generate the spline on the cpu
    int BSplines(atom_t *atomPts,D3DXVECTOR3 **detailedRef,int atomPtsLength);
    int DetailPtGen(D3DXVECTOR3 *point1,D3DXVECTOR3 *point2,D3DXVECTOR3 *point3,D3DXVECTOR3 *point4,D3DXVECTOR3 **retResult);
    D3DXVECTOR3 Interpolate(D3DXVECTOR3 *point1,D3DXVECTOR3 *point2,D3DXVECTOR3 *point3,D3DXVECTOR3 *point4,float amount);

    //other helper methods
    void GetCaAtoms(const uint **ca_atoms_idx_ref);
    void GetOAtoms(const uint **o_atoms_idx_ref);

    //other initialization methods
    void InitShaders();
    void InitTubes();
    void InitSplines();



    //temp vars for debugging
    glutils::bufobj_ptr_t m_spline_dir_pts_bo;
    glutils::bufobj_ptr_t m_helix_imposter_bo;
    int num_helix_imposter_bo_pts;

};
#endif
