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

#ifndef TESTMODEL_H_INCLUDED
#define TESTMODEL_H_INCLUDED

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

class test_model_t
{
private:
    //vbos to hold the structures
    glutils::bufobj_ptr_t m_tube_Pts_bo; //vb holds all points on the spline

    GLSLProgram *s_tubeShader;

public:
    test_model_t();
    void Render();

private:
    void InitShaders();

};
#endif
