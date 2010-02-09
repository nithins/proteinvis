/***************************************************************************
 *   Copyright (C) 2009 by Nithin Shivashankar,   *
 *   nithin@gauss   *
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

#include <string>
#include <climits>
#include <vector>

#include <boost/regex.hpp>

#include <logutil.h>
#include <cpputils.h>

#include <framework.h>

#include <proteinModel.h>


using namespace std;

static const boost::regex protein_file_re ( "(-pf ([[:alnum:]\\./_]+))" );
static const std::string  protein_file_replace ( "(())" );

static const boost::regex surface_file_re ( "(-sf ([[:alnum:]\\./_]+))" );
static const std::string  surface_file_replace ( "(())" );

static const boost::regex alpha_complex_file_re ( "(-acf ([[:alnum:]\\./_]+))" );
static const std::string  alpha_complex_file_replace ( "(())" );

static const boost::regex pocket_file_re ( "(-pocf ([[:alnum:]\\./_]+))" );
static const std::string  pocket_file_replace ( "(())" );

static const boost::regex tet_file_re ( "(-tetf ([[:alnum:]\\./_]+))" );
static const std::string  tet_file_replace ( "(())" );


int main ( int argc, char *argv[] )
{

  try
  {

    if ( argc == 1 )
    {
      cout
          << "usage: " << argv[0]
          << " -pf <crd/pdb file> [-sf <off file>] [-acf <acf file>] [-tetf <tet file> -pocf <poc file>]"
          << endl ;
      exit ( 0 );
    }

    string cmdline ( argv[1] );

    for ( uint i = 2 ; i < ( uint ) argc;i++ )
    {
      cmdline += " ";
      cmdline += argv[i];
    }

    vector<string> protein_def_str_vec;

    split_string ( cmdline, protein_def_str_vec, "-pf" );

    vector<boost::shared_ptr<ProteinModel> > protein_models;

    boost::shared_ptr<IFramework> framework =  IFramework::Create (argc,argv);

    for ( uint i = 0 ; i < protein_def_str_vec.size(); ++i )
    {
      string protein_def_str = string ( "-pf" ) + protein_def_str_vec[i];

      string protein_filename, surface_filename, alpha_complex_filename;

      string pocket_filename, tet_filename;

      boost::smatch matches;

      if ( regex_search ( protein_def_str, matches, protein_file_re ) )
      {
        protein_filename.assign ( matches[2].first, matches[2].second );
      }
      else
      {
        _ERROR ( "no protein file name specified" );
        exit ( 0 );
      }

      if ( regex_search ( protein_def_str, matches, surface_file_re ) )
      {
        surface_filename.assign ( matches[2].first, matches[2].second );
      }

      if ( regex_search ( protein_def_str, matches, alpha_complex_file_re ) )
      {
        alpha_complex_filename.assign ( matches[2].first, matches[2].second );
      }

      if ( regex_search ( protein_def_str, matches, pocket_file_re ) )
      {
        pocket_filename.assign ( matches[2].first, matches[2].second );
      }

      if ( regex_search ( protein_def_str, matches, tet_file_re ) )
      {
        tet_filename.assign ( matches[2].first, matches[2].second );
      }

      if ( tet_filename.size() != 0 && pocket_filename.size() == 0 )
      {
        _LOG ( "tet file specified without pockets .. tet file will not be used" );
        tet_filename.clear();
      }

      if ( tet_filename.size() == 0 && pocket_filename.size() != 0 )
      {
        _LOG ( "pocket file specified without tets .. pocket file will not be used" );
        pocket_filename.clear();
      }


        boost::shared_ptr<ProteinModel> pm(
              new ProteinModel
              ( protein_filename,
                surface_filename ,
                alpha_complex_filename,
                pocket_filename,
                tet_filename )
              );


      protein_models.push_back ( pm );
    }

    for ( uint i = 0 ; i < protein_models.size();++i )
    {
      framework->AddModel ( protein_models[i] );
    }

    framework->Exec();

  }
  catch(const char *str)
  {
    _ERROR(str);
  }

  catch(string str)
  {
    _ERROR(str);
  }

  catch(std::exception e)
  {
    _LOG(e.what());
  }

  return 0;

}
