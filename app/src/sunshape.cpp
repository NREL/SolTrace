
/*******************************************************************************************************
*  Copyright 2018 Alliance for Sustainable Energy, LLC
*
*  NOTICE: This software was developed at least in part by Alliance for Sustainable Energy, LLC
*  ("Alliance") under Contract No. DE-AC36-08GO28308 with the U.S. Department of Energy and the U.S.
*  The Government retains for itself and others acting on its behalf a nonexclusive, paid-up,
*  irrevocable worldwide license in the software to reproduce, prepare derivative works, distribute
*  copies to the public, perform publicly and display publicly, and to permit others to do so.
*
*  Redistribution and use in source and binary forms, with or without modification, are permitted
*  provided that the following conditions are met:
*
*  1. Redistributions of source code must retain the above copyright notice, the above government
*  rights notice, this list of conditions and the following disclaimer.
*
*  2. Redistributions in binary form must reproduce the above copyright notice, the above government
*  rights notice, this list of conditions and the following disclaimer in the documentation and/or
*  other materials provided with the distribution.
*
*  3. The entire corresponding source code of any redistribution, with or without modification, by a
*  research entity, including but not limited to any contracting manager/operator of a United States
*  National Laboratory, any institution of higher learning, and any non-profit organization, must be
*  made publicly available under this license for as long as the redistribution is made available by
*  the research entity.
*
*  4. Redistribution of this software, without modification, must refer to the software by the same
*  designation. Redistribution of a modified version of this software (i) may not refer to the modified
*  version by the same designation, or by any confusingly similar designation, and (ii) must refer to
*  the underlying software originally provided by Alliance as "SolTrace". Except to comply with the 
*  foregoing, the term "SolTrace", or any confusingly similar designation may not be used to refer to 
*  any modified version of this software or any modified version of the underlying software originally 
*  provided by Alliance without the prior written consent of Alliance.
*
*  5. The name of the copyright holder, contributors, the United States Government, the United States
*  Department of Energy, or any of their employees may not be used to endorse or promote products
*  derived from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
*  IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
*  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER,
*  CONTRIBUTORS, UNITED STATES GOVERNMENT OR UNITED STATES DEPARTMENT OF ENERGY, NOR ANY OF THEIR
*  EMPLOYEES, BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
*  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
*  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
*  IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
*  THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*******************************************************************************************************/


#include <wx/stattext.h>
#include <wx/checkbox.h>
#include <wx/sizer.h>
#include <wx/button.h>
#include <wx/statbox.h>
#include <wx/filedlg.h>
#include <wx/msgdlg.h>

#include <wex/csv.h>
#include <wex/extgrid.h>
#include <wex/numeric.h>
#include <wex/plot/plplotctrl.h>
#include <wex/plot/pllineplot.h>

#include <wex/extgrid.h>
#include <wex/radiochoice.h>

#include "sunshape.h"
#include "soltrace.h"

enum { _id_first = wxID_HIGHEST+123,
	ID_X, ID_Y, ID_Z, ID_SIGMA, ID_HALFWIDTH, ID_POINTSRC,
	ID_SHAPE, ID_XYZLDH, ID_NPOINTS, ID_USERDATA,
	ID_COPY, ID_PASTE, ID_LOAD, ID_SAVE, ID_PREDEFINED };

BEGIN_EVENT_TABLE( SunShapeForm, wxPanel )
	EVT_NUMERIC( ID_X, SunShapeForm::OnCommand )
	EVT_NUMERIC( ID_Y, SunShapeForm::OnCommand )
	EVT_NUMERIC( ID_Z, SunShapeForm::OnCommand )
	EVT_NUMERIC( ID_SIGMA, SunShapeForm::OnCommand )
	EVT_NUMERIC( ID_HALFWIDTH, SunShapeForm::OnCommand )
	EVT_CHECKBOX( ID_POINTSRC, SunShapeForm::OnCommand )
	EVT_RADIOBUTTON( ID_SHAPE, SunShapeForm::OnCommand )
	EVT_RADIOBUTTON( ID_XYZLDH, SunShapeForm::OnCommand )
	EVT_NUMERIC( ID_NPOINTS, SunShapeForm::OnCommand )
	EVT_GRID_CMD_CELL_CHANGED( ID_USERDATA, SunShapeForm::OnGridCellChange )
	EVT_BUTTON( ID_COPY, SunShapeForm::OnCommand )
	EVT_BUTTON( ID_PASTE, SunShapeForm::OnCommand )
	EVT_BUTTON( ID_LOAD, SunShapeForm::OnCommand )
	EVT_BUTTON( ID_SAVE, SunShapeForm::OnCommand )
	EVT_BUTTON( ID_PREDEFINED, SunShapeForm::OnCommand )
END_EVENT_TABLE()

SunShapeForm::SunShapeForm( wxWindow *parent, Project &prj )
	: wxPanel( parent ), m_prj(prj)
{
	m_pointSrc = new wxCheckBox( this, ID_POINTSRC, "Point source at finite distance" );
	
	wxStaticBoxSizer *pos_sizer = new wxStaticBoxSizer( wxVERTICAL, this, "Sun position" );
	
	m_xyzOrLdh = new wxRadioChoice( pos_sizer->GetStaticBox(), ID_XYZLDH );
	m_xyzOrLdh->Add( "Global coordinates" );
	m_xyzOrLdh->Add( "Latitude, day, hour" );

	m_lblX = new wxStaticText( pos_sizer->GetStaticBox(), wxID_ANY, "X:" );
	m_lblY = new wxStaticText( pos_sizer->GetStaticBox(), wxID_ANY, "Y:" );
	m_lblZ = new wxStaticText( pos_sizer->GetStaticBox(), wxID_ANY, "Z:" );

	m_x = new wxNumericCtrl( pos_sizer->GetStaticBox(), ID_X );
	m_y = new wxNumericCtrl( pos_sizer->GetStaticBox(), ID_Y );
	m_z = new wxNumericCtrl( pos_sizer->GetStaticBox(), ID_Z );

	wxGridSizer *xyz_sizer = new wxGridSizer( 3, 3, 3 );
	xyz_sizer->Add( m_lblX );
	xyz_sizer->Add( m_lblY );
	xyz_sizer->Add( m_lblZ );
	xyz_sizer->Add( m_x );
	xyz_sizer->Add( m_y );
	xyz_sizer->Add( m_z );
	
	wxBoxSizer *box_sizer1 = new wxBoxSizer( wxHORIZONTAL );
	box_sizer1->Add( m_xyzOrLdh, 0, wxLEFT|wxRIGHT|wxEXPAND, 10 );
	box_sizer1->Add( xyz_sizer, 0, wxLEFT|wxEXPAND, 20 );

	pos_sizer->Add( box_sizer1, 0, wxALL|wxEXPAND, 0 );
	m_dirLabel = new wxStaticText( pos_sizer->GetStaticBox(), wxID_ANY, "Note: global X axis points west, Y to zenith, Z to north.");
	pos_sizer->Add( m_dirLabel, 0, wxALL, 8 );


	wxStaticBoxSizer *shp_sizer = new wxStaticBoxSizer( wxHORIZONTAL, this, "Sun shape" );
	m_grpSunShape = shp_sizer->GetStaticBox();

	m_shape = new wxRadioChoice( m_grpSunShape, ID_SHAPE );
	m_shape->LayoutEvenly( true );
	m_shape->Add( "Gaussian" );
	m_shape->Add( "Pillbox" );
	m_shape->Add( "User defined, # of points:" );
	
	m_sigma = new wxNumericCtrl( m_grpSunShape, ID_SIGMA );
	m_halfWidth = new wxNumericCtrl( m_grpSunShape, ID_HALFWIDTH );
	m_userPoints = new wxNumericCtrl( m_grpSunShape, ID_NPOINTS, wxNUMERIC_INTEGER );

	wxBoxSizer *box2 = new wxBoxSizer( wxVERTICAL );
	box2->Add( m_sigma, 0, wxALL|wxEXPAND, 2 );
	box2->Add( m_halfWidth, 0, wxALL|wxEXPAND, 2 );
	box2->Add( m_userPoints, 0, wxALL|wxEXPAND, 2 );

	wxBoxSizer *box1 = new wxBoxSizer( wxHORIZONTAL );
	box1->Add( m_shape, 0, wxALL|wxEXPAND, 10 );
	box1->Add( box2, 0, wxALL, 10 );

	wxBoxSizer *btn_sizer = new wxBoxSizer( wxHORIZONTAL );
	btn_sizer->Add( m_predefined = new wxButton( m_grpSunShape, ID_PREDEFINED, "Predefined...", wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT ) );
	btn_sizer->Add( m_load = new wxButton( m_grpSunShape, ID_LOAD, "Load", wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT ) );
	btn_sizer->Add( m_save = new wxButton( m_grpSunShape, ID_SAVE, "Save", wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT ) );
	btn_sizer->Add( m_copy = new wxButton( m_grpSunShape, ID_COPY, "Copy", wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT ) );
	btn_sizer->Add( m_paste = new wxButton( m_grpSunShape, ID_PASTE, "Paste", wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT ) );
	
	m_userData = new wxExtGridCtrl( m_grpSunShape, ID_USERDATA );
	m_userData->CreateGrid( 1, 2 );	
	m_userData->SetColLabelValue(0, "Angle (mrad)");
	m_userData->SetColLabelValue(1, "Intensity" );

	wxBoxSizer *left_sizer = new wxBoxSizer( wxVERTICAL );
	left_sizer->Add( box1 );
	left_sizer->Add( btn_sizer, 0, wxTOP, 2 );
	left_sizer->Add( m_userData, 0, wxTOP|wxBOTTOM|wxEXPAND, 5 );

	m_plot = new wxPLPlotCtrl( m_grpSunShape, wxID_ANY, wxDefaultPosition, wxSize(400,300) );

	shp_sizer->Add( left_sizer );
	shp_sizer->Add( m_plot );
	
	wxBoxSizer *top_sizer = new wxBoxSizer( wxVERTICAL );
	top_sizer->Add( m_pointSrc, 0, wxALL, 10 );
	top_sizer->Add( pos_sizer, 0, wxALL, 10 );
	top_sizer->Add( shp_sizer, 0, wxALL, 10 );
	SetSizer( top_sizer );

	UpdateFromData();
}

void SunShapeForm::UpdateFromData()
{
	SunShape &S = m_prj.Sun;

	switch(S.Shape)
	{
	case SunShape::GAUSSIAN:
		m_shape->SetSelection(0);
		break;
	case SunShape::PILLBOX:
		m_shape->SetSelection(1);
		break;
	case SunShape::USER_DEFINED:
		m_shape->SetSelection(2);
		break;
	}

	UpdateWidgetsEnabled();
	UpdatePlot();

	m_sigma->SetValue( S.Sigma );
	m_halfWidth->SetValue( S.HalfWidth );

	m_userPoints->SetValue( S.UserShapeData.size() );
	UpdateTable();

	m_xyzOrLdh->SetSelection( S.UseLDHSpec == false || S.PointSource ? 0 : 1 );
	m_xyzOrLdh->Enable( 1, !S.PointSource );

	UpdateCoordinates();

	m_pointSrc->SetValue( S.PointSource );
	UpdatePointSourceWidgets();
}

void SunShapeForm::UpdateCoordinates()
{
	// save the values to set because
	// modifying the ranges of the inputs
	// may cause a signal to be emitted
	// and thus overwrite the values in the Sun structure

	if ( m_xyzOrLdh->GetSelection() == 0 )
	{
		double x = m_prj.Sun.X;
		double y = m_prj.Sun.Y;
		double z = m_prj.Sun.Z;

		m_dirLabel->Show( true );
		m_lblX->SetLabel("X");
		m_lblY->SetLabel("Y");
		m_lblZ->SetLabel("Z");
				
		m_x->ClearRange();
		m_y->ClearRange();
		m_z->ClearRange();

		m_x->SetValue( x );
		m_y->SetValue( y );
		m_z->SetValue( z );

	}
	else
	{
		double l = m_prj.Sun.Latitude;
		double d = m_prj.Sun.Day;
		double h = m_prj.Sun.Hour;

		m_dirLabel->Show( true );
		m_lblX->SetLabel("Latitude");
		m_lblY->SetLabel("Day");
		m_lblZ->SetLabel("Hour");

		m_x->SetRange( -90, 90 );
		m_y->SetRange( 1, 365 );
		m_z->SetRange( 0, 24 );

		m_x->SetValue( l );
		m_y->SetValue( d );
		m_z->SetValue( h );
	}

	Layout();
}

void SunShapeForm::UpdatePointSourceWidgets()
{
	bool en = !m_pointSrc->GetValue();
	m_grpSunShape->Show(en);
	m_plot->Show(en);
	Layout();
}

void SunShapeForm::UpdateWidgetsEnabled()
{
	int shp = m_shape->GetSelection();

	m_sigma->Enable( shp==0 );
	m_halfWidth->Enable( shp==1 );

	bool bUserDefined = shp==2;

	m_userPoints->Enable( bUserDefined );
	m_userData->Enable( bUserDefined );
	m_predefined->Enable( bUserDefined );
	m_load->Enable( bUserDefined );
	m_save->Enable( bUserDefined );
	m_copy->Enable( bUserDefined );
	m_paste->Enable( bUserDefined );
}

void SunShapeForm::UpdatePlot()
{
	m_plot->DeleteAllPlots();

	if ( m_pointSrc->GetValue() )
	{
		m_plot->SetTitle( wxEmptyString );
		m_plot->Refresh();
		return;
	}


	if ( 0 == m_shape->GetSelection() )
	{
		double sigma = m_prj.Sun.Sigma;
		int npoints = 50;
		double thetax = -sigma*3;
		double thetainc = 6*sigma/npoints;

		std::vector<wxRealPoint> data;

		for (int i=0;i<npoints;i++)
		{
			data.push_back( wxRealPoint( thetax, 1.0/exp(thetax*thetax/(2*sigma*sigma)) ) );
			thetax += thetainc;
		}

		m_plot->AddPlot( new wxPLLinePlot( data, "Sun shape", *wxRED ) );
		m_plot->GetXAxis1()->SetWorld( -3.3*sigma, 3.3*sigma );
		m_plot->GetYAxis1()->SetWorld(  0, 1.2 );
	}
	else if ( 1 == m_shape->GetSelection() )
	{
		double hw = m_prj.Sun.HalfWidth;
		std::vector<wxRealPoint> data;
		data.push_back( wxRealPoint( -hw, 0 ) );
		data.push_back( wxRealPoint( -hw, 1 ) );
		data.push_back( wxRealPoint( hw, 1 ) );
		data.push_back( wxRealPoint( hw, 0 ) );

		m_plot->AddPlot( new wxPLLinePlot( data, "Sun shape", *wxRED ) );
		m_plot->GetXAxis1()->SetWorld( -3*hw, 3*hw );
		m_plot->GetYAxis1()->SetWorld(  0, 1.2 );
	}
	else if ( 2 == m_shape->GetSelection() )
	{
		SunShape &S = m_prj.Sun;
		int i;
		int npoints = S.UserShapeData.size();
		if (npoints >= 2)
		{
			std::vector<wxRealPoint> data;

			data.resize( 2*npoints - 1 );

			double maxx=0, maxy=0;
			for(i=0;i<npoints;i++)
			{
				data[npoints+i-1] = wxRealPoint( S.UserShapeData[i].x, S.UserShapeData[i].y );
				if (data[npoints+i-1].x > maxx) maxx = data[npoints+i-1].x;
				if (data[npoints+i-1].y > maxy) maxy = data[npoints+i-1].y;
			}
			for (i=0;i<npoints-1;i++)
			{
				data[i] = wxRealPoint( -S.UserShapeData[npoints-i-1].x, S.UserShapeData[npoints-i-1].y );
			}
			
			m_plot->AddPlot( new wxPLLinePlot( data, "Sun shape", *wxRED ) );
			m_plot->GetXAxis1()->SetWorld( -1.3*maxx, 1.3*maxx  );
			m_plot->GetYAxis1()->SetWorld(  0, 1.2 );
		}

	}


	if ( m_plot->GetXAxis1() ) m_plot->GetXAxis1()->SetLabel( "Angle from center (mrad)" );
	if ( m_plot->GetYAxis1() ) m_plot->GetYAxis1()->SetLabel( "Intensity" );
	m_plot->SetTitle( "Sun shape profile" );
	m_plot->ShowLegend( false );
	m_plot->Refresh();
}

void SunShapeForm::UpdateTable()
{
	int npoints = m_userPoints->AsInteger();
	if ( npoints <= 0 )
		return;

	m_userData->ResizeGrid( npoints, 2 );

	SunShape &S = m_prj.Sun;

	int curlen = S.UserShapeData.size();

	if (curlen != npoints)
	{
		S.UserShapeData.resize( npoints );
		for (int i=curlen;i<npoints;i++)
		{
			S.UserShapeData[i].x = 0.0;
			S.UserShapeData[i].y = 0.0;
		}

		Modified();
	}

	for (int r=0;r<npoints;r++)
	{
		m_userData->SetCellValue( r, 0, wxString::Format("%lg", S.UserShapeData[r].x) );
		m_userData->SetCellValue( r, 1, wxString::Format("%lg", S.UserShapeData[r].y) );
	}

	Layout();
	UpdatePlot();
}



void SunShapeForm::OnCommand( wxCommandEvent &evt )
{
	SunShape &S = m_prj.Sun;

	switch( evt.GetId() )
	{
	case ID_POINTSRC:
		S.PointSource = m_pointSrc->GetValue();
		if ( S.PointSource )
		{
			m_xyzOrLdh->SetLabel( 0, "Global coordinates (Note: point source cannot be at origin 0,0,0)");
			if ( S.UseLDHSpec )
			{
				S.UseLDHSpec = false;
				m_xyzOrLdh->SetSelection( 0 );
				UpdateCoordinates();
			}
			
			m_xyzOrLdh->Enable( 1, false );
		}
		else
		{
			m_xyzOrLdh->SetLabel( 0, "Global coordinates" );
			m_xyzOrLdh->Enable( 1, true );
		}

		UpdatePointSourceWidgets();
		Modified();

		break;

	case ID_XYZLDH:		
		S.UseLDHSpec = (m_xyzOrLdh->GetSelection() == 1);
		UpdateCoordinates();
		Modified();
		break;

	case ID_X:
		if ( !S.UseLDHSpec )
			S.X = m_x->Value();
		else
			S.Latitude = m_x->Value();
		Modified();
		break;

	case ID_Y:
		if ( !S.UseLDHSpec )
			S.Y = m_y->Value();
		else
			S.Day = m_y->Value();
		Modified();
		break;

	case ID_Z:
		if ( !S.UseLDHSpec )
			S.Z = m_z->Value();
		else
			S.Hour = m_z->Value();
		Modified();
		break;

	case ID_SHAPE:
		if ( 0 == m_shape->GetSelection() )
			S.Shape = SunShape::GAUSSIAN;
		else if ( 1 == m_shape->GetSelection() )
			S.Shape = SunShape::PILLBOX;
		else
			S.Shape = SunShape::USER_DEFINED;

		UpdateWidgetsEnabled();
		UpdatePlot();
		Modified();
		break;

	case ID_SIGMA:
		S.Sigma = m_sigma->Value();
		UpdatePlot();
		Modified();
		break;

	case ID_HALFWIDTH:
		S.HalfWidth = m_halfWidth->Value();
		UpdatePlot();
		Modified();
		break;

	case ID_NPOINTS:
		UpdateTable();
		break;

	case ID_PREDEFINED:
		wxMessageBox("Currently there are no predefined sun shape files available." );
		break;

	case ID_LOAD:
	{
		wxFileDialog dlg( this, "Load sun shape data", wxEmptyString, wxEmptyString, "Comma-separated values (*.csv)|*.csv", wxFD_OPEN );
		if ( wxID_OK == dlg.ShowModal() )
		{
			wxCSVData csv;
			if ( csv.ReadFile( dlg.GetPath() ) && csv.NumRows() > 0 )
			{
				int nr = csv.NumRows();
				m_userPoints->SetValue( nr );
				S.UserShapeData.resize( nr );
				for( int i=0;i<nr;i++ )
					S.UserShapeData[i] = PointF( wxAtof( csv(i,0) ), wxAtof( csv(i,1) ) );
				UpdateTable();
				Modified();
			}
			else				
				wxMessageBox("Failed to open input file " + dlg.GetPath() + " for reading." );
		}
	}
		break;
	case ID_SAVE:
	{
		wxFileDialog dlg( this, "Save sun shape data", wxEmptyString, wxEmptyString, "Comma-separated values (*.csv)|*.csv", wxFD_SAVE|wxFD_OVERWRITE_PROMPT );
		if ( wxID_OK == dlg.ShowModal() )
		{
			if ( FILE *fp = fopen( dlg.GetPath().c_str(), "w" ) )
			{
				for( size_t i=0;i<S.UserShapeData.size();i++ )
					fprintf( fp, "%lg, %lg\n", S.UserShapeData[i].x, S.UserShapeData[i].y );

				fclose(fp);
			}
			else
				wxMessageBox("Failed to open output file " + dlg.GetPath() + " for writing." );
		}
	}
		break;

	case ID_COPY:
		m_userData->Copy();
		break;

	case ID_PASTE:
		m_userData->Paste( wxExtGridCtrl::PASTE_ALL_RESIZE_ROWS );
		break;
	}
}

void SunShapeForm::OnGridCellChange( wxGridEvent &evt )
{
	SunShape &S = m_prj.Sun;
	int row = evt.GetRow();
	int col = evt.GetCol();

	if ( row < 0 || col < 0 )
	{
		int nrows = m_userData->GetNumberRows();
		if ( nrows <= 0 ) return;

		m_userPoints->SetValue( nrows );
		S.UserShapeData.resize( nrows );
		for( int i=0;i<nrows;i++ )
		{
			S.UserShapeData[i].x = wxAtof( m_userData->GetCellValue( i, 0 ) );
			S.UserShapeData[i].y = wxAtof( m_userData->GetCellValue( i, 1 ) );
		}

		UpdateTable();
	}
	else
	{
		double val = wxAtof( m_userData->GetCellValue( row, col ) );
		m_userData->SetCellValue( row, col, wxString::Format( "%lg", val ) );

		if ( row < S.UserShapeData.size() )
		{
			if ( col == 0 )
				S.UserShapeData[row].x = val;
			else
				S.UserShapeData[row].y = val;

			Modified();
			UpdatePlot();
		}
	}
}

void SunShapeForm::Modified() 
{
	MainWindow::Instance().SetModified( true );
}
