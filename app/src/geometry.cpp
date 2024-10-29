
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


#include <wx/busyinfo.h>
#include <wx/clipbrd.h>
#include <wx/tokenzr.h>
#include <wx/statline.h>
#include <wx/splitter.h>
#include <wx/grid.h>
#include <wx/notebook.h>
#include <wx/filename.h>
#include <wx/tokenzr.h>
#include <wx/stdpaths.h>
#include <wx/gbsizer.h>

#include <wex/extgrid.h>
#include <wex/exttext.h>
#include <wex/numeric.h>
#include <wex/utils.h>

#include "soltrace.h"
#include "geometry.h"


struct TypeInfo
{
	char ap_t;
	int nparams;
	const char *param_names;
	const char *image;
	const char *label;
};

static TypeInfo aperture_types[] = {
	{ 'c', 1, "D", "circle", "Circular" },
	{ 'h', 1, "D", "hexagon", "Hexagonal" },
	{ 't', 1, "D", "triangle", "Triangular" },
	{ 'r', 2, "W H", "rectangle", "Rectangular" },
	{ 'l', 3, "X1 X2 L", "troughsection", "Single Axis Curvature Section" },
	{ 'a', 3, "R1 R2 Theta", "annulus", "Annular" },
	{ 'i', 6, "X1 Y1 X2 Y2 X3 Y3", "irreg_tri", "Irregular Triangle" },
	{ 'q', 8, "X1 Y1 X2 Y2 X3 Y3 X4 Y4", "irreg_quad", "Irregular Quadrilateral" },
	{ 0, 0, NULL, NULL } };

static TypeInfo surface_types[] = {
	{ 'm', -1, "", "zernike", "Zernike Series File" },
	{ 'v', -1, "", "vshot", "VSHOT Data File" },
	{ 'r', -1, "", "polynomial", "Rotationally Symmetric Polynomial File" },
	{ 'i', -1, "", "cubicspline", "Rotationally Symmetric Cubic Spline File" },
	{ 'e', -1, "", "finiteelement", "Finite Element Data File" }, 
	{ 'p', 2, "Cx Cy", "parabola", "Parabolic" },
	{ 's', 1, "C", "sphere", "Spherical" },
	{ 'o', 2, "C Kappa", "hyperhemi", "Hyperboloids and Hemiellipsoids" },
/*	{ 'g', 8, "Cx Cy Kappa Alpha1 Alpha2 Alpha3 Alpha4 Alpha5", "spencermurty", "General Spencer Murty Description" }, */
	{ 'f', 0, "", "flat", "Flat" },
	{ 'c', 1, "Theta", "cone", "Conical" },
	{ 't', 1, "1/R", "cylinder", "Cylindrical" },
/*	{ 'd', 2, "Ra Rs", "torus", "Torus" }, */
	{ 0, 0, NULL, NULL } };

enum { ID_TYPE = wxID_HIGHEST+491, ID_SELECT_FILE };

class ApertureSurfaceDialog : public wxDialog
{
	const static size_t NMAXPARAMS = 8;
	TypeInfo *m_typeInfo;
	wxChoice *m_type;
	wxStaticBitmap *m_bitmap;
	wxStaticText *m_labels[ NMAXPARAMS ];
	wxNumericCtrl *m_params[ NMAXPARAMS ];
	wxTextCtrl *m_fileName;
	wxButton *m_btnFileSel;
public:
	enum { APERTURE, SURFACE };

	ApertureSurfaceDialog( wxWindow *parent, int mode, const wxString &title )
		: wxDialog( parent, wxID_ANY, title,
			wxDefaultPosition, wxSize(600,470), wxDEFAULT_DIALOG_STYLE | wxRESIZE_BORDER ),
		m_typeInfo( (mode==APERTURE)?aperture_types:surface_types )
	{
		m_type = new wxChoice( this, ID_TYPE );
		
		size_t i=0;
		while( m_typeInfo[i].ap_t != 0 )
			m_type->Append( m_typeInfo[i++].label );

		m_type->SetSelection( 0 );

		m_bitmap = new wxStaticBitmap( this, wxID_ANY, wxNullBitmap, wxDefaultPosition, wxSize(370,250) );

		wxFlexGridSizer *grid = new wxFlexGridSizer( 2, wxSize(2,2) );

		for( i=0;i<NMAXPARAMS;i++ )
		{
			m_labels[i] = new wxStaticText( this, wxID_ANY, wxEmptyString );
			m_params[i] = new wxNumericCtrl( this );

			grid->Add( m_labels[i], 0, wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT|wxLEFT|wxRIGHT, 3 );
			grid->Add( m_params[i], 0, wxALL|wxALIGN_CENTER_VERTICAL, 0 );
		}

		wxBoxSizer *bitgridbox = new wxBoxSizer( wxHORIZONTAL );
		bitgridbox->Add( m_bitmap, 1, wxALL|wxEXPAND, 10 );
		bitgridbox->Add( grid, 0, wxALL, 10 );

		m_fileName = new wxTextCtrl( this, wxID_ANY );
		m_btnFileSel = new wxButton( this, ID_SELECT_FILE, "...", wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT );

		wxBoxSizer *filebox = new wxBoxSizer( wxHORIZONTAL );
		filebox->Add( m_fileName, 1, wxALL|wxALIGN_CENTER_VERTICAL, 10 );
		filebox->Add( m_btnFileSel, 0, wxRIGHT|wxALIGN_CENTER_VERTICAL, 10 );


		wxBoxSizer *sizer = new wxBoxSizer( wxVERTICAL );
		sizer->Add( m_type, 0, wxALL, 10 );
		sizer->Add( bitgridbox, 1, wxALL|wxEXPAND, 0 );
		sizer->Add( filebox, 0, wxALL|wxEXPAND, 0 );
		sizer->Add( CreateButtonSizer( wxOK|wxCANCEL ), 0, wxALL|wxEXPAND, 15 );

		SetSizer( sizer );

		UpdateImageAndLabels();
	}

	void UpdateImageAndLabels()
	{
		int sel = m_type->GetSelection();
		if ( sel < 0 ) return;

		TypeInfo &ty = m_typeInfo[sel];
		wxArrayString parts( wxStringTokenize( ty.param_names ) );

		for( size_t i=0;i<NMAXPARAMS;i++ )
		{
			if ( i < parts.size() )
			{
				m_labels[i]->SetLabel( parts[i] );
				m_labels[i]->Show( true );
				m_params[i]->Show( true );
			}
			else
			{
				m_labels[i]->Show( false );
				m_params[i]->Show( false );
			}
		}

		m_fileName->Show( ty.nparams < 0 );
		m_btnFileSel->Show( ty.nparams < 0 );

		m_bitmap->SetBitmap( wxBitmap( MainWindow::Instance().GetAppDataDir() + "/images/" + wxString(ty.image) + ".png", wxBITMAP_TYPE_PNG ) );
		
		// avoid some funky display issues when changing bitmaps
		Layout();
		m_bitmap->Refresh();
		Refresh();
	}

	void OnCommand( wxCommandEvent &evt )
	{
		switch( evt.GetId() )
		{
		case ID_TYPE:
			UpdateImageAndLabels();
			break;
		case ID_SELECT_FILE:
			{
				int sel = m_type->GetSelection();
				if ( sel < 0 ) return;
				TypeInfo &ty = m_typeInfo[sel];
				wxFileDialog dlg( this, "Select " + wxString(ty.label), wxPathOnly( m_fileName->GetValue() ), m_fileName->GetValue(), "All Files (*.*)|*.*", wxFD_OPEN );
				if ( wxID_OK == dlg.ShowModal() )
				{
					wxString file( dlg.GetPath() );
					if ( wxPathOnly(MainWindow::Instance().GetFileName()) == wxPathOnly(file)
							|| MainWindow::Instance().GetWorkDir() == wxPathOnly(file) )
						m_fileName->ChangeValue( wxFileNameFromPath( file ) );
					else
						m_fileName->ChangeValue( file );
				}
			}
			break;
		}
	}

	void Set( const wxString &apstr )
	{
		if ( apstr.Len() < 3 ) return;
		char ap_t = tolower( (char)apstr[0] );
		
		int it=-1;

		size_t k=0;
		while( m_typeInfo[k].ap_t != 0 )
		{
			if ( m_typeInfo[k].ap_t == ap_t )
			{
				it = (int)k;
				break;
			}
			k++;
		}

		if ( it < 0 ) return; // invalid aperture or surface

		m_type->SetSelection(it);

		if (m_typeInfo[it].nparams < 0)
		{
			m_fileName->ChangeValue( apstr.Mid(2) );
		}
		else
		{
			wxArrayString parts( wxStringTokenize( apstr.Mid(2), "," ) );
			for (size_t i=0;i<NMAXPARAMS;i++)
			{
				if (i < parts.size())
					m_params[i]->SetValue( wxAtof(parts[i]) );
				else
					m_params[i]->SetValue( 0.0 );
			}
		}

		UpdateImageAndLabels();
	}

	wxString Get() const
	{
		int sel = m_type->GetSelection();
		if ( sel < 0 ) return wxEmptyString;

		wxString apstr( m_typeInfo[ sel ].ap_t );

		if ( m_typeInfo[sel].nparams < 0)
		{
			apstr += ":" + m_fileName->GetValue();
		}
		else
		{
			apstr += "-";
			for (int i=0;i<NMAXPARAMS;i++)
			{
				if (i<m_typeInfo[sel].nparams)
					apstr += wxString::Format("%lg", m_params[i]->Value());
				else
					apstr += "0";

				if (i < NMAXPARAMS-1) apstr += ",";
			}
		}

		return apstr;
	}

	DECLARE_EVENT_TABLE();
};

BEGIN_EVENT_TABLE( ApertureSurfaceDialog, wxDialog )
	EVT_CHOICE( ID_TYPE, ApertureSurfaceDialog::OnCommand )
	EVT_BUTTON( ID_SELECT_FILE, ApertureSurfaceDialog::OnCommand )
END_EVENT_TABLE()



class ZRotationDialog : public wxDialog
{
	double m_x, m_y, m_z, m_ax, m_ay, m_az;
	wxNumericCtrl *m_zrot, *m_step;
	wxNumericCtrl *m_alpha, *m_beta, *m_gamma;
	wxNumericCtrl *m_Xi, *m_Xj, *m_Xk;
	wxNumericCtrl *m_Yi, *m_Yj, *m_Yk;
	wxNumericCtrl *m_xAngle1, *m_xAngle2, *m_xAngle3;
	wxNumericCtrl *m_yAngle1, *m_yAngle2, *m_yAngle3;
	wxStaticText *m_object, *m_unitVectorInfo, *m_error;
	wxStaticText *m_xAngle1Label, *m_xAngle2Label, *m_xAngle3Label;
	wxStaticText *m_yAngle1Label, *m_yAngle2Label, *m_yAngle3Label;
public:
	enum { _1=wxID_HIGHEST+912, ID_ZROT_ANGLE, ID_UP, ID_DOWN };

	ZRotationDialog( wxWindow *parent, const wxString &title )
		: wxDialog( parent, wxID_ANY, title, wxDefaultPosition, wxSize( 700, 500 ),
			wxDEFAULT_DIALOG_STYLE | wxRESIZE_BORDER )
	{
		m_x = m_y = m_z = m_ax = m_ay = m_az = 0.0;

		wxBoxSizer *sizer1 = new wxBoxSizer( wxHORIZONTAL );
		sizer1->Add( new wxStaticText( this, wxID_ANY, "Z rotation:"), 0, wxLEFT|wxRIGHT|wxALIGN_CENTER_VERTICAL, 3 );
		sizer1->Add( m_zrot = new wxNumericCtrl( this, ID_ZROT_ANGLE ), 0, wxALIGN_CENTER_VERTICAL, 0 );
		sizer1->Add( new wxButton( this, ID_UP, "Up" ), 0, wxALIGN_CENTER_VERTICAL|wxALL, 3 );
		sizer1->Add( new wxButton( this, ID_DOWN, "Down" ), 0, wxALIGN_CENTER_VERTICAL|wxALL, 3 );
		sizer1->AddStretchSpacer();
		sizer1->Add( new wxStaticText( this, wxID_ANY, "Up/down step size:"), 0, wxLEFT|wxRIGHT|wxALIGN_CENTER_VERTICAL, 3 );
		sizer1->Add( m_step = new wxNumericCtrl( this, wxID_ANY, 0.005 ), 0, wxALL|wxALIGN_CENTER_VERTICAL, 0 );

		wxBoxSizer *sizer2 = new wxBoxSizer( wxVERTICAL );
		sizer2->Add( m_object = new wxStaticText( this, wxID_ANY, "obj" ), 0, wxALIGN_CENTER_VERTICAL|wxALL, 10 );
		sizer2->Add( sizer1, 0, wxALL|wxEXPAND, 10 );

		wxFont font( m_object->GetFont() );
		font.SetWeight( wxFONTWEIGHT_BOLD );
		m_object->SetFont( font );

		wxFlexGridSizer *grid1 = new wxFlexGridSizer( 4, wxSize(2, 2) );
			
		grid1->AddStretchSpacer();
		grid1->Add( new wxStaticText( this, wxID_ANY, "Alpha"), 0, wxALIGN_CENTER_VERTICAL|wxLEFT|wxRIGHT, 5 );
		grid1->Add( new wxStaticText( this, wxID_ANY, "Beta" ), 0, wxALIGN_CENTER_VERTICAL|wxLEFT|wxRIGHT, 5 );
		grid1->Add( new wxStaticText( this, wxID_ANY, "Gamma"), 0, wxALIGN_CENTER_VERTICAL|wxLEFT|wxRIGHT, 5 );
			
		grid1->Add( new wxStaticText( this, wxID_ANY, "Euler angles"), 0, wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT|wxLEFT|wxRIGHT, 5 );
		grid1->Add( m_alpha = NewIndicator(), 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );
		grid1->Add( m_beta  = NewIndicator(), 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );
		grid1->Add( m_gamma = NewIndicator(), 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );

		grid1->AddStretchSpacer();
		grid1->Add( new wxStaticText( this, wxID_ANY, "i"), 0, wxALIGN_CENTER_VERTICAL|wxLEFT|wxRIGHT, 5 );
		grid1->Add( new wxStaticText( this, wxID_ANY, "j" ), 0, wxALIGN_CENTER_VERTICAL|wxLEFT|wxRIGHT, 5 );
		grid1->Add( new wxStaticText( this, wxID_ANY, "k"), 0, wxALIGN_CENTER_VERTICAL|wxLEFT|wxRIGHT, 5 );
			
		grid1->Add( new wxStaticText( this, wxID_ANY, "X-axis unit vector"), 0, wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT|wxLEFT|wxRIGHT, 5 );
		grid1->Add( m_Xi = NewIndicator(), 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );
		grid1->Add( m_Xj = NewIndicator(), 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );
		grid1->Add( m_Xk = NewIndicator(), 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );

		grid1->Add( new wxStaticText( this, wxID_ANY, "Y-axis unit vector"), 0, wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT|wxLEFT|wxRIGHT, 5 );
		grid1->Add( m_Yi = NewIndicator(), 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );
		grid1->Add( m_Yj = NewIndicator(), 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );
		grid1->Add( m_Yk = NewIndicator(), 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );
	
		grid1->AddStretchSpacer();
		grid1->AddStretchSpacer();
		grid1->AddStretchSpacer();
		grid1->Add( m_unitVectorInfo = new wxStaticText( this, wxID_ANY, "in global coordinate system"), 0, wxALIGN_CENTER_VERTICAL|wxLEFT|wxRIGHT, 5 );
		
	
		wxFlexGridSizer *grid2 = new wxFlexGridSizer( 3, wxSize(2,2) );
				
		grid2->Add( new wxStaticText( this, wxID_ANY, "X-axis makes angle of"), 0, wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT|wxLEFT|wxRIGHT, 5 );
		grid2->Add( m_xAngle1 = NewIndicator(), 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );
		grid2->Add( m_xAngle1Label = new wxStaticText( this, wxID_ANY, wxEmptyString ), 0, wxALIGN_CENTER_VERTICAL|wxLEFT|wxRIGHT, 5 );
		
		grid2->AddStretchSpacer();
		grid2->Add( m_xAngle2 = NewIndicator(), 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );
		grid2->Add( m_xAngle2Label = new wxStaticText( this, wxID_ANY, wxEmptyString ), 0, wxALIGN_CENTER_VERTICAL|wxLEFT|wxRIGHT, 5 );
		
		grid2->AddStretchSpacer();
		grid2->Add( m_xAngle3 = NewIndicator(), 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );
		grid2->Add( m_xAngle3Label = new wxStaticText( this, wxID_ANY, wxEmptyString ), 0, wxALIGN_CENTER_VERTICAL|wxLEFT|wxRIGHT, 5 );
					
		grid2->Add( new wxStaticText( this, wxID_ANY, "Y-axis makes angle of"), 0, wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT|wxLEFT|wxRIGHT, 5 );
		grid2->Add( m_yAngle1 = NewIndicator(), 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );
		grid2->Add( m_yAngle1Label = new wxStaticText( this, wxID_ANY, wxEmptyString ), 0, wxALIGN_CENTER_VERTICAL|wxLEFT|wxRIGHT, 5 );
		
		grid2->AddStretchSpacer();
		grid2->Add( m_yAngle2 = NewIndicator(), 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );
		grid2->Add( m_yAngle2Label = new wxStaticText( this, wxID_ANY, wxEmptyString ), 0, wxALIGN_CENTER_VERTICAL|wxLEFT|wxRIGHT, 5 );
		
		grid2->AddStretchSpacer();
		grid2->Add( m_yAngle3 = NewIndicator(), 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );
		grid2->Add( m_yAngle3Label = new wxStaticText( this, wxID_ANY, wxEmptyString ), 0, wxALIGN_CENTER_VERTICAL|wxLEFT|wxRIGHT, 5 );
		
		sizer2->Add( new wxStaticLine( this ), 0, wxALL|wxEXPAND, 3 );
		sizer2->Add( grid1, 1, wxALL, 10 );
		sizer2->Add( grid2, 1, wxALL, 10 );
		sizer2->Add( m_error = new wxStaticText( this, wxID_ANY, wxEmptyString ), 0, wxALL|wxEXPAND, 10 );
		m_error->SetForegroundColour( *wxRED );
		
		sizer2->Add( new wxStaticLine( this ), 0, wxALL|wxEXPAND, 3 );
		sizer2->Add( CreateButtonSizer( wxOK|wxCANCEL ), 0, wxALL|wxEXPAND, 10 );

		SetSizer( sizer2 );
	}

	wxNumericCtrl *NewIndicator()
	{
		wxNumericCtrl *num = new wxNumericCtrl( this );
		num->SetEditable( false );
		num->SetForegroundColour( *wxBLUE );
		return num;
	}


	enum CoordSys { STAGE, ELEMENT };

	void Setup( CoordSys mode, double zrot, double x, double y, double z, 
		double ax, double ay, double az, const wxString &obj )
	{
		m_x = x;
		m_y = y;
		m_z = z;
		m_ax = ax;
		m_ay = ay;
		m_az = az;

		m_zrot->SetValue(zrot);

		m_object->SetLabel( obj + wxString::Format(":  origin (%lg, %lg, %lg) aim point (%lg, %lg, %lg)", x, y, z, ax, ay, az ) );
		
		if (mode == STAGE)
		{
			m_unitVectorInfo->SetLabel("in global coordinate system");
			m_xAngle1Label->SetLabel("with XY plane of global coordinate system");
			m_xAngle2Label->SetLabel("with YZ plane of global coordinate system");
			m_xAngle3Label->SetLabel("with XZ plane of global coordinate system");
			m_yAngle1Label->SetLabel("with XY plane of global coordinate system");
			m_yAngle2Label->SetLabel("with YZ plane of global coordinate system");
			m_yAngle3Label->SetLabel("with XZ plane of global coordinate system");
		}
		else
		{
			m_unitVectorInfo->SetLabel("in stage coordinate system");
			m_xAngle1Label->SetLabel("with XY plane of stage coordinate system");
			m_xAngle2Label->SetLabel("with YZ plane of stage coordinate system");
			m_xAngle3Label->SetLabel("with XZ plane of stage coordinate system");
			m_yAngle1Label->SetLabel("with XY plane of stage coordinate system");
			m_yAngle2Label->SetLabel("with YZ plane of stage coordinate system");
			m_yAngle3Label->SetLabel("with XZ plane of stage coordinate system");
		}

		Fit();
		Refresh();

		Calculate();
	}

	void Calculate()
	{
		double zrot = m_zrot->Value();

		double CosRefZ[3],Euler[3],CosLoc[3],CosRefX[3],CosRefY[3],Alpha,Beta,Gamma,
			CosAlpha,CosBeta,CosGamma,
			SinAlpha,SinBeta,SinGamma,
			RRefToLoc[3][3],RLocToRef[3][3];

		double dx = m_ax - m_x;
		double dy = m_ay - m_y;
		double dz = m_az - m_z;
		double dtot = sqrt(dx*dx + dy*dy + dz*dz);

		if (dtot == 0.0)
		{
			m_error->SetLabel( "Stage geometry as defined not possible. Check inputs." );
			Layout();
			return;
		}
		else
		{
			bool noerr = m_error->GetLabel().IsEmpty();
			m_error->SetLabel( wxEmptyString );
			if ( !noerr ) Layout();
		}

		dx /= dtot;
		dy /= dtot;
		dz /= dtot;

		CosRefZ[0] = dx;
		CosRefZ[1] = dy;
		CosRefZ[2] = dz;

		Euler[0]  = atan2(dx,dz);
		Euler[1] = asin(dy);
		Euler[2] = zrot*(acos(-1.0)/180.0);

		Alpha = Euler[0];
		Beta  = Euler[1];
		Gamma = Euler[2];
		CosAlpha = cos(Alpha);
		CosBeta  = cos(Beta);
		CosGamma = cos(Gamma);
		SinAlpha = sin(Alpha);
		SinBeta  = sin(Beta);
		SinGamma = sin(Gamma);

		//{Fill in elements of the transformation matrix as per Spencer and Murty paper
		// page 673 equation (2)}
		RRefToLoc[0][0] = CosAlpha*CosGamma + SinAlpha*SinBeta*SinGamma;
		RRefToLoc[0][1] = -CosBeta*SinGamma;
		RRefToLoc[0][2] = -SinAlpha*CosGamma + CosAlpha*SinBeta*SinGamma;
		RRefToLoc[1][0] = CosAlpha*SinGamma - SinAlpha*SinBeta*CosGamma;
		RRefToLoc[1][1] = CosBeta*CosGamma;
		RRefToLoc[1][2] = -SinAlpha*SinGamma - CosAlpha*SinBeta*CosGamma;
		RRefToLoc[2][0] = SinAlpha*CosBeta;
		RRefToLoc[2][1] = SinBeta;
		RRefToLoc[2][2] = CosAlpha*CosBeta;

		m_alpha->SetValue(Alpha*(180.0/acos(-1.0)));
		m_beta->SetValue(Beta*(180.0/acos(-1.0)));
		m_gamma->SetValue(Gamma*(180.0/acos(-1.0)));

		//{Transpose the matrix to get the inverse transformation matrix for going back
		// to reference system.  This is used by the TransformToReference procedure.}
		::st_matrix_transpose( RRefToLoc, RLocToRef );

		//first x-axis
		CosLoc[0] = 1;
		CosLoc[1] = 0;
		CosLoc[2] = 0;

		::st_matrix_vector_mult(RLocToRef, CosLoc, CosRefX);

		m_Xi->SetValue(CosRefX[0]);
		m_Xj->SetValue(CosRefX[1]);
		m_Xk->SetValue(CosRefX[2]);

		//second y-axis
		CosLoc[0] = 0;
		CosLoc[1] = 1;
		CosLoc[2] = 0;
		::st_matrix_vector_mult(RLocToRef, CosLoc, CosRefY);

		m_Yi->SetValue(CosRefY[0]);
		m_Yj->SetValue(CosRefY[1]);
		m_Yk->SetValue(CosRefY[2]);

		m_xAngle1->SetValue(asin(CosRefX[2])*(180.0/acos(-1.0)));
		m_xAngle2->SetValue(asin(CosRefX[0])*(180.0/acos(-1.0)));
		m_xAngle3->SetValue(asin(CosRefX[1])*(180.0/acos(-1.0)));

		m_yAngle1->SetValue(asin(CosRefY[2])*(180.0/acos(-1.0)));
		m_yAngle2->SetValue(asin(CosRefY[0])*(180.0/acos(-1.0)));
		m_yAngle3->SetValue(asin(CosRefY[1])*(180.0/acos(-1.0)));
	}

	double GetZRot()
	{
		return m_zrot->Value();
	}

	void OnCommand( wxCommandEvent &evt )
	{
		switch( evt.GetId() )
		{
		case ID_UP:
			m_zrot->SetValue( m_zrot->Value() + m_step->Value() );
			break;
		case ID_DOWN:
			m_zrot->SetValue( m_zrot->Value() - m_step->Value() );
			break;
		}
			Calculate();
	}

	DECLARE_EVENT_TABLE();
};

BEGIN_EVENT_TABLE( ZRotationDialog, wxDialog )
	EVT_NUMERIC( ZRotationDialog::ID_ZROT_ANGLE, ZRotationDialog::OnCommand )
	EVT_BUTTON( ZRotationDialog::ID_UP, ZRotationDialog::OnCommand )
	EVT_BUTTON( ZRotationDialog::ID_DOWN, ZRotationDialog::OnCommand )
END_EVENT_TABLE()

enum Columns
{
    Col_Enable,
    Col_X,
    Col_Y,
    Col_Z,
    Col_AX,
    Col_AY,
    Col_AZ,
	Col_ZRot,
	Col_Aperture,
	Col_Surface,
	Col_Interaction,
	Col_Optics,
	Col_Comment,
    Col_Max
};

static const char *headers[Col_Max] =
{
    "En.",
    "X-Coord.",
    "Y-Coord.",
    "Z-Coord.",
    "X-AimPt.",
    "Y-AimPt.",
    "Z-AimPt.",
    "Z-Rot.",
	"Aperture",
	"Surface",
	"Interaction",
	"Optics",
	"Comment"	
};

static int colwidths[Col_Max] =
{
	24,
	70,
	70,
	70,
	70,
	70,
	70,
	70,
	170,
	170,
	95,
	90,
	110
};

class ElementTable : public wxGridTableBase
{
	Stage *m_stage;
	GeometryForm *m_geoForm;

	int m_startRow;
	int m_endRow;

public:
    ElementTable(GeometryForm *gf)
	{
		m_geoForm = gf;
		m_stage = 0;
	}

	void SetStage(Stage *s)
	{
		m_stage = s;
		m_startRow = m_endRow = -1;
	}

	void SetEditRangeRows(int s, int e)
	{
		m_startRow = s;
		m_endRow = e;
	}

	void GetEditRange(int *s, int *e)
	{
		if (s) *s = m_startRow;
		if (e) *e = m_endRow;
	}

    virtual int GetNumberRows()
	{
		if (!m_stage) return 0;
		else return m_stage->ElementList.size();
	}

	virtual int GetNumberCols() { return Col_Max; }

    virtual bool IsEmptyCell( int row, int col )
	{
		if (!m_stage || row >= m_stage->ElementList.size()
			|| col >= Col_Max) return true;

		return false;
	}

    virtual wxString GetValue( int row, int col )
	{
		Element *elm = GetElement(row);
		if (!elm) return wxEmptyString;

		wxString val;
		switch(col)
		{
		case Col_Enable:  val.Printf("%d", (int) elm->Enabled ); break;
		case Col_X:  val.Printf("%lg", elm->X); break;
		case Col_Y:  val.Printf("%lg", elm->Y); break;
		case Col_Z:  val.Printf("%lg", elm->Z); break;
		case Col_AX:  val.Printf("%lg", elm->AX); break;
		case Col_AY:  val.Printf("%lg", elm->AY); break;
		case Col_AZ:  val.Printf("%lg", elm->AZ); break;
		case Col_ZRot: val.Printf("%lg", elm->ZRot); break;

		case Col_Aperture:
			val.Printf("%c-%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg",
				elm->ApertureIndex,
				elm->ApertureParams[0],
				elm->ApertureParams[1],
				elm->ApertureParams[2],
				elm->ApertureParams[3],
				elm->ApertureParams[4],
				elm->ApertureParams[5],
				elm->ApertureParams[6],
				elm->ApertureParams[7] );
			break;
		case Col_Surface:
			if (elm->SurfaceFile.IsEmpty())
			{
				val.Printf("%c-%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg",
					elm->SurfaceIndex,
					elm->SurfaceParams[0],
					elm->SurfaceParams[1],
					elm->SurfaceParams[2],
					elm->SurfaceParams[3],
					elm->SurfaceParams[4],
					elm->SurfaceParams[5],
					elm->SurfaceParams[6],
					elm->SurfaceParams[7] );
			}
			else
			{
				val = wxString(elm->SurfaceIndex) + ":" + elm->SurfaceFile;
			}
			break;
		case Col_Interaction:
			if (elm->InteractionType == Element::REFLECTION)
				val = "Reflection";
			else
				val = "Refraction";
			break;
		case Col_Optics:
			val = elm->OpticName;
			break;
		case Col_Comment:
			val = elm->Comment;
			break;
		default:
			val = "<ERR>";
		}

		return val;
	}

    virtual void SetValue( int row, int col, const wxString& value )
	{
		int row0 = row, row1 = row;
		if (m_startRow >= 0 && m_endRow >= 0 && m_endRow > m_startRow)
		{
			row0 = m_startRow;
			row1 = m_endRow;
		}

		for (row=row0;row<=row1;row++)
		{
			Element *elm = GetElement(row);
			if (!elm) continue;

			switch(col)
			{
			case Col_Enable:  elm->Enabled = (value!="0"); break;
			case Col_X:  elm->X = atof(value.c_str()); break;
			case Col_Y:  elm->Y = atof(value.c_str()); break;
			case Col_Z:  elm->Z = atof(value.c_str()); break;
			case Col_AX:  elm->AX = atof(value.c_str()); break;
			case Col_AY:  elm->AY = atof(value.c_str()); break;
			case Col_AZ:  elm->AZ = atof(value.c_str()); break;
			case Col_ZRot: elm->ZRot = atof(value.c_str()); break;
			case Col_Aperture:
				{
					if (value.Len() < 2) return;
					elm->ApertureIndex = tolower(value[0]);
					if (strchr( "chtrlaiq", elm->ApertureIndex ) == NULL)
						elm->ApertureIndex = 'c';

					wxArrayString parts = wxStringTokenize(value.Mid(2), ",");
					for (int i=0;i<8;i++)
					{
						if (i < parts.Count())
							elm->ApertureParams[i] = atof( parts[i].c_str() );
						else
							elm->ApertureParams[i] = 0.0;
					}
				}
				break;
			case Col_Surface:
				{
					if (value.Len() < 2) return;
					elm->SurfaceIndex = value[0];
					if (value[1] == ':' && strchr("mvrie", elm->SurfaceIndex) != NULL)
					{
						elm->SurfaceFile = value.Mid(2);
					}
					else
					{
						if (strchr("psogfctd", elm->SurfaceIndex) == NULL)
							elm->SurfaceIndex = 't';

						wxArrayString parts = wxStringTokenize(value.Mid(2), ",");
						elm->SurfaceFile.Empty();
						for (int j=0;j<8;j++)
						{
							if (j<parts.Count())
								elm->SurfaceParams[j] = atof( parts[j].c_str() );
							else
								elm->SurfaceParams[j] = 0.0;
						}
					}
				}
				break;
			case Col_Interaction:
				if (value.Lower().Left(5) == "refra")
					elm->InteractionType = Element::REFRACTION;
				else
					elm->InteractionType = Element::REFLECTION;
				break;
			case Col_Optics:
				elm->OpticName = value;
				break;
			case Col_Comment:
				elm->Comment = value;
				break;
			}
		}

		m_geoForm->Modified();
	}

    virtual wxString GetColLabelValue( int col )
	{
		return headers[col];
	}

    virtual wxString GetTypeName( int row, int col )
	{
		switch ( col )
		{
		case Col_Enable:
			return wxGRID_VALUE_BOOL;
		case Col_X:
		case Col_Y:
		case Col_Z:
		case Col_AX:
		case Col_AY:
		case Col_AZ:
		case Col_ZRot:
		case Col_Aperture:
		case Col_Surface:
			return wxGRID_VALUE_STRING;
		case Col_Interaction:
			return wxString::Format("%s:Refraction,Reflection", wxGRID_VALUE_CHOICE);
		case Col_Optics:
		case Col_Comment:
			return wxGRID_VALUE_STRING;
		default:
			return wxEmptyString;
		}
	}

    virtual bool CanGetValueAs( int row, int col, const wxString& typeName )
	{
		if (typeName == wxGRID_VALUE_BOOL) 
			return (col == Col_Enable);
	
		if (typeName == wxGRID_VALUE_STRING)
			return (col==Col_X || col==Col_Y || col==Col_Z 
				|| col==Col_AX || col==Col_AY || col==Col_AZ || col==Col_ZRot
				|| col==Col_Aperture||col==Col_Surface
				|| col==Col_Interaction || col==Col_Optics ||col==Col_Comment);
	
		return false;
	}

    virtual bool CanSetValueAs( int row, int col, const wxString& typeName )
	{
		return CanGetValueAs(row,col,typeName);
	}

	virtual bool GetValueAsBool( int row, int col )
	{
		Element *elm = GetElement(row);
		if (!elm) return false;
		if (col == Col_Enable)
			return elm->Enabled;
		else
			return false;
	}

	virtual void SetValueAsBool( int row, int col, bool value )
	{
		int row0 = row, row1 = row;
		if (m_startRow >= 0 && m_endRow >= 0 && m_endRow > m_startRow)
		{
			row0 = m_startRow;
			row1 = m_endRow;
		}

		for (row=row0;row<=row1;row++)
		{
			Element *elm = GetElement(row);
			if (!elm) continue;
			if (col == Col_Enable)
				elm->Enabled = value;
		}
	}
	
	virtual bool AppendRows(size_t nrows)
	{	
		if (!m_stage||nrows<1) return true;
		for (int i=0;i<(int)nrows;i++)
			m_stage->ElementList.push_back( new Element );

		m_geoForm->Modified();
	
		if ( GetView() )
		{
			wxGridTableMessage msg( this,
									wxGRIDTABLE_NOTIFY_ROWS_APPENDED,
									nrows );

			GetView()->ProcessTableMessage( msg );
		}
		return true;
	}

	virtual bool InsertRows(size_t pos, size_t nrows)
	{
		if (!m_stage||nrows<1) return true;
		m_geoForm->Modified();

		if (pos < 0) pos = 0;
		if (pos > m_stage->ElementList.size()) pos = m_stage->ElementList.size();
		for (int i=0;i<(int)nrows;i++)
		{
			m_stage->ElementList.insert( m_stage->ElementList.begin() + pos, new Element );
		}
	
		if ( GetView() )
		{
			wxGridTableMessage msg( this,
									wxGRIDTABLE_NOTIFY_ROWS_INSERTED,
									pos,
									nrows );

			GetView()->ProcessTableMessage( msg );
		}

		return true;
	}

	virtual bool DeleteRows(size_t pos, size_t nrows)
	{
		if (!m_stage||nrows<1) return true;
		
		m_geoForm->Modified();

		if (nrows > m_stage->ElementList.size() - pos)
			nrows = m_stage->ElementList.size() - pos;

		int nremaining = nrows;
		while (nremaining>0)
		{
			if (pos < m_stage->ElementList.size())
			{
				delete m_stage->ElementList[pos];
				m_stage->ElementList.erase( m_stage->ElementList.begin() + pos);
			}
			nremaining--;
		}
			
		if ( GetView() )
		{
			wxGridTableMessage msg( this,
									wxGRIDTABLE_NOTIFY_ROWS_DELETED,
									pos,
									nrows );

			GetView()->ProcessTableMessage( msg );
		}

		return true;
	}

private:
	Element *GetElement(int id)
	{
		if (m_stage && id >= 0 && id < m_stage->ElementList.size())
			return m_stage->ElementList[id];
		else
			return NULL;
	}

};


enum { ID_ELEMENT_GRID = wxID_HIGHEST + 938,
	ID_STAGE_NAME, ID_VIRTUAL_STAGE, ID_MULTIPLE_HITS, ID_TRACE_THROUGH,
	ID_X, ID_Y, ID_Z, ID_AX, ID_AY, ID_AZ, ID_ZROT,
	ID_ELEMENT_INSERT, ID_ELEMENT_ZROT, ID_ELEMENT_APPEND, ID_ELEMENT_APERTURE,
	ID_ELEMENT_DELETE, ID_ELEMENT_SURFACE, ID_ELEMENT_DELETE_ALL, ID_ELEMENT_OPTICS,
	ID_STAGE_ZROT,
	ID_IMPORT, ID_EXPORT };

BEGIN_EVENT_TABLE( StageForm, wxPanel )

	EVT_GRID_CMD_CELL_CHANGED( ID_ELEMENT_GRID, StageForm::OnGridCellChange)
	EVT_GRID_CMD_SELECT_CELL( ID_ELEMENT_GRID, StageForm::OnGridCellSelect)
	EVT_GRID_CMD_RANGE_SELECTED( ID_ELEMENT_GRID, StageForm::OnGridRangeSelect)  //<< use for wx3.5+
	//EVT_GRID_CMD_RANGE_SELECT( ID_ELEMENT_GRID, StageForm::OnGridRangeSelect)  //<< use for wx3.1-3.4
	EVT_GRID_CMD_CELL_LEFT_DCLICK( ID_ELEMENT_GRID, StageForm::OnGridCellDClick )
	EVT_GRID_CMD_CELL_RIGHT_CLICK( ID_ELEMENT_GRID, StageForm::OnGridCellRightClick )

	EVT_TEXT_ENTER( ID_STAGE_NAME, StageForm::OnCommand )
	EVT_CHECKBOX( ID_VIRTUAL_STAGE, StageForm::OnCommand )
	EVT_CHECKBOX( ID_MULTIPLE_HITS, StageForm::OnCommand )
	EVT_CHECKBOX( ID_TRACE_THROUGH, StageForm::OnCommand )
	EVT_NUMERIC( ID_X, StageForm::OnCommand )
	EVT_NUMERIC( ID_Y, StageForm::OnCommand )
	EVT_NUMERIC( ID_Z, StageForm::OnCommand )
	EVT_NUMERIC( ID_AX, StageForm::OnCommand )
	EVT_NUMERIC( ID_AY, StageForm::OnCommand )
	EVT_NUMERIC( ID_AZ, StageForm::OnCommand )
	EVT_NUMERIC( ID_ZROT, StageForm::OnCommand )

	EVT_BUTTON( ID_ELEMENT_INSERT, StageForm::OnCommand )
	EVT_BUTTON( ID_ELEMENT_ZROT, StageForm::OnCommand )
	EVT_BUTTON( ID_ELEMENT_APPEND, StageForm::OnCommand )
	EVT_BUTTON( ID_ELEMENT_APERTURE,StageForm::OnCommand )
	EVT_BUTTON( ID_ELEMENT_DELETE, StageForm::OnCommand )
	EVT_BUTTON( ID_ELEMENT_SURFACE, StageForm::OnCommand )
	EVT_BUTTON( ID_ELEMENT_DELETE_ALL, StageForm::OnCommand )
	EVT_BUTTON( ID_ELEMENT_OPTICS,StageForm::OnCommand )
	EVT_BUTTON( ID_STAGE_ZROT, StageForm::OnCommand )
	EVT_BUTTON( ID_IMPORT, StageForm::OnCommand )
	EVT_BUTTON( ID_EXPORT, StageForm::OnCommand )
	EVT_BUTTON( wxID_COPY, StageForm::OnCommand )
	EVT_BUTTON( wxID_PASTE, StageForm::OnCommand )
END_EVENT_TABLE()

StageForm::StageForm(wxWindow *parent, GeometryForm *geo, Project &prj, Stage *s)
	: wxPanel(parent), m_prj(prj)
{
	m_stage = s;
	m_geoForm = geo;
	
	wxStaticBoxSizer *sizer1 = new wxStaticBoxSizer( wxVERTICAL, this, "Stage properties" );
	m_stageName = new wxExtTextCtrl( sizer1->GetStaticBox(), ID_STAGE_NAME, "New stage" );
	m_virtualStage = new wxCheckBox( sizer1->GetStaticBox(), ID_VIRTUAL_STAGE, "Virtual stage" );
	m_multipleHits = new wxCheckBox( sizer1->GetStaticBox(), ID_MULTIPLE_HITS, "Multiple hits per ray" );
	m_traceThrough = new wxCheckBox( sizer1->GetStaticBox(), ID_TRACE_THROUGH, "Trace through" );

	wxBoxSizer *name_sizer = new wxBoxSizer( wxHORIZONTAL );
	name_sizer->Add( new wxStaticText( sizer1->GetStaticBox(), wxID_ANY, "Name" ), 0, wxLEFT|wxRIGHT|wxALIGN_CENTER_VERTICAL, 5 );
	name_sizer->Add( m_stageName, 0, wxALL|wxEXPAND, 5 );

	sizer1->Add( name_sizer, 0, wxALL|wxEXPAND, 0 );
	sizer1->Add( m_virtualStage, 0, wxALL, 5 );
	sizer1->Add( m_multipleHits, 0, wxALL, 5 );
	sizer1->Add( m_traceThrough, 0, wxALL, 5 );

	wxStaticBoxSizer *sizer2 = new wxStaticBoxSizer( wxVERTICAL, this, "Global coordinates" );

	wxFlexGridSizer *sizer3 = new wxFlexGridSizer( 4 );
	
	sizer3->AddStretchSpacer();
	sizer3->Add( new wxStaticText( sizer2->GetStaticBox(), wxID_ANY, "X"), 0, wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER, 2 );
	sizer3->Add( new wxStaticText( sizer2->GetStaticBox(), wxID_ANY, "Y"), 0, wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER, 2 );
	sizer3->Add( new wxStaticText( sizer2->GetStaticBox(), wxID_ANY, "Z"), 0, wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER, 2 );

	sizer3->Add( new wxStaticText( sizer2->GetStaticBox(), wxID_ANY, "Origin"), 0, wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 2 );
	sizer3->Add( m_x = new wxNumericCtrl( sizer2->GetStaticBox(), ID_X ), 0, wxALL, 2 );
	sizer3->Add( m_y = new wxNumericCtrl( sizer2->GetStaticBox(), ID_Y ), 0, wxALL, 2 );
	sizer3->Add( m_z = new wxNumericCtrl( sizer2->GetStaticBox(), ID_Z ), 0, wxALL, 2 );
	
	sizer3->Add( new wxStaticText( sizer2->GetStaticBox(), wxID_ANY, "Aim point"), 0, wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 2 );
	sizer3->Add( m_ax = new wxNumericCtrl( sizer2->GetStaticBox(), ID_AX ), 0, wxALL, 2 );
	sizer3->Add( m_ay = new wxNumericCtrl( sizer2->GetStaticBox(), ID_AY ), 0, wxALL, 2 );
	sizer3->Add( m_az = new wxNumericCtrl( sizer2->GetStaticBox(), ID_AZ, 1.0 ), 0, wxALL, 2 );
	
	sizer3->Add( new wxStaticText( sizer2->GetStaticBox(), wxID_ANY, "Z rotation"), 0, wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 2 );
	sizer3->Add( m_zrot = new wxNumericCtrl( sizer2->GetStaticBox(), ID_ZROT ), 0, wxALL, 2 );
	sizer3->Add( new wxStaticText( sizer2->GetStaticBox(), wxID_ANY, "(deg)"), 0, wxALIGN_CENTER_VERTICAL|wxALIGN_LEFT, 2 );
	sizer3->Add( new wxButton( sizer2->GetStaticBox(), ID_STAGE_ZROT, "Edit Z rotation...", wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT ), 0, wxALIGN_CENTER_VERTICAL, 0 );

	sizer2->Add( sizer3, 1, wxLEFT|wxRIGHT, 10 );

	wxStaticBoxSizer *sizer5 = new wxStaticBoxSizer( wxVERTICAL, this, "Element editing" );
	wxFlexGridSizer *elmbtn_sizer = new wxFlexGridSizer( 2, wxSize(1,1) );
	elmbtn_sizer->Add( new wxButton( sizer5->GetStaticBox(), ID_ELEMENT_INSERT, "Insert..." ) );
	elmbtn_sizer->Add( new wxButton( sizer5->GetStaticBox(), ID_ELEMENT_ZROT, "Z rot..." ) );
	elmbtn_sizer->Add( new wxButton( sizer5->GetStaticBox(), ID_ELEMENT_APPEND, "Append..." ) );
	elmbtn_sizer->Add( new wxButton( sizer5->GetStaticBox(), ID_ELEMENT_APERTURE, "Aperture..." ) );
	elmbtn_sizer->Add( new wxButton( sizer5->GetStaticBox(), ID_ELEMENT_DELETE, "Delete..." ) );
	elmbtn_sizer->Add( new wxButton( sizer5->GetStaticBox(), ID_ELEMENT_SURFACE, "Surface..." ) );
	elmbtn_sizer->Add( new wxButton( sizer5->GetStaticBox(), ID_ELEMENT_DELETE_ALL, "Delete all" ) );
	elmbtn_sizer->Add( new wxButton( sizer5->GetStaticBox(), ID_ELEMENT_OPTICS, "Optics..." ) );
	sizer5->Add( elmbtn_sizer, 1, wxLEFT|wxRIGHT, 10 );

	wxStaticBoxSizer *sizer6 = new wxStaticBoxSizer( wxVERTICAL, this, "Data transfer" );
	wxFlexGridSizer *datbtn_sizer = new wxFlexGridSizer( 1, wxSize(1,1) );
	datbtn_sizer->Add( new wxButton( sizer6->GetStaticBox(), ID_IMPORT, "Import..." ) );
	datbtn_sizer->Add( new wxButton( sizer6->GetStaticBox(), ID_EXPORT, "Export..." ) );
	datbtn_sizer->Add( new wxButton( sizer6->GetStaticBox(), wxID_COPY, "Copy all" ) );
	datbtn_sizer->Add( new wxButton( sizer6->GetStaticBox(), wxID_PASTE, "Paste all" ) );
	sizer6->Add( datbtn_sizer, 1, wxLEFT|wxRIGHT, 10 );

	m_gridTable = NULL;
	m_grid = new wxExtGridCtrl( this, ID_ELEMENT_GRID );
	m_grid->DisableDragCell();
	m_grid->DisableDragRowSize();
	m_grid->DisableDragColMove();
	m_grid->DisableDragGridSize();
	m_grid->SetRowLabelAlignment( wxALIGN_LEFT, wxALIGN_CENTER );

	wxClientDC dc(this);
	dc.SetFont( *wxNORMAL_FONT );
	int height = dc.GetCharHeight() + 4;
	m_grid->SetRowLabelSize( height );
	m_grid->SetColLabelSize( height );
		

	wxBoxSizer *prop_sizer = new wxBoxSizer( wxHORIZONTAL );
	prop_sizer->Add( sizer1, 0, wxALL|wxEXPAND, 5 );
	prop_sizer->Add( sizer2, 0, wxALL|wxEXPAND, 5 );
	prop_sizer->Add( sizer5, 0, wxALL|wxEXPAND, 5 );
	prop_sizer->Add( sizer6, 0, wxALL|wxEXPAND, 5 );
	prop_sizer->AddStretchSpacer();	

	wxBoxSizer *sizer = new wxBoxSizer( wxVERTICAL );
	sizer->Add( prop_sizer, 0, wxALL|wxEXPAND, 0 );
	sizer->Add( m_grid, 1, wxALL|wxEXPAND, 0 );
	SetSizer( sizer );

	UpdateFromData();
}

void StageForm::UpdateProperties()
{
	if ( !m_stage ) return;

	m_stageName->ChangeValue( m_stage->Name );
	m_virtualStage->SetValue( m_stage->Virtual );
	m_multipleHits->SetValue( m_stage->MultiHit );
	m_traceThrough->SetValue( m_stage->TraceThrough );
	m_x->SetValue( m_stage->X );
	m_y->SetValue( m_stage->Y );
	m_z->SetValue( m_stage->Z );
	m_ax->SetValue( m_stage->AX );
	m_ay->SetValue( m_stage->AY );
	m_az->SetValue( m_stage->AZ );
	m_zrot->SetValue( m_stage->ZRot );
}

void StageForm::UpdateFromData()
{
	if (m_gridTable) m_gridTable->SetStage(NULL);

	m_grid->SetTable( NULL );

	m_gridTable = new ElementTable( m_geoForm );
	m_gridTable->SetStage( m_stage );
	m_gridTable->SetAttrProvider( new wxExtGridCellAttrProvider );

	m_grid->SetTable( m_gridTable, true );
	
	for (int i=0;i<Col_Max;i++)
	{
		int W = (int)( colwidths[i] * wxGetScreenHDScale() );
		m_grid->SetColMinimalWidth(i, W);
		m_grid->SetColSize(i, W);
	}

	m_grid->SetRowLabelSize( (int)(60.0 * wxGetScreenHDScale()) );
	m_grid->ForceRefresh();

	UpdateProperties();
}

void StageForm::OnGridCellChange(wxGridEvent &evt)
{
	m_gridTable->SetEditRangeRows(-1,-1);
	evt.Skip();
}

void StageForm::OnGridCellSelect(wxGridEvent &evt)
{
	if (evt.ControlDown())
	{
		evt.Veto();
		return;
	}

	if ( !m_grid->IsCellEditControlShown() )
		m_gridTable->SetEditRangeRows(-1, -1);

	evt.Skip();
}

void StageForm::OnGridRangeSelect(wxGridRangeSelectEvent &evt)
{
	if (evt.Selecting() && !m_grid->IsCellEditControlShown() )
	{
		int top, bottom, left, right;
		m_grid->GetSelRange( &top, &bottom, &left, &right );
		if ( left == right )
			m_gridTable->SetEditRangeRows( top, bottom );
	}
	evt.Skip();
}

void StageForm::OnGridCellRightClick(wxGridEvent &evt)
{
	int start, end;
	m_gridTable->GetEditRange(&start, &end);
	if (start >= 0 && end >= 0)
	{
		int col = m_grid->GetGridCursorCol();
		wxString value = m_grid->GetCellValue(m_grid->GetGridCursorRow(), col);
		for (int r=start;r<=end;r++)
			m_grid->SetCellValue( r, col, value );
	}
	else
		evt.Skip();
}

void StageForm::EditStageZRot()
{
	if (!m_stage) return;


	ZRotationDialog dlg(this, "Edit Stage " + m_stage->Name + " Z Rotation");
	dlg.Setup( ZRotationDialog::STAGE, 
		m_zrot->Value(),
		m_x->Value(),
		m_y->Value(),
		m_z->Value(),
		m_ax->Value(),
		m_ay->Value(),
		m_az->Value(),
		"Stage " + m_stage->Name);

	dlg.CenterOnParent();

	if (dlg.ShowModal() == wxID_OK)
	{
		m_stage->ZRot = dlg.GetZRot();
		m_zrot->SetValue(m_stage->ZRot);
		Modified();
	}
}

void StageForm::EditZRot(int idx)
{
	if (idx < 0) idx = m_grid->GetGridCursorRow();
	if (idx < 0) return;
	
	ZRotationDialog dlg(this, "Edit Element " + wxString::Format("%d",idx+1) + " Z Rotation");
	dlg.Setup( ZRotationDialog::ELEMENT, 
		wxAtof( m_grid->GetCellValue(idx,Col_ZRot) ),
		wxAtof( m_grid->GetCellValue(idx,Col_X) ),
		wxAtof( m_grid->GetCellValue(idx,Col_Y) ),
		wxAtof( m_grid->GetCellValue(idx,Col_Z) ),	
		wxAtof( m_grid->GetCellValue(idx,Col_AX) ),
		wxAtof( m_grid->GetCellValue(idx,Col_AY) ),
		wxAtof( m_grid->GetCellValue(idx,Col_AZ) ),
		wxString::Format("Element %d", idx+1));
	
	dlg.CenterOnParent();
	if (dlg.ShowModal() == wxID_OK)
		m_grid->SetCellValue( idx, Col_ZRot, wxString::Format("%lg",dlg.GetZRot()) );

}

void StageForm::EditAperture(int idx)
{
	if (idx < 0) idx = m_grid->GetGridCursorRow();
	if (idx < 0) return;

	ApertureSurfaceDialog dlg(this, ApertureSurfaceDialog::APERTURE, "Edit element " + wxString::Format("%d",(idx+1)) + " aperture");
	dlg.CenterOnParent();
	dlg.Set( m_grid->GetCellValue( idx, Col_Aperture ) );
	if ( dlg.ShowModal() == wxID_OK )
		m_grid->SetCellValue( idx, Col_Aperture, dlg.Get() );
}

void StageForm::EditSurface(int idx)
{
	if (idx < 0) idx = m_grid->GetGridCursorRow();
	if (idx < 0) return;
	

	ApertureSurfaceDialog dlg( this, ApertureSurfaceDialog::SURFACE, "Edit element " + wxString::Format("%d",(idx+1)) + " surface");
	dlg.CenterOnParent();
	dlg.Set( m_grid->GetCellValue( idx, Col_Surface ) );
	if ( dlg.ShowModal() == wxID_OK )
		m_grid->SetCellValue( idx, Col_Surface, dlg.Get() );
}

void StageForm::EditOptics(int idx)
{
	if ( m_prj.OpticsList.size() == 0 )
	{
		wxMessageBox("No optical properties defined.");
		return;
	}

	if (idx < 0) idx = m_grid->GetGridCursorRow();
	if (idx < 0) return;

	wxArrayString list;
	for (int i=0;i< m_prj.OpticsList.size();i++)
		list.Add( m_prj.OpticsList[i]->Name );

	wxString sel = ::wxGetSingleChoice("Choose optical properties:", "Element " + wxString::Format("%d",idx+1) + " Optics", list, this );
	if (!sel.IsEmpty())
		m_grid->SetCellValue( idx, Col_Optics, sel );
}

void StageForm::OnGridCellDClick(wxGridEvent &evt)
{
	int col = evt.GetCol();
	if (col == Col_Aperture)
		EditAperture( evt.GetRow() );
	else if (col == Col_Surface)
		EditSurface( evt.GetRow() );
	else if (col == Col_Optics)
		EditOptics( evt.GetRow() );
	else
		evt.Skip();
}

void StageForm::Import()
{
	if (!m_stage) return;
	
	wxFileDialog dlg(this, "Import stage geometry");
	if ( dlg.ShowModal()==wxID_OK )
	{
		wxString file = dlg.GetPath();
		FILE *fp = fopen(file.c_str(), "r");
		if (fp)
		{
			Stage *ss = new Stage;
			if ( ss->Read( fp ) )
			{
				// free memory associated with current stage elements
				m_stage->ClearElements();

				// copy over everything, incl. element pointers
				*m_stage = *ss; 

				// erase element pointer list without deleting element data, since pointers
				// to newly read in elements are part of m_stage
				ss->ElementList.clear();

				UpdateFromData();
				m_geoForm->UpdateStageNames();
				Modified();
			}
			else
				wxMessageBox("Error reading stage data from:\n\n" + file );
			
			delete ss;
			fclose( fp );			
		}
		else
			wxMessageBox("Could not open file for reading:\n\n" + file );
	}
}

void StageForm::Export()
{
	if (!m_stage) return;

	wxFileDialog dlg(this, "Export geometry for stage " + m_stageName->GetValue(), wxEmptyString, wxEmptyString, 
		wxFileSelectorDefaultWildcardStr, wxFD_SAVE|wxFD_OVERWRITE_PROMPT);

	if ( wxID_OK == dlg.ShowModal() )
	{
		if ( FILE *fp = fopen(dlg.GetPath().c_str(), "w") )
		{
			m_stage->Write( fp );
			fclose(fp);
		}
		else
			wxMessageBox( "Could not export element data to file:\n\n" + dlg.GetPath() );
	}
}

void StageForm::Clear()
{
	m_grid->DeleteRows(0, m_grid->GetNumberRows() );
}

void StageForm::Delete(int idx)
{
	if (idx >= 0)
	{
		m_grid->DeleteRows(idx, 1);
	}
	else
	{
		if ( !m_stage ) return;

		wxString s = wxGetTextFromUser("Enter rows to delete (e.g. 1,4-5,9), or '*' to delete all:", 
			"Delete elements", wxEmptyString, this);
		
		if (s.IsEmpty()) return;

		if (s.Trim() == "*")
		{
			wxBusyInfo info("Deleting all elements...");
			m_grid->DeleteRows(0, m_stage->ElementList.size());
			return;
		}

		wxArrayString parts = wxStringTokenize(s,",");
		int nremoved = 0;
		for (int i=0;i<(int)parts.Count();i++)
		{
			s = parts[i];
			int hpos = s.Find('-');
			if (hpos<0)
			{
				long row;
				if (s.ToLong(&row) && row-nremoved >= 1 && row-nremoved <= (int)m_stage->ElementList.size())
				{
					m_grid->DeleteRows(row-nremoved-1);
					nremoved++;
				}
				else
					wxMessageBox("Invalid single row specification or value: " + parts[i]);
			}
			else
			{
				long start,end=0;
				if (s.Mid(0,hpos).ToLong(&start)
					&& s.Mid(hpos+1).ToLong(&end)
					&& end-nremoved >= start
					&& start-nremoved >= 1 )
				{
					int ndelete = end-start+1;
					wxBusyInfo *info = NULL;
					if (ndelete > 100)
						info = new wxBusyInfo(wxString::Format("Deleting %d elements...", ndelete));

					m_grid->DeleteRows( start-nremoved-1, ndelete );

					if (info) delete info;

					nremoved += ndelete;
				}
				else
					wxMessageBox("Invalid row specification: " + parts[i]);
			}
		}
	}
}

void StageForm::Insert()
{
	wxString s = wxGetTextFromUser("Enter row position and number of rows to insert (e.g. 3,10):", 
		"Insert elements", wxEmptyString, this );
	if (s.IsEmpty()) return;
	
	wxBusyInfo *info = NULL;
	int hpos = s.Find(',');
	if (hpos < 0)
	{
		long nrows;
		if (s.ToLong(&nrows) && nrows > 0)
		{
			if (nrows > 100) info = new wxBusyInfo(wxString::Format("Inserting %d elements...", nrows));
			m_grid->InsertRows(0, nrows);
		}
		else
			wxMessageBox("Bad row specification for insert: " + s );
	}
	else
	{
		long start, nrows=0;
		if (s.Mid(0,hpos).ToLong(&start)
			&& s.Mid(hpos+1).ToLong(&nrows) 
			&& start >= 1 && nrows > 0)
		{
			if (nrows > 100) info = new wxBusyInfo(wxString::Format("Inserting %d elements...", nrows));
			m_grid->InsertRows(start-1, nrows);
		}
		else
			wxMessageBox("Bad row specification for insert: " + s );
	}

	if (info) delete info;
}

void StageForm::Append(int n)
{
	if (n > 0)
	{
		m_grid->AppendRows( n );
	}
	else
	{
		wxString srows = wxGetTextFromUser("Enter number of rows to append:", "Append Rows", "1", this);

		if (srows.IsEmpty()) return;
		long nrows;
		if (srows.ToLong(&nrows) && nrows > 0)
		{
			wxBusyInfo *info = NULL;
			if (nrows > 100) info = new wxBusyInfo(wxString::Format("Appending %d elements...", nrows));
			m_grid->AppendRows(nrows);
			if (info) delete info;
		}
		else
			wxMessageBox("Invalid number format or number less than 1. No elements appended.");
	}
}

void StageForm::OnCommand( wxCommandEvent &evt )
{
	if ( !m_stage ) return;

	switch( evt.GetId() )
	{
	case ID_STAGE_NAME:
		m_stage->Name = m_stageName->GetValue();
		NameModified();
		return;
		
	case ID_STAGE_ZROT: EditStageZRot(); return;
	case ID_ELEMENT_ZROT: EditZRot(); return;
	case ID_ELEMENT_APERTURE: EditAperture(); return;
	case ID_ELEMENT_SURFACE: EditSurface(); return;
	case ID_ELEMENT_OPTICS: EditOptics(); return;

	case ID_ELEMENT_INSERT: Insert(); return;
	case ID_ELEMENT_APPEND: Append(); return;
	case ID_ELEMENT_DELETE: Delete(); return;
	case ID_ELEMENT_DELETE_ALL:
		if ( wxYES == wxMessageBox("Are you sure you want to delete all elements in the current stage?\n\nThis action is not undo-able.","Query", wxICON_EXCLAMATION|wxYES_NO) )
			Clear();
		return;

	case ID_VIRTUAL_STAGE: m_stage->Virtual = m_virtualStage->GetValue(); break;
	case ID_MULTIPLE_HITS: m_stage->MultiHit = m_multipleHits->GetValue(); break;
	case ID_TRACE_THROUGH: m_stage->TraceThrough = m_traceThrough->GetValue(); break;
	case ID_X: m_stage->X = m_x->Value(); break;
	case ID_Y: m_stage->Y = m_y->Value(); break;
	case ID_Z: m_stage->Z = m_z->Value(); break;
	case ID_AX: m_stage->AX = m_ax->Value(); break;
	case ID_AY: m_stage->AY = m_ay->Value(); break;
	case ID_AZ: m_stage->AZ = m_az->Value(); break;
	case ID_ZROT: m_stage->ZRot = m_zrot->Value(); break;


	case ID_IMPORT: Import(); break;
	case ID_EXPORT: Export(); break;
	case wxID_COPY: m_grid->Copy( true ); break;
	case wxID_PASTE: m_grid->Paste( wxExtGridCtrl::PASTE_ALL_RESIZE_ROWS ); break;
	}

	Modified();
}

void StageForm::Modified()
{
	m_geoForm->Modified();
}

void StageForm::NameModified()
{
	m_geoForm->Modified();
	m_geoForm->UpdateStageNames();
}


enum {
	ID_ADDSTAGE=wxID_HIGHEST+931, ID_INSERTSTAGE, ID_REMOVESTAGE, 
	ID_STAGEBOOK, ID_CLEARALL };

BEGIN_EVENT_TABLE(GeometryForm, wxPanel)
	EVT_BUTTON( ID_ADDSTAGE, GeometryForm::OnButton )
	EVT_BUTTON( ID_INSERTSTAGE, GeometryForm::OnButton )
	EVT_BUTTON( ID_REMOVESTAGE, GeometryForm::OnButton )
	EVT_BUTTON( ID_CLEARALL, GeometryForm::OnButton )
END_EVENT_TABLE()

GeometryForm::GeometryForm(wxWindow *parent, Project &prj)
	: wxPanel(parent), m_prj( prj )
{
	wxBoxSizer *buttonsizer = new wxBoxSizer(wxHORIZONTAL);
	buttonsizer->Add( new wxButton(this, ID_ADDSTAGE, "Add stage..."), 0, wxALL|wxEXPAND, 2);
	buttonsizer->Add( new wxButton(this, ID_INSERTSTAGE, "Insert stage..."), 0, wxALL|wxEXPAND, 2);
	buttonsizer->Add( new wxButton(this, ID_REMOVESTAGE, "Remove stage"), 0, wxALL|wxEXPAND, 2);
	buttonsizer->Add( new wxButton(this, ID_CLEARALL, "Remove all"), 0, wxALL|wxEXPAND, 2);
	buttonsizer->AddStretchSpacer();

	m_stageTabs = new wxNotebook( this,  ID_STAGEBOOK, 
		wxDefaultPosition, wxDefaultSize, wxNB_NOPAGETHEME );

	wxBoxSizer *szvert = new wxBoxSizer(wxVERTICAL);
	szvert->Add( buttonsizer, 0, wxALL|wxEXPAND, 2);
	szvert->Add( m_stageTabs, 1, wxALL|wxEXPAND, 0);
	SetSizer( szvert );
}


void GeometryForm::UpdateStageNames()
{
	for (int i=0;i<m_prj.StageList.size() && i < m_stageTabs->GetPageCount();i++)
		m_stageTabs->SetPageText(i, m_prj.StageList[i]->Name );
}

void GeometryForm::UpdateStageFlagsAndGeom(int stage_num)
{
	if ( StageForm *ss = dynamic_cast<StageForm*>( m_stageTabs->GetPage(stage_num) ) )
		ss->UpdateFromData();
}

void GeometryForm::NewStage( const wxString &name, int pos )
{
	Stage *s = new Stage;
	s->Name = name;
	StageForm *ss = new StageForm( m_stageTabs, this, m_prj, s );

	if ( pos >= 0 )
	{	
		m_prj.StageList.insert( m_prj.StageList.begin() + pos, s);
		m_stageTabs->InsertPage( pos, ss, name, true );
	}
	else
	{
		m_prj.StageList.push_back(s);
		m_stageTabs->AddPage( ss, name, true );
	}

	Modified();
}

StageForm *GeometryForm::GetStageForm(int idx)
{
	if (idx >= 0 && idx < m_stageTabs->GetPageCount()) return dynamic_cast<StageForm*>( m_stageTabs->GetPage(idx) );
	else return NULL;
}

StageForm *GeometryForm::GetStageForm( Stage *stg )
{
	for( size_t i=0;i<m_prj.StageList.size();i++ )
		if ( m_prj.StageList[i] == stg 
			&& i < m_stageTabs->GetPageCount() )
			return dynamic_cast<StageForm*>( m_stageTabs->GetPage(i) );

	return NULL;
}

void GeometryForm::DeleteStage(int sel, bool quiet)
{
	if (sel >= 0 && sel < m_prj.StageList.size())
	{
		if(!quiet && wxMessageBox("Really delete stage: " + m_prj.StageList[sel]->Name ,"Query",wxYES_NO)==wxNO)
			return;
	
		m_stageTabs->DeletePage(sel);
		delete m_prj.StageList[sel];
		m_prj.StageList.erase( m_prj.StageList.begin() + sel );
		Modified();
	}
}

void GeometryForm::ClearStages(bool quiet)
{
	//if (quiet || wxMessageBox("Really delete all stages?","Query",wxYES_NO)==wxYES)
	//{
		m_stageTabs->DeleteAllPages();
		m_prj.ClearStages();
		Modified();
	//}
}


void GeometryForm::UpdateForm()
{
	m_stageTabs->DeleteAllPages();
	for (int i=0;i<m_prj.StageList.size();i++)
	{
		StageForm *ss = new StageForm( m_stageTabs, this, m_prj, m_prj.StageList[i] );
		m_stageTabs->AddPage( ss, m_prj.StageList[i]->Name );
	}
}


void GeometryForm::OnButton(wxCommandEvent &evt)
{
	int sel;
	wxString name;
	switch(evt.GetId())
	{
	case ID_ADDSTAGE:
		NewStage( "New stage" );		
		break;
	case ID_INSERTSTAGE:
		NewStage( "New stage" ,m_stageTabs->GetSelection());		
		break;
	case ID_REMOVESTAGE:
		sel = m_stageTabs->GetSelection();
		if (sel>=0)
			DeleteStage(sel);
		break;
	case ID_CLEARALL:
		if (wxMessageBox("Really remove all stages?", "Query", wxYES_NO,this)==wxYES)
		{
			m_prj.ClearStages();
			UpdateForm();
		}
		break;
	}
}

void GeometryForm::Modified() 
{
	MainWindow::Instance().SetModified( true );
}