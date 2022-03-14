
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


#include <wx/wx.h>
#include <wx/busyinfo.h>
#include <wx/clipbrd.h>
#include <wx/tokenzr.h>
#include <wx/statline.h>
#include <wx/splitter.h>
#include <wx/notebook.h>
#include <wx/sizer.h>
#include <wx/gbsizer.h>
#include <wx/dialog.h>

#include <wex/numeric.h>
#include <wex/exttext.h>
#include <wex/extgrid.h>

#include "optics.h"
#include "soltrace.h"

class AngularInputDialog : public wxDialog
{
	enum { ID_NPOINTS = wxID_HIGHEST + 239, ID_GRID };

	wxExtGridCtrl *m_grid;
	wxNumericCtrl *m_npoints;

public:
	AngularInputDialog( wxWindow *parent, const char* quantity_str = "Reflectivity" )
		: wxDialog( parent, wxID_ANY, "Edit angular dependence table", wxDefaultPosition,
			wxSize( 450, 400 ), wxDEFAULT_DIALOG_STYLE|wxRESIZE_BORDER )
	{
		m_grid = new wxExtGridCtrl( this, ID_GRID );
		m_grid->CreateGrid( 3, 2 );
		m_grid->SetColLabelValue( 0, "Incidence angle\n(mrad)" );
		m_grid->SetColLabelValue(1, wxString::Format("%s \n(0..1)", quantity_str));
		m_grid->AutoSizeColumns();

		m_npoints = new wxNumericCtrl( this, ID_NPOINTS, 3, wxNUMERIC_INTEGER );

		wxBoxSizer *btn_sizer = new wxBoxSizer( wxVERTICAL );
		btn_sizer->Add( new wxStaticText( this, wxID_ANY, "Number of points:") );
		btn_sizer->Add( m_npoints );
		btn_sizer->AddSpacer( 4 );
		btn_sizer->Add( new wxButton( this, wxID_COPY ), 0, wxTOP|wxBOTTOM, 2 );
		btn_sizer->Add( new wxButton( this, wxID_PASTE ), 0, wxTOP|wxBOTTOM, 2 );

		wxBoxSizer *h_sizer = new wxBoxSizer( wxHORIZONTAL );
		h_sizer->Add( m_grid, 1, wxALL|wxEXPAND, 10 );
		h_sizer->Add( btn_sizer, 0, wxALL, 10 );

		wxBoxSizer *sizer = new wxBoxSizer( wxVERTICAL );
		sizer->Add( h_sizer, 1, wxALL|wxEXPAND, 10 );
		sizer->Add( CreateButtonSizer( wxOK|wxCANCEL ), 0, wxALL|wxEXPAND, 10 );
		SetSizer( sizer );
	}

	void SetData( const std::vector<PointF> &data )
	{
		if ( data.size() < 1 ) return;

		m_grid->ResizeGrid( data.size(), 2 );
		for( size_t i=0;i<data.size();i++ )
		{
			m_grid->SetCellValue( i, 0, wxString::Format("%lg", data[i].x ) );
			m_grid->SetCellValue( i, 1, wxString::Format("%lg", data[i].y ) );
		}

		m_npoints->SetValue( data.size() );
	}

	std::vector<PointF> GetData()
	{
		size_t nr = (size_t)m_grid->GetNumberRows();
		std::vector<PointF> data( nr );
		for( size_t i=0;i<nr;i++ )
		{
			data[i].x = wxAtof( m_grid->GetCellValue( i, 0 ) );
			data[i].y = wxAtof( m_grid->GetCellValue( i, 1 ) );
		}
		return data;
	}

	void OnCommand( wxCommandEvent &evt )
	{
		switch( evt.GetId() )
		{
		case ID_NPOINTS:
			if ( m_npoints->AsInteger() > 0 )
			{
				m_grid->ResizeGrid( m_npoints->AsInteger(), 2 );
				SetData( GetData() );
			}
			break;

		case wxID_COPY: m_grid->Copy(); break;
		case wxID_PASTE: m_grid->Paste( wxExtGridCtrl::PASTE_ALL_RESIZE_ROWS ); break;
		}
	}

	void OnCellChange( wxGridEvent &evt )
	{
		int r = evt.GetRow();
		int c = evt.GetCol();

		if ( r < 0 || c < 0 )
		{
			m_npoints->SetValue( m_grid->GetNumberRows() );
			for( int i=0;i<m_grid->GetNumberRows();i++ )
				for( int j=0;j<m_grid->GetNumberCols();j++ )
					m_grid->SetCellValue( i, j, wxString::Format("%lg", wxAtof(m_grid->GetCellValue(i,j))) );
		}
		else
			m_grid->SetCellValue( r, c, wxString::Format("%lg", wxAtof(m_grid->GetCellValue(r,c))) );
	}

	DECLARE_EVENT_TABLE();
};

BEGIN_EVENT_TABLE( AngularInputDialog, wxDialog )
	EVT_BUTTON( wxID_COPY, AngularInputDialog::OnCommand )
	EVT_BUTTON( wxID_PASTE, AngularInputDialog::OnCommand )
	EVT_NUMERIC( AngularInputDialog::ID_NPOINTS, AngularInputDialog::OnCommand )
	EVT_GRID_CMD_CELL_CHANGED( AngularInputDialog::ID_GRID, AngularInputDialog::OnCellChange )
END_EVENT_TABLE()

enum { ID_OPTIC_SURF_NUMBER = wxID_HIGHEST+696,
	ID_APER_STOP_OR_GRATING, ID_DIFFRACTION_ORDER,
	ID_REFRACT_REAL, ID_REFRACT_IMAG,
	ID_GRATING_SPACING, ID_REFLECTIVITY, ID_TRANSMISSIVITY,
	ID_SLOPE_ERROR, ID_SPECULARITY_ERROR, ID_ERROR_TYPE, 
	ID_USE_REFLECTIVITY_TABLE, ID_USE_TRANSMISSIVITY_TABLE, 
	ID_EDIT_REFLECTIVITY_TABLE, ID_EDIT_TRANSMISSIVITY_TABLE };

class OpticalPropertyForm : public wxPanel
{
	wxNumericCtrl *m_opticSurfNumber,
		*m_aperStopOrGrating,
		*m_diffractionOrder,
		*m_refractReal, *m_refractImag,
		*m_gratingSpacing[4],
		*m_reflectivity,
		*m_transmissivity,
		*m_slopeError,
		*m_specularityError;

	wxChoice *m_errorType;
	wxCheckBox* m_useReflectivityTable;
	wxCheckBox *m_useTransmissivityTable;
	wxButton *m_editReflectivityTable;
	wxButton* m_editTransmissivityTable;

	SurfaceOptic *m_surf;
	OpticsForm *m_optForm;


public:
	OpticalPropertyForm( wxWindow *parent ) : wxPanel( parent )
	{
		m_surf = 0;
		m_optForm = 0;

		wxFlexGridSizer *sizer1 = new wxFlexGridSizer( 2 );
		sizer1->AddGrowableCol( 1 );
		sizer1->Add( new wxStaticText( this, wxID_ANY, "Optical surface number" ), 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 2 );
		sizer1->Add( m_opticSurfNumber = new wxNumericCtrl( this, ID_OPTIC_SURF_NUMBER, 1, wxNUMERIC_INTEGER ), 0, wxALL, 2 );
		sizer1->Add( new wxStaticText( this, wxID_ANY, "Aperture stop or grating type" ), 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 2 );
		sizer1->Add( m_aperStopOrGrating = new wxNumericCtrl( this, ID_APER_STOP_OR_GRATING, 3, wxNUMERIC_INTEGER ), 0, wxALL, 2 );
		sizer1->Add( new wxStaticText( this, wxID_ANY, "Diffraction order" ), 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 2 );
		sizer1->Add( m_diffractionOrder = new wxNumericCtrl( this, ID_DIFFRACTION_ORDER, 4, wxNUMERIC_INTEGER ), 0, wxALL, 2 );


		wxStaticBoxSizer *sizer2 = new wxStaticBoxSizer( wxHORIZONTAL, this, "Refraction indicies" );
		sizer2->Add( new wxStaticText( sizer2->GetStaticBox(), wxID_ANY, "Real" ), 0, wxALL|wxALIGN_CENTER_VERTICAL, 4 );
		sizer2->Add( m_refractReal = new wxNumericCtrl( sizer2->GetStaticBox(), ID_REFRACT_REAL, 1.1 ), 0, wxALL, 4 );
		sizer2->Add( new wxStaticText( sizer2->GetStaticBox(), wxID_ANY, "Imag" ), 0, wxALL|wxALIGN_CENTER_VERTICAL, 4 );
		sizer2->Add( m_refractImag = new wxNumericCtrl( sizer2->GetStaticBox(), ID_REFRACT_IMAG, 1.2 ), 0, wxALL, 4 );

		wxStaticBoxSizer *sizer3 = new wxStaticBoxSizer( wxVERTICAL, this, "Grating spacing coefficients" );		
		const char *labels[4] = { "1st", "2nd", "3rd", "4th" };
		wxFlexGridSizer *sizer4 = new wxFlexGridSizer(2 );
		sizer4->AddGrowableCol(1);
		for( int i=0;i<4;i++ )
		{
			sizer4->Add( new wxStaticText( sizer3->GetStaticBox(), wxID_ANY, labels[i] ), 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 2 );
			sizer4->Add( m_gratingSpacing[i] = new wxNumericCtrl( sizer3->GetStaticBox(), ID_GRATING_SPACING, 1.1 + i*0.1 ) , 0, wxALL, 2 );
			m_gratingSpacing[i]->Enable( false );
		}
		sizer3->Add( sizer4, 0, wxALL, 5 );

		m_opticSurfNumber->Enable( false );
		m_aperStopOrGrating->Enable( false );
		m_diffractionOrder->Enable( false );
		m_refractImag->Enable( false );


		wxBoxSizer *r_sizer = new wxBoxSizer( wxVERTICAL );
		r_sizer->Add( sizer1, 0, wxALL, 5 );
		r_sizer->Add( sizer2, 0, wxALL, 5 );
		r_sizer->Add( sizer3, 0, wxALL, 5 );


		wxFlexGridSizer *l_sizer = new wxFlexGridSizer( 4 );
		l_sizer->AddGrowableCol( 1 );
		
		l_sizer->Add( new wxStaticText( this, wxID_ANY, "Reflectivity" ), 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 2 );
		l_sizer->Add( m_reflectivity = new wxNumericCtrl( this, ID_REFLECTIVITY, 0.96 ), 0, wxALL, 2 );
		l_sizer->Add(m_useReflectivityTable = new wxCheckBox(this, ID_USE_REFLECTIVITY_TABLE, L"\u03c1(\u03b8)"), 0, wxLEFT | wxALIGN_CENTER_VERTICAL | wxALIGN_RIGHT, 25);
		l_sizer->Add(m_editReflectivityTable = new wxButton(this, ID_EDIT_REFLECTIVITY_TABLE, "...", wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT), 0, wxALL | wxALIGN_CENTER_VERTICAL, 0);
		
		l_sizer->Add( new wxStaticText( this, wxID_ANY, "Transmissivity" ), 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 2 );
		l_sizer->Add( m_transmissivity = new wxNumericCtrl( this, ID_TRANSMISSIVITY, 1.0 ), 0, wxALL, 2 );
		l_sizer->Add(m_useTransmissivityTable= new wxCheckBox( this, ID_USE_TRANSMISSIVITY_TABLE, L"\u03c4(\u03b8)" ), 0, wxLEFT|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 25 );
		l_sizer->Add(m_editTransmissivityTable = new wxButton( this, ID_EDIT_TRANSMISSIVITY_TABLE, "...", wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT ), 0, wxALL|wxALIGN_CENTER_VERTICAL, 0 );
		
		l_sizer->Add( new wxStaticText( this, wxID_ANY, "Slope error" ), 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 2 );
		l_sizer->Add( m_slopeError = new wxNumericCtrl( this, ID_SLOPE_ERROR, 0.95 ), 0, wxALL, 2 );
		l_sizer->Add( new wxStaticText( this, wxID_ANY, "mrad" ), 0, wxALL|wxALIGN_CENTER_VERTICAL, 2 );
		l_sizer->AddStretchSpacer();
		
		l_sizer->Add( new wxStaticText( this, wxID_ANY, "Specularity error" ), 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 2 );
		l_sizer->Add( m_specularityError = new wxNumericCtrl( this, ID_SPECULARITY_ERROR, 0.2 ), 0, wxALL, 2 );
		l_sizer->Add( new wxStaticText( this, wxID_ANY, "mrad" ), 0, wxALL|wxALIGN_CENTER_VERTICAL, 2 );
		l_sizer->AddStretchSpacer();
		
		l_sizer->Add( new wxStaticText( this, wxID_ANY, "Error type" ), 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 2 );
		l_sizer->Add( m_errorType = new wxChoice( this, ID_ERROR_TYPE ), 0, wxALL, 2 );
		l_sizer->AddStretchSpacer();
		l_sizer->AddStretchSpacer();
		m_errorType->Append( "Gaussian" );
		m_errorType->Append( "Pillbox" );
		m_errorType->Append( "Diffuse" );
		m_errorType->SetSelection(0);
		
		wxBoxSizer *sizer = new wxBoxSizer( wxHORIZONTAL );
		sizer->Add( l_sizer, 0, wxALL, 5 );
		sizer->Add( r_sizer, 0, wxALL, 5 );
		SetSizer( sizer );
	}

	void InitForm( OpticsForm *o, SurfaceOptic *surf )
	{
		if ( !o || !surf ) {
			m_optForm = 0;
			m_surf = 0;
			Enable( false );
			return;
		}

		m_surf = surf;
		m_optForm = o;

		m_reflectivity->Enable( !surf->UseReflectivityTable );
		m_editReflectivityTable->Enable( surf->UseReflectivityTable );
		m_transmissivity->Enable(!surf->UseTransmissivityTable);
		m_editTransmissivityTable->Enable(surf->UseTransmissivityTable);
		m_opticSurfNumber->SetValue( surf->OpticalSurfaceNumber );
		m_aperStopOrGrating->SetValue( surf->ApertureStopOrGratingType );
		m_diffractionOrder->SetValue( surf->DiffractionOrder );
		m_refractReal->SetValue( surf->RefractionIndexReal );
		m_refractImag->SetValue( surf->RefractionIndexImag );
		for( int i=0;i<4;i++ )
			m_gratingSpacing[i]->SetValue( surf->GratingCoeffs[i] );
		m_reflectivity->SetValue( surf->Reflectivity );
		m_transmissivity->SetValue( surf->Transmissivity );
		m_slopeError->SetValue( surf->RMSSlope );
		m_specularityError->SetValue( surf->RMSSpecularity );
		//m_errorType->SetSelection( surf->ErrorDistribution=='g' ? 0 : 1 );
		if (surf->ErrorDistribution == 'g')
		{
			m_errorType->SetSelection(0);
		}
		if (surf->ErrorDistribution == 'p')
		{
			m_errorType->SetSelection(1);
		}
		if (surf->ErrorDistribution == 'f')
		{
			m_errorType->SetSelection(2);
		}

		m_useReflectivityTable->SetValue(surf->UseReflectivityTable);
		m_useTransmissivityTable->SetValue(surf->UseTransmissivityTable);

		Enable( true );
	}

	void Modified() {
		if ( m_optForm ) m_optForm->Modified();
	}

	void HandleEvent( wxCommandEvent &evt )
	{
		int i;
		if ( !m_surf ) return;

		switch( evt.GetId() )
		{
		case ID_OPTIC_SURF_NUMBER: m_surf->OpticalSurfaceNumber = m_opticSurfNumber->AsInteger(); break;
		case ID_APER_STOP_OR_GRATING: m_surf->ApertureStopOrGratingType = m_aperStopOrGrating->AsInteger(); break;
		case ID_DIFFRACTION_ORDER: m_surf->DiffractionOrder = m_diffractionOrder->AsInteger(); break;
		case ID_REFRACT_REAL: m_surf->RefractionIndexReal = m_refractReal->Value(); break;
		case ID_REFRACT_IMAG: m_surf->RefractionIndexImag = m_refractImag->Value(); break;
		case ID_GRATING_SPACING:
			for( i=0;i<4;i++ )
				m_surf->GratingCoeffs[i] = m_gratingSpacing[i]->Value();
			break;
		case ID_REFLECTIVITY: m_surf->Reflectivity = m_reflectivity->Value(); break;
		case ID_TRANSMISSIVITY: m_surf->Transmissivity = m_transmissivity->Value(); break;
		case ID_SLOPE_ERROR: m_surf->RMSSlope = m_slopeError->Value(); break;
		case ID_SPECULARITY_ERROR: m_surf->RMSSpecularity = m_specularityError->Value(); break;
		case ID_ERROR_TYPE:
			//m_surf->ErrorDistribution = (m_errorType->GetSelection()==0) ? 'g' : 'p';
			if (m_errorType->GetSelection() == 0)
			{
				m_surf->ErrorDistribution = 'g';
			}
			if (m_errorType->GetSelection() == 1)
			{
				m_surf->ErrorDistribution = 'p';
			}
			if (m_errorType->GetSelection() == 2)
			{
				m_surf->ErrorDistribution = 'f';
			}			
			break;
		case ID_USE_REFLECTIVITY_TABLE:
			m_surf->UseReflectivityTable = m_useReflectivityTable->GetValue();
			m_editReflectivityTable->Enable( m_useReflectivityTable->GetValue() );
			m_reflectivity->Enable( !m_useReflectivityTable->GetValue() );
			break;
		case ID_USE_TRANSMISSIVITY_TABLE:
			m_surf->UseTransmissivityTable = m_useTransmissivityTable->GetValue();
			m_editTransmissivityTable->Enable(m_useTransmissivityTable->GetValue());
			m_transmissivity->Enable(!m_useTransmissivityTable->GetValue());
			break;
		case ID_EDIT_REFLECTIVITY_TABLE:
			{
				AngularInputDialog dlg( this, "Reflectivity");
				dlg.CenterOnParent();
				dlg.SetData( m_surf->ReflectivityTable );
				if (  wxID_OK == dlg.ShowModal() )
					m_surf->ReflectivityTable = dlg.GetData();
				else
					return; // skip modifiying the project if refl. table not changed.
			}
			break;
		case ID_EDIT_TRANSMISSIVITY_TABLE:
		{
			AngularInputDialog dlg(this, "Transmissivity");
			dlg.CenterOnParent();
			dlg.SetData(m_surf->TransmissivityTable);
			if (wxID_OK == dlg.ShowModal())
				m_surf->TransmissivityTable = dlg.GetData();
			else
				return; // skip modifiying the project if refl. table not changed.
		}
		break;
		}

		Modified();
	}

	DECLARE_EVENT_TABLE();
};

BEGIN_EVENT_TABLE( OpticalPropertyForm, wxPanel )
	EVT_NUMERIC( ID_OPTIC_SURF_NUMBER, OpticalPropertyForm::HandleEvent )
	EVT_NUMERIC( ID_APER_STOP_OR_GRATING, OpticalPropertyForm::HandleEvent )
	EVT_NUMERIC( ID_DIFFRACTION_ORDER, OpticalPropertyForm::HandleEvent )
	EVT_NUMERIC( ID_REFRACT_REAL, OpticalPropertyForm::HandleEvent )
	EVT_NUMERIC( ID_REFRACT_IMAG, OpticalPropertyForm::HandleEvent )
	EVT_NUMERIC( ID_GRATING_SPACING, OpticalPropertyForm::HandleEvent )
	EVT_NUMERIC( ID_REFLECTIVITY, OpticalPropertyForm::HandleEvent )
	EVT_NUMERIC( ID_TRANSMISSIVITY, OpticalPropertyForm::HandleEvent )
	EVT_NUMERIC( ID_SLOPE_ERROR, OpticalPropertyForm::HandleEvent )
	EVT_NUMERIC( ID_SPECULARITY_ERROR, OpticalPropertyForm::HandleEvent )
	EVT_CHOICE( ID_ERROR_TYPE, OpticalPropertyForm::HandleEvent )
	EVT_CHECKBOX(ID_USE_REFLECTIVITY_TABLE, OpticalPropertyForm::HandleEvent)
	EVT_BUTTON( ID_EDIT_REFLECTIVITY_TABLE, OpticalPropertyForm::HandleEvent )
	EVT_CHECKBOX( ID_USE_TRANSMISSIVITY_TABLE, OpticalPropertyForm::HandleEvent )
	EVT_BUTTON(ID_EDIT_TRANSMISSIVITY_TABLE, OpticalPropertyForm::HandleEvent)
END_EVENT_TABLE()

enum { ID_ADDOPTIC= wxID_HIGHEST+239, ID_REMOVEOPTIC, ID_OPTNAME, ID_OPTICLIST, ID_IMPORT, ID_EXPORT, ID_CLEARALL };

BEGIN_EVENT_TABLE(OpticsForm, wxPanel)
	EVT_BUTTON( ID_ADDOPTIC, OpticsForm::OnButton )
	EVT_BUTTON( ID_REMOVEOPTIC, OpticsForm::OnButton )
	EVT_BUTTON( ID_IMPORT, OpticsForm::OnButton )
	EVT_BUTTON( ID_EXPORT, OpticsForm::OnButton )
	EVT_BUTTON( ID_CLEARALL, OpticsForm::OnButton )
	EVT_TEXT_ENTER( ID_OPTNAME, OpticsForm::OnNameChange )
	EVT_LISTBOX( ID_OPTICLIST, OpticsForm::OnListSelect)
END_EVENT_TABLE()

OpticsForm::OpticsForm( wxWindow *parent, Project &p)
	: wxPanel( parent ), m_prj(p)
{
	wxBoxSizer *btn_sizer = new wxBoxSizer(wxHORIZONTAL);
	btn_sizer->Add( new wxButton(this, ID_ADDOPTIC, "Add optical properties..."), 0, wxALL|wxEXPAND, 2);
	btn_sizer->Add( new wxButton(this, ID_REMOVEOPTIC, "Remove"), 0, wxALL|wxEXPAND, 2);
	btn_sizer->Add( new wxButton(this, ID_CLEARALL, "Remove all"), 0, wxALL|wxEXPAND, 2);
	btn_sizer->Add( new wxButton(this, ID_IMPORT, "Import..."), 0, wxALL|wxEXPAND, 2);
	btn_sizer->Add( new wxButton(this, ID_EXPORT, "Export..."), 0, wxALL|wxEXPAND, 2);
	btn_sizer->AddStretchSpacer();
	
	m_opticList = new wxListBox(this, ID_OPTICLIST);

	wxPanel *tabpanel = new wxPanel(this );
	m_optName = new wxExtTextCtrl( tabpanel, ID_OPTNAME );
	m_optTabs = new wxNotebook( tabpanel, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxNB_NOPAGETHEME );
	
	wxBoxSizer *tabszh = new wxBoxSizer(wxHORIZONTAL);
	tabszh->Add(new wxStaticText(tabpanel, wxID_ANY, "Name:"), 0, wxLEFT|wxRIGHT|wxALIGN_CENTER_VERTICAL, 10 );
	tabszh->Add( m_optName );

	wxBoxSizer *tabsz = new wxBoxSizer(wxVERTICAL);
	tabsz->Add( tabszh, 0, wxALL, 5 );
	tabsz->Add( m_optTabs, 0, wxALL, 0);
	tabpanel->SetSizer(tabsz);
	
	m_frontOpt = new OpticalPropertyForm( m_optTabs );
	m_optTabs->AddPage( m_frontOpt, "Front side" );

	m_backOpt = new OpticalPropertyForm( m_optTabs );
	m_optTabs->AddPage( m_backOpt, "Back side" );

	m_optName->Enable(false);
	m_frontOpt->Enable(false);
	m_backOpt->Enable(false);
	

	wxBoxSizer *h_sizer = new wxBoxSizer( wxHORIZONTAL );
	h_sizer->Add( m_opticList, 0, wxALL|wxEXPAND, 3 );
	h_sizer->Add( tabpanel, 1, wxALL, 3 );

	wxBoxSizer *szvert = new wxBoxSizer(wxVERTICAL);
	szvert->Add( btn_sizer, 0, wxALL|wxEXPAND, 0);
	szvert->Add( h_sizer, 1, wxALL|wxEXPAND, 0);

	SetSizer( szvert );
}

void OpticsForm::UpdateList(int sel)
{
	m_opticList->Clear();
	std::vector<Optical*> list = m_prj.OpticsList;
	for (size_t i=0;i<list.size();i++)
		m_opticList->Append( list[i]->Name );

	if (sel >= 0 && sel < m_opticList->GetCount())
	{
		m_opticList->SetSelection(sel);
		UpdateOptForms();
	}
}

void OpticsForm::Modified()
{
	MainWindow::Instance().SetModified( true );
}

void OpticsForm::UpdateOptForms(int idx)
{
	if (idx == -1) idx = m_opticList->GetSelection();
	if (idx >= 0 && idx < m_prj.OpticsList.size())
	{
		m_optName->SetValue( m_prj.OpticsList[idx]->Name );
		m_frontOpt->InitForm( this, &m_prj.OpticsList[idx]->Front );
		m_backOpt->InitForm( this, &m_prj.OpticsList[idx]->Back );

		m_frontOpt->Enable(true);
		m_backOpt->Enable(true);
		m_optName->Enable(true);
	}
	else
	{
		m_frontOpt->InitForm(NULL, NULL);
		m_frontOpt->Enable(false);
		m_backOpt->InitForm(NULL, NULL);
		m_backOpt->Enable(false);
		m_optName->SetValue(wxEmptyString);
		m_optName->Enable(false);
	}
}

void OpticsForm::AddOptic(const wxString &name)
{
	Optical *opt = new Optical;
	opt->Name = name;
	m_prj.OpticsList.push_back(opt);
	UpdateList( m_prj.OpticsList.size()-1 );

	Modified();
}

void OpticsForm::DeleteOptic(int sel)
{
	if (sel >= 0 && sel < m_prj.OpticsList.size())
	{
		m_frontOpt->InitForm(NULL, NULL);
		m_backOpt->InitForm(NULL, NULL);
		delete m_prj.OpticsList[sel];
		m_prj.OpticsList.erase( m_prj.OpticsList.begin() + sel );

		if (sel-1 >= 0)
			UpdateList( sel-1 );
		else if ( m_prj.OpticsList.size() > 0)
			UpdateList( 0 );
		else
		{
			UpdateList();
			UpdateOptForms(-2);
		}
		
		Modified();
	}
}

void OpticsForm::ClearOptics()
{
	m_prj.ClearOptics();
	UpdateList();
	UpdateOptForms(-2);
	Modified();
}

void OpticsForm::OnButton(wxCommandEvent &evt)
{
	switch(evt.GetId())
	{
	case ID_ADDOPTIC:
		AddOptic();
		break;
	case ID_REMOVEOPTIC:
		{
			int sel = m_opticList->GetSelection();
			DeleteOptic(sel);
		}
		break;
	case ID_IMPORT:
		{
			wxFileDialog dlg(this, "Import optical properties file");
			if (dlg.ShowModal() != wxID_OK) return;
			FILE *fp = fopen(dlg.GetPath().c_str(), "r");
			if (!fp)
			{
				wxMessageBox("Could not open the specified file for reading.");
				return;
			}

			Optical *opt = new Optical;
			if (opt->Read( fp, wxFileNameFromPath(dlg.GetPath())))
			{
				if (opt->Name.IsEmpty())
					opt->Name = wxFileNameFromPath(dlg.GetPath());

				m_prj.OpticsList.push_back(opt);
				UpdateList( m_prj.OpticsList.size()-1 );
				Modified();
			}
			else
			{
				delete opt;
				wxMessageBox("Failed to read optical properties file format.");
			}

			fclose(fp);
		}
		break;
	case ID_EXPORT:
		{
			int sel = m_opticList->GetSelection();
			if (sel >= 0 && sel < m_prj.OpticsList.size())
			{
				wxFileDialog dlg(this, "Export optical properties file", wxEmptyString, wxEmptyString, wxFileSelectorDefaultWildcardStr, wxFD_SAVE|wxFD_OVERWRITE_PROMPT);
				if (dlg.ShowModal() != wxID_OK) return;

				FILE *fp = fopen(dlg.GetPath().c_str(), "w");
				if (!fp)
				{
					wxMessageBox("Could not open file for writing.");
					return;
				}

				m_prj.OpticsList[sel]->Write(fp);
				fclose(fp);
			}
			else
			{
				wxMessageBox("No optical properties selected.");
			}
		}
		break;
	case ID_CLEARALL:
		{
			ClearOptics();
		}
		break;
	}
}

void OpticsForm::OnNameChange(wxCommandEvent &evt)
{
	int sel = m_opticList->GetSelection();
	if (sel >= 0 && sel < m_prj.OpticsList.size())
	{
		if ( !m_optName->GetValue().IsEmpty())
		{
			m_prj.OpticsList[sel]->Name = m_optName->GetValue();
			UpdateList( m_opticList->GetSelection() );
		}
		else
			m_optName->SetValue( m_prj.OpticsList[sel]->Name );
		
		Modified();
	}
}

void OpticsForm::OnListSelect(wxCommandEvent &evt)
{
	UpdateOptForms();
}
