
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
#include <wx/progdlg.h>

#include <wex/extgrid.h>

#include <stapi.h>
#include "raydata.h"

class RayDataTable : public wxGridTableBase
{
	Project *m_prj;
	std::vector<size_t> m_indices;
	int m_coordSys;

public:
	static const int NCOLS = 9;

	RayDataTable( Project *prj = 0, int stage = -1, int coords = RayData::COORD_GLOBAL )
	{
		SetAttrProvider( new wxExtGridCellAttrProvider );

		m_prj = prj;
		m_indices.clear();
		m_coordSys = coords;

		if ( !m_prj ) return;

		RayData &r = m_prj->Results;

		size_t nindices = 0;
		for (size_t i = 0; i < r.Length; i++ )
			if (stage < 0 || stage == r.StageMap[i]-1)
				nindices++;

		try {
			m_indices.resize( nindices );
		} catch( std::exception &e )
		{
			wxMessageBox("Exception during generation of results data table: " + wxString(e.what()));
			m_indices.clear();
			return;
		}

		size_t j=0;
		for (size_t i=0; i < r.Length; i++)
			if (stage < 0 || stage == r.StageMap[i]-1)
				m_indices[j++] = i;
	}

	virtual ~RayDataTable()
	{
		m_prj = 0;
	}
	

	int GetNumberRows()
	{
		return m_indices.size();
	}

	int GetNumberCols()
	{
		return NCOLS;
	}

	bool IsEmptyCell( int WXUNUSED(row), int WXUNUSED(col) )
	{
		return false;
	}

	wxString GetValue( int row, int col )
	{
		wxString val;
		size_t nindices = m_indices.size();

		double origin[3], posin[3], cosin[3], posout[3], cosout[3];
		if ( 0 == m_prj || nindices == 0 || ((size_t)row) >= nindices ) return val;

		size_t i = m_indices[ (size_t)row ];
		if (i >= m_prj->Results.Length) return val;

		posout[0] = posin[0] = m_prj->Results.Xi[i];
		posout[1] = posin[1] = m_prj->Results.Yi[i];
		posout[2] = posin[2] = m_prj->Results.Zi[i];

		cosout[0] = cosin[0] = m_prj->Results.Xc[i];
		cosout[1] = cosin[1] = m_prj->Results.Yc[i];
		cosout[2] = cosin[2] = m_prj->Results.Zc[i];

		if ( RayData::COORD_GLOBAL == m_coordSys )
		{
			Stage *stage = m_prj->GetStage( m_prj->Results.StageMap[i] - 1 );
			if (!stage) return "<error>";

			origin[0] = stage->X;
			origin[1] = stage->Y;
			origin[2] = stage->Z;

			::st_transform_to_reference( posin, cosin,
											origin, stage->RLocToRef,
											posout, cosout );
		}
		else if ( RayData::COORD_ELEMENT == m_coordSys )
		{
			if ( Element *element = m_prj->GetElement( m_prj->Results.StageMap[i]-1,
														abs(m_prj->Results.ElementMap[i])-1 ) )
			{
				origin[0] = element->X;
				origin[1] = element->Y;
				origin[2] = element->Z;

				::st_transform_to_local( posin, cosin,
											origin, element->RRefToLoc,
											posout, cosout );
			}
			else if ( col < 7 ) return "missed"; // still display stage, and ray numbers
		}

		switch( col )
		{
		case 0: val.Printf("%lg", posout[0]); break;
		case 1: val.Printf("%lg", posout[1]); break;
		case 2: val.Printf("%lg", posout[2]); break;
		case 3: val.Printf("%lg", cosout[0]); break;
		case 4: val.Printf("%lg", cosout[1]); break;
		case 5: val.Printf("%lg", cosout[2]); break;
		case 6: val.Printf("%d", m_prj->Results.ElementMap[i] ); break;
		case 7: val.Printf("%d", m_prj->Results.StageMap[i] ); break;
		case 8: val.Printf("%d", m_prj->Results.RayNumbers[i] ); break;
		}
		
		return val;
	}

	void SetValue( int WXUNUSED(row), int WXUNUSED(col), const wxString &)
	{
		/* nothing to do */
	}

	wxString GetColLabelValue( int col )
	{
		switch (col)
		{
		case 0: return "Pos X";
		case 1: return "Pos Y";
		case 2: return "Pos Z";
		case 3: return "Cos X";
		case 4: return "Cos Y";
		case 5: return "Cos Z";
		case 6: return "Element";
		case 7: return "Stage";
		case 8: return "Ray";
		default:
			return wxEmptyString;
		}
	}

	wxString GetTypeName( int WXUNUSED(row), int WXUNUSED(col) )
	{
		return wxGRID_VALUE_STRING;
	}
};

enum { ID_COORDINATES = wxID_HIGHEST + 931,
	ID_STAGE };

BEGIN_EVENT_TABLE( RayDataForm, wxPanel )
	EVT_BUTTON( wxID_COPY, RayDataForm::OnCommand )
	EVT_BUTTON( wxID_SAVE, RayDataForm::OnCommand )
	EVT_CHOICE( ID_COORDINATES, RayDataForm::OnCommand )
	EVT_CHOICE( ID_STAGE, RayDataForm::OnCommand )
END_EVENT_TABLE()

RayDataForm::RayDataForm( wxWindow *parent, Project &prj )
	: wxPanel( parent ), m_prj(prj)
{
	m_coords = new wxChoice( this, ID_COORDINATES );
	m_coords->Append( "Global" );
	m_coords->Append( "Stage" );
	m_coords->Append( "Element" );
	m_coords->SetSelection( 0 );
	
	m_stage = new wxChoice( this, ID_STAGE );
	m_stage->Append( "All stages" );
	m_stage->SetSelection( 0 );

	wxBoxSizer *tool_sizer = new wxBoxSizer(wxHORIZONTAL);
	tool_sizer->Add( new wxStaticText( this, wxID_ANY, "Coordinate system:" ), 0, wxLEFT|wxALIGN_CENTER_VERTICAL, 10 );
	tool_sizer->Add( m_coords, 0, wxLEFT|wxRIGHT|wxALIGN_CENTER_VERTICAL, 5 );
	tool_sizer->Add( new wxStaticText( this, wxID_ANY, "Stage:" ), 0, wxLEFT|wxALIGN_CENTER_VERTICAL, 10 );
	tool_sizer->Add( m_stage, 0, wxLEFT|wxRIGHT|wxALIGN_CENTER_VERTICAL, 5 );
	tool_sizer->Add( new wxButton(this, wxID_COPY, "Copy to clipboard"), 0, wxALL|wxEXPAND, 2);
	tool_sizer->Add( new wxButton(this, wxID_SAVE, "Save as CSV..."), 0, wxALL|wxEXPAND, 2);
	tool_sizer->AddStretchSpacer();

	m_grid = new wxExtGridCtrl( this, wxID_ANY );
	m_grid->DisableDragRowSize();
	m_grid->DisableDragColMove();
	m_grid->DisableDragGridSize();
	m_grid->SetRowLabelSize(23);
	m_grid->SetColLabelSize(23);
	m_grid->SetDefaultCellAlignment( wxALIGN_RIGHT, wxALIGN_CENTER );
	m_grid->SetRowLabelAlignment( wxALIGN_LEFT, wxALIGN_CENTER );
	
	m_grid->SetTable( new RayDataTable(), true );

	wxBoxSizer *sizer = new wxBoxSizer( wxVERTICAL );
	sizer->Add( tool_sizer, 0, wxALL|wxEXPAND, 0 );
	sizer->Add( m_grid, 1, wxALL|wxEXPAND, 0 );
	SetSizer(sizer);

	UpdateDataDisplay();
}

void RayDataForm::UpdateDataDisplay()
{
	int coord = RayData::COORD_GLOBAL;
	if ( m_coords->GetSelection() == 1 ) coord = RayData::COORD_STAGE;
	else if ( m_coords->GetSelection() == 2 ) coord = RayData::COORD_ELEMENT;

	int istage = m_stage->GetSelection();
	m_stage->Clear();
	m_stage->Append( "All stages" );
	for( size_t i=0;i<m_prj.StageList.size();i++ )
		m_stage->Append( m_prj.StageList[i]->Name );

	if ( istage < 0 ) istage = 0;
	if ( istage >= (int)m_stage->GetCount() ) istage = m_stage->GetCount()-1;

	m_stage->SetSelection( istage );

	Layout();

	m_grid->SetTable( new RayDataTable( &m_prj, istage-1, coord ), true );
	
	for (int i=0;i<RayDataTable::NCOLS;i++)
	{
		m_grid->SetColMinimalWidth(i, 100);
		m_grid->SetColSize(i, 100);
	}

	m_grid->SetRowLabelSize( 70 );
	m_grid->ForceRefresh();
}

void RayDataForm::CopyDataToClipboard()
{
	size_t nrecords = m_prj.Results.PrepareExport( m_coords->GetSelection(), m_stage->GetSelection() );

	size_t bytes_estimate = 24 * RayDataTable::NCOLS * (nrecords+1);
	if ( bytes_estimate > 25e6 && wxNO == wxMessageBox( "This will place more than 25 MB of data on the system clipboard.  "
								"You are recommended to export the data as CSV or binary format instead.\n\nContinue anyway?", 
								"Query", wxYES_NO ) )
			return;

	
	wxBusyInfo info( "Copying to clipboard...", this );
	wxYield();

	wxString text;
	try {
		// preallocate a big string
		text.Alloc( bytes_estimate );
	} catch( std::exception &ex ) {
		wxMessageBox("Error copying to the clipboard.  The data is probably too large.\n\n" + wxString(ex.what()) );
		return;
	}

	double Pos[3], Cos[3];
	int Elm, Stg, Ray;
	
	text = "Pos X\tPos Y\tPos Z\tCos X\tCos Y\tCos Z\tElement\tStage\tRay Number\n";
	while( m_prj.Results.GetNextExport( m_prj, Pos, Cos, Elm, Stg, Ray ) )
	{
		text += wxString::Format("%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%d\t%d\t%d\n",
			Pos[0], Pos[1], Pos[2],
			Cos[0], Cos[1], Cos[2],
			Elm, Stg, Ray );
	}
		
	if (wxTheClipboard->Open())
	{
		wxTheClipboard->SetData( new wxTextDataObject(text) );
		wxTheClipboard->Close();
	}
}

void RayDataForm::SaveDataAsCSV()
{
	wxFileDialog dlg( this, "Save ray data as CSV", wxEmptyString, 
		"raydata.csv", "CSV Files (*.csv)|*.csv", wxFD_SAVE|wxFD_OVERWRITE_PROMPT);

	if (dlg.ShowModal() != wxID_OK)
		return;

	FILE *fp = fopen( dlg.GetPath().c_str(), "w" );
	if ( !fp )
	{
		wxMessageBox("Could not open file for writing:\n\n" + dlg.GetPath());
		return;
	}
	
	wxProgressDialog pd( "CSV data file export", "Preparing...", 100, 
		this, wxPD_APP_MODAL|wxPD_REMAINING_TIME|wxPD_SMOOTH|wxPD_CAN_ABORT|wxPD_AUTO_HIDE );
#ifdef __WXMSW__
	pd.SetIcon( wxICON( appicon ) );
#endif
	pd.Show();

	wxYield();


	// by default, all ray data is in stage coordinates
	RayData &rd = m_prj.Results;
	
	int nwrite = rd.PrepareExport( m_coords->GetSelection(), m_stage->GetSelection() );
	
	// progress update every 0.5 %
	int nupdate = nwrite / 200;
	size_t nwr = 0;
	size_t bytes = 0;
	
	double Pos[3], Cos[3];
	int Elm, Stg, Ray;

	fputs( "Pos X,Pos Y,Pos Z,Cos X,Cos Y,Cos Z,Element,Stage,Ray Number\n", fp );
	while( rd.GetNextExport( m_prj, Pos, Cos, Elm, Stg, Ray ) )
	{
		if ( nwr % nupdate == 0 )
		{
			bool proceed = pd.Update( (int)( 100*((double)nwr)/((double)nwrite) ), 
					"Writing data to " + dlg.GetPath() + wxString::Format(" (%.2lf MB)", bytes*0.000001 ) );
			if( !proceed )
				break;
		}

		int nb = fprintf( fp, "%lg,%lg,%lg,%lg,%lg,%lg,%d,%d,%d\n",
			Pos[0], Pos[1], Pos[2],
			Cos[0], Cos[1], Cos[2],
			Elm, Stg, Ray );

		if ( nb < 0 ) break; // file error, disk space issue?

		bytes += nb;
		nwr++;
	}

	if ( ferror( fp ) )
		wxMessageBox(wxString("An error occurred exporting the CSV data file: ") + strerror(errno));

	fclose( fp );
}

void RayDataForm::OnCommand( wxCommandEvent &evt )
{
	switch( evt.GetId() )
	{
	case ID_COORDINATES:
	case ID_STAGE:
		UpdateDataDisplay();
		break;

	case wxID_SAVE:
		SaveDataAsCSV();
		break;

	case wxID_COPY:
		CopyDataToClipboard();
		break;
	}
}
