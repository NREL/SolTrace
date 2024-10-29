
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
#include <wx/filedlg.h>

#include <wex/numeric.h>
#include <wex/plot/plplotctrl.h>
#include <wex/plot/plcontourplot.h>
#include <wex/plot/plcolourmap.h>

#include "elementlist.h"
#include "fluxmap.h"
#include "soltrace.h"

enum{ ID_ELEMENT_LIST = wxID_HIGHEST+198, 
	ID_XBINS, ID_YBINS,
	ID_AUTO_EXTENT,
	ID_MIN_X, ID_MAX_X,
	ID_MIN_Y, ID_MAX_Y,
	ID_DNI, ID_FINAL_ONLY,
	ID_CONTOUR_LEVELS,
	ID_COLOR_SCHEME,
	ID_EXPORT_TECPLOT
};

BEGIN_EVENT_TABLE( FluxMapForm, wxPanel )
	EVT_NUMERIC( ID_XBINS, FluxMapForm::OnCommand )
	EVT_NUMERIC( ID_YBINS, FluxMapForm::OnCommand )
	EVT_LISTBOX( ID_ELEMENT_LIST, FluxMapForm::OnCommand )
	EVT_CHECKBOX( ID_AUTO_EXTENT, FluxMapForm::OnCommand )
	EVT_NUMERIC( ID_MIN_X, FluxMapForm::OnCommand ) 
	EVT_NUMERIC( ID_MAX_X, FluxMapForm::OnCommand )
	EVT_NUMERIC( ID_MIN_Y, FluxMapForm::OnCommand )
	EVT_NUMERIC( ID_MAX_Y, FluxMapForm::OnCommand )
	EVT_NUMERIC( ID_DNI, FluxMapForm::OnCommand )
	EVT_NUMERIC( ID_CONTOUR_LEVELS, FluxMapForm::OnCommand )
	EVT_CHECKBOX( ID_FINAL_ONLY, FluxMapForm::OnCommand )
	EVT_BUTTON( ID_EXPORT_TECPLOT, FluxMapForm::OnCommand )
	EVT_CHOICE( ID_COLOR_SCHEME, FluxMapForm::OnCommand )
END_EVENT_TABLE()

FluxMapForm::FluxMapForm( wxWindow *parent, Project &p )
	: wxPanel( parent ), m_prj(p), m_es(p)
{

	m_elementList = new ElementListBox( m_prj, this, ID_ELEMENT_LIST );


	wxFlexGridSizer *sizer2 = new wxFlexGridSizer( 4, wxSize(2,2) );
	sizer2->Add( new wxStaticText( this, wxID_ANY, "X bins:" ), 0, wxLEFT|wxRIGHT|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 5 );
	sizer2->Add( m_numXBins = new wxNumericCtrl( this, ID_XBINS, 30, wxNUMERIC_INTEGER ) );
	sizer2->Add( new wxStaticText( this, wxID_ANY, "Y bins:" ), 0, wxLEFT|wxRIGHT|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 5 );
	sizer2->Add( m_numYBins = new wxNumericCtrl( this, ID_YBINS, 30, wxNUMERIC_INTEGER ) );
	
	sizer2->Add( m_autoExtents = new wxCheckBox( this, ID_AUTO_EXTENT, "Autoscale") , 0, wxALIGN_CENTER_VERTICAL|wxALL, 2 );
	sizer2->AddStretchSpacer();
	sizer2->AddStretchSpacer();
	sizer2->AddStretchSpacer();

	sizer2->Add( new wxStaticText( this, wxID_ANY, "Min X:"), 0, wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 0 );
	sizer2->Add( m_minX = new wxNumericCtrl( this, ID_MIN_X ) );
	sizer2->Add( new wxStaticText( this, wxID_ANY, "Max X:"), 0, wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 0  );
	sizer2->Add( m_maxX = new wxNumericCtrl( this, ID_MAX_X ) );
	sizer2->Add( new wxStaticText( this, wxID_ANY, "Min Y:"), 0, wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 0  );
	sizer2->Add( m_minY = new wxNumericCtrl( this, ID_MIN_Y ) );
	sizer2->Add( new wxStaticText( this, wxID_ANY, "Max Y:"), 0, wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 0  );
	sizer2->Add( m_maxY = new wxNumericCtrl( this, ID_MAX_Y ) );

	
	sizer2->Add( new wxStaticText( this, wxID_ANY, "DNI:"), 0, wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 0  );
	sizer2->Add( m_dni = new wxNumericCtrl( this, ID_DNI, 1000.0 ) );
	sizer2->AddStretchSpacer();
	sizer2->Add( m_finalOnly = new wxCheckBox(this, ID_FINAL_ONLY, "Final only"), 0, wxALIGN_CENTER_VERTICAL|wxALL, 2 );
	
	sizer2->Add( new wxStaticText( this, wxID_ANY, "Contours:"), 0, wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 0  );
	sizer2->Add( m_contourLevels = new wxNumericCtrl( this, ID_CONTOUR_LEVELS, 15, wxNUMERIC_UNSIGNED ) );

	wxArrayString schemes;
	schemes.Add( "Jet" );
	schemes.Add( "Parula" );
	schemes.Add( "Grayscale" );

	sizer2->Add( new wxStaticText( this, wxID_ANY, "Scheme:" ), 0, wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 0 );
	sizer2->Add( m_colorScheme = new wxChoice( this, ID_COLOR_SCHEME, wxDefaultPosition, wxDefaultSize, schemes ), 0, wxALL|wxALIGN_CENTER_VERTICAL, 0 );
	m_colorScheme->SetSelection( 1 );

	wxBoxSizer *sizer_left = new wxBoxSizer( wxVERTICAL );
	sizer_left->Add( new wxStaticText( this, wxID_ANY, "Available elements (flat or cylindrical):" ), 0, wxALL, 2 );
	sizer_left->Add( m_elementList, 1, wxALL|wxEXPAND, 0 );
	sizer_left->Add( sizer2, 0, wxALL|wxEXPAND, 5 );
	sizer_left->Add( new wxButton( this, ID_EXPORT_TECPLOT, "Write Tecplot data files (.tec and .flx)...", wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT ) );


	m_plot = new wxPLPlotCtrl( this, wxID_ANY );
	//m_summary = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxTE_READONLY|wxTE_MULTILINE );
	m_summary = new wxStaticText( this, wxID_ANY, "Summary" );
	m_summary->SetFont( wxFont( 10, wxFONTFAMILY_MODERN, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL, false ) );
	m_summary->SetBackgroundColour( *wxWHITE );
	m_summary->SetForegroundColour( "Navy" );
	
	wxBoxSizer *sizer_right = new wxBoxSizer( wxVERTICAL );
	sizer_right->Add( m_plot, 1, wxALL|wxEXPAND, 0 );
	sizer_right->Add( m_summary, 0, wxALL|wxEXPAND, 0 );

	wxBoxSizer *sizer = new wxBoxSizer( wxHORIZONTAL );
	sizer->Add( sizer_left, 0, wxALL|wxEXPAND, 3 );
	sizer->Add( sizer_right, 2, wxALL|wxEXPAND, 0 );
	SetSizer( sizer );

	m_autoExtents->SetValue( true );
	m_minX->Enable( false );
	m_maxX->Enable( false );
	m_minY->Enable( false );
	m_maxY->Enable( false );
}

FluxMapForm::~FluxMapForm()
{
	/* nothing to do */
}

void FluxMapForm::UpdateList()
{
	m_elementList->Clear();
	for (size_t i = 0; i < m_prj.StageList.size(); i++)
	{
		for (size_t j = 0; j < m_prj.StageList[i]->ElementList.size(); j++)
		{
			if (m_prj.StageList[i]->ElementList[j]->SurfaceIndex == 't' || m_prj.StageList[i]->ElementList[j]->SurfaceIndex == 'f')
				m_elementList->Add(i, j);
			if (m_prj.StageList[i]->ElementList[j]->SurfaceIndex == 's' && m_prj.StageList[i]->ElementList[j]->ApertureIndex == 'l')
				m_elementList->Add(i, j);
		}
	}
	
	m_elementList->Invalidate(); // update scrollbar	
	m_elementList->Refresh();
}

void FluxMapForm::GetCurrentSelection( int *stageIdx, int *elementIdx, char *surfidx )
{
	*stageIdx = *elementIdx = -1;

	Element *elm = 0;
	for( size_t i=0;i<m_prj.StageList.size();i++ )
	{
		for (size_t j=0;j<m_prj.StageList[i]->ElementList.size();j++)
		{
			if ( m_elementList->IsSelected( i,j ) )
			{
				elm = m_prj.StageList[i]->ElementList[j];
				*stageIdx = i;
				*elementIdx = j;
				if ( surfidx ) *surfidx = elm->SurfaceIndex;

				break;
			}
		}
	}
}

void FluxMapForm::UpdatePlot()
{
	m_plot->SetTitle( "" );
	m_plot->DeleteAllPlots();
	m_plot->SetSideWidget( NULL );

	
	int stageIdx=-1, elementIdx=-1;
	char surfIdx = '0';
	GetCurrentSelection( &stageIdx, &elementIdx, &surfIdx );
	if ( stageIdx < 0 || elementIdx < 0 )
		return;


	double minx = m_minX->Value();
	double maxx = m_maxX->Value();
	double miny = m_minY->Value();
	double maxy = m_maxY->Value();

	RayData &rd = m_prj.Results;
	
	int nbinsx = m_numXBins->AsInteger();
	int nbinsy = m_numYBins->AsInteger();

	if (nbinsx <= 1 || nbinsy <= 1
		|| maxx < minx || maxy < miny )
	{
		m_summary->SetLabel("Error in flux map plot parameters.");
		return;
	}

	wxString xlabel, ylabel, title;
	xlabel = "X";
	ylabel = "Y";
	title = "Flat Surface Flux Plot";

	if ( surfIdx =='t' || surfIdx == 's')
	{
		xlabel = "Circumference";
		ylabel = "Length";
		title = "Cylindrical Flux Plot";
	}

	if (!m_es.Compute(
				stageIdx,
				elementIdx,
				nbinsx,
				nbinsy,
				m_autoExtents->GetValue(),
				m_finalOnly->GetValue(),
				m_dni->Value(),
				minx, miny, maxx, maxy ))
	{
		m_summary->SetLabel("Error in flux map plot parameters.");
		return;
	}

	double xextend = 0;
	double yextend = 0;

	if ( m_autoExtents->GetValue() )
	{
		// increase the view area slightly
		// so that all rays are plotted when autoscaling
		// but don't change min/max for the binning routines
		xextend = (maxx-minx)*0.01;
		yextend = (maxy-miny)*0.01;

		m_minX->SetValue( minx - xextend );
		m_maxX->SetValue( maxx + xextend );
		m_minY->SetValue( miny - yextend );
		m_maxY->SetValue( maxy + yextend );

	}

	wxString details;

	details += wxString::Format("Sun ray count: %d, over box of dimensions: %lg x %lg\n",
			    rd.SunRayCount, rd.SunXMax - rd.SunXMin, rd.SunYMax - rd.SunYMin );
	details += wxString::Format("Power per ray: %lg\n", m_es.PowerPerRay);
	details += wxString::Format("Peak flux: %lg\n", m_es.PeakFlux);
	details += wxString::Format("Peak flux uncertainty: +/- %lg %%\n", m_es.PeakFluxUncertainty);
	details += wxString::Format("Min flux: %lg\n", m_es.MinFlux);
	details += wxString::Format("Sigma flux: %lg\n",  m_es.SigmaFlux);
	details += wxString::Format("Avg. flux: %lg\n", m_es.AveFlux);
	details += wxString::Format("Avg. flux uncertainty: +/- %lg %%\n", m_es.AveFluxUncertainty);
	details += wxString::Format("Uniformity: %lg\n", m_es.Uniformity);
	details += wxString::Format("Power of plotted rays: %lg\n", m_es.NumberOfRays*m_es.PowerPerRay);
	details += wxString::Format("Centroid: ( %lg, %lg, %lg )\n", m_es.Centroid[0], m_es.Centroid[1], m_es.Centroid[2]);

	m_summary->SetLabel(details);
	Layout();


	wxMatrix<double> xx, yy, zz;
	wxPLContourPlot::MeshGrid( minx, maxx, nbinsx, miny, maxy, nbinsy, xx, yy );

	double zmin=0, zmax=0;
	zz.Resize( xx.Rows(), xx.Cols() );
	for( int i=0;i<nbinsy;i++ )
	{
		for( int j=0;j<nbinsx;j++ )
		{
			double z = m_es.zScale * m_es.fluxGrid(j,i);
			zz(i,j) = z;
			if ( z < zmin ) zmin = z;
			if ( z > zmax ) zmax = z;
		}
	}

	int nlev = m_contourLevels->AsInteger();
	if ( nlev < 2 ) nlev = 2;
	if ( nlev > 50 ) nlev = 50;

	if ( nlev != m_contourLevels->AsInteger() )
		m_contourLevels->SetValue( nlev );

	wxPLColourMap *cmap = 0;
	switch( m_colorScheme->GetSelection() )
	{
	case 1: cmap = new wxPLParulaColourMap( zmin, zmax ); break;
	case 2: cmap = new wxPLGrayscaleColourMap( zmin, zmax ); break;
	//case 3: cmap = new wxPLCoarseRainbowColourMap( zmin, zmax ); break;
	default:
		cmap = new wxPLJetColourMap( zmin, zmax );
	}

	wxPLContourPlot *cntr = new wxPLContourPlot( xx, yy, zz, true, "", nlev, cmap );
	
	m_plot->SetSideWidget( cmap );
	m_plot->AddPlot( cntr );
	m_plot->SetTitle( "Flux Intensity" );
	m_plot->Refresh();
}

void FluxMapForm::OnCommand( wxCommandEvent &evt )
{
	switch( evt.GetId() )
	{
	case ID_AUTO_EXTENT:
		m_minX->Enable( !m_autoExtents->GetValue() );
		m_maxX->Enable( !m_autoExtents->GetValue() );
		m_minY->Enable( !m_autoExtents->GetValue() );
		m_maxY->Enable( !m_autoExtents->GetValue() );

	case ID_CONTOUR_LEVELS:
	case ID_DNI:
	case ID_XBINS:
	case ID_YBINS:
	case ID_MIN_X:
	case ID_MIN_Y:
	case ID_MAX_X:
	case ID_MAX_Y:
	case ID_FINAL_ONLY:
	case ID_COLOR_SCHEME:
		UpdatePlot();
		break;

	case ID_ELEMENT_LIST:
		{
			int row = evt.GetSelection();
			unsigned int s, e;
			if ( m_elementList->GetStageElementIndices( row, &s, &e ) )
			{
				bool sel = m_elementList->IsSelected( s, e );
				m_elementList->ClearSelections();
				m_elementList->Select( s, e, !sel );
				m_elementList->Refresh();
				UpdatePlot();
			}
		}
		break;

	case ID_EXPORT_TECPLOT:
		WriteTecFlx();
		break;
	}
}


void FluxMapForm::WriteTecFlx()
{
	int stageIdx, elementIdx;
	GetCurrentSelection( &stageIdx, &elementIdx );

	if ( stageIdx < 0 || elementIdx < 0 )
	{
		wxMessageBox( "An element must be selected for export data." );
		return;
	}

	wxString prjfile = MainWindow::Instance().GetFileName();

	wxFileDialog dlg( this, "Export Tecplot Data Files", 
		wxPathOnly( prjfile ), prjfile, 
		"All Files (*.*)|*.*", wxFD_OPEN );
				

	if ( dlg.ShowModal() != wxID_OK )
		return;

	wxString file = dlg.GetPath();

	double dni = m_es.DNI;
	if (dni == 0) dni = 1;

	int nbinsx = m_numXBins->AsInteger();
	int nbinsy = m_numYBins->AsInteger();

	bool create_flx = true;
	if ( wxFileExists( file+".flx")
		&& wxNO == wxMessageBox( "The FLX file already exists, overwrite?", "Query", wxYES_NO ))
		create_flx = false;

	if (create_flx)
	{
		FILE *fp = fopen( wxString(file+".flx").c_str(), "w");
		if (fp)
		{
			fprintf(fp, "Stage %d Element %d Flux\n", stageIdx+1, elementIdx+1);
			fprintf(fp, "# Bins in X, # Bins in Y\n");
			fprintf(fp, "                     %d                     %d\n", nbinsx, nbinsy);
			fprintf(fp, "                    ");
			for (int i=0;i<m_es.xValues.size();i++)
				fprintf(fp, "   %.10le", m_es.xValues[i]);

			fprintf(fp, "\n");

			for (int j=0;j<m_es.yValues.size();j++)
			{
				fprintf(fp,"%.10le     ",m_es.yValues[ m_es.yValues.size() - j - 1]);
				for (int i=0;i<m_es.xValues.size();i++)
					fprintf(fp, "   %.10le", m_es.fluxGrid.at(i,j)*m_es.PowerPerRay/(m_es.binszx*m_es.binszy)/dni);

				fprintf(fp, "\n");
			}

			fclose(fp);
		}
	}

	bool create_tec = true;
	if ( wxFileExists(file+".tec")
		&& wxNO == wxMessageBox( "The TEC file already exists, overwrite it?", "Query", wxYES_NO ) )
		create_tec = false;


	if (create_tec)
	{
		FILE *fp = fopen(wxString(file+".tec").c_str(), "w");
		if (fp)
		{
			wxString title = wxString::Format("Stage %d Element %d Flux Map", stageIdx+1, elementIdx+1);

			fprintf(fp, "TITLE=\"%s\"\n", (const char*)title.c_str());
			fprintf(fp, "VARIABLES= \"X\", \"Y\", \"Z\", \"Flux, kW/m<sup>2\"\n");
			fprintf(fp, "ZONE T=\"%s\"\n", (const char*)title.c_str());
			fprintf(fp, "I=%d, J=%d, K=1\n", (int)m_es.xValues.size(), (int)m_es.yValues.size());
			fprintf(fp, "DATAPACKING=POINT\n");
			fprintf(fp, "DT=(SINGLE, SINGLE, SINGLE, SINGLE)\n");

			for (int i=0;i<m_es.xValues.size();i++)
			{
				for (int j=0;j<m_es.yValues.size();j++)
				{
					fprintf(fp, "%.10le %.10le %.10le %.10le\n",
							m_es.xValues[i],
							m_es.yValues[m_es. yValues.size() - j - 1],
							0.0,
							m_es.fluxGrid.at(i,j)*m_es.PowerPerRay/(m_es.binszx*m_es.binszy)/dni);
				}
			}

			fclose(fp);
		}
	}
}