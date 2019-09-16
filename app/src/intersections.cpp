
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


#include <vector>
#include <algorithm>
#include <unordered_map>

#include <wx/panel.h>
#include <wx/dcbuffer.h>
#include <wx/button.h>
#include <wx/checkbox.h>
#include <wx/checklst.h>
#include <wx/slider.h>
#include <wx/choice.h>
#include <wx/sizer.h>
#include <wx/dcclient.h>
#include <wx/stattext.h>
#include <wx/statline.h>
#include <wx/textdlg.h>
#include <wx/msgdlg.h>
#include <wx/busyinfo.h>
#include <wx/vlbox.h>
#include <wx/renderer.h>

#include <wex/numeric.h>

#include <wex/plot/plaxis.h>
#include <wex/plot/plplot.h>
#include <wex/utils.h>

#include "intersections.h"
#include "elementlist.h"

IntersectionViewer::IntersectionViewer( wxWindow *parent, Project &prj )
	: wxGLEasyCanvas( parent, wxID_ANY ), m_prj(prj)
{
	m_showAxes = true;
	m_showTicks = true;
	m_vertexListValid = false;

	m_coordSys = RayData::COORD_GLOBAL;
	m_pointColors = 0;
	m_finalOnly = false;
	m_includeMissedRays = true;
	
	m_centroid[0] = m_centroid[1] = m_centroid[2] = 0.0;
	m_nplotted = 0;

	m_colorList.push_back( *wxBLACK );
	m_colorList.push_back( *wxRED );
	m_colorList.push_back( *wxBLUE );
	m_colorList.push_back( "Forest Green" );
	m_colorList.push_back( "Magenta" );
	m_colorList.push_back( "Orange" );
	m_colorList.push_back( "Sea Green" );
	m_colorList.push_back( "Grey" );
	m_colorList.push_back( "Salmon" );

}

IntersectionViewer::~IntersectionViewer()
{
	if ( m_vertexListValid )
		glDeleteLists( m_vertexListId, 1 );
}

void IntersectionViewer::SetupAxes( bool show, bool ticks )
{
	m_showAxes = show;
	m_showTicks = ticks;
}


void IntersectionViewer::Invalidate( ElementListBox *lb )
{
	bool fit = ( m_nplotted == 0 );
	
	RebuildGeometry( lb );
	
	if ( m_nplotted > 0 && fit )
		FitView();

	Refresh(false);
}

void IntersectionViewer::FitView()
{
	float dX = std::max(fabs( m_max.x ), fabs( m_min.x ));
	float dY = std::max(fabs( m_max.y ), fabs( m_min.y )); 
	float dZ = std::max(fabs( m_max.z ), fabs( m_min.z ));

	wxSize client( GetClientSize() );
	float cdim = std::min( client.x, client.y );

	float view = GetViewSize(); // default 100
	float geom = std::max( std::max( dX, dY ), dZ ); // say, 30

	float zoom = geom > 0  ? view/geom * 0.5f : 1.0;

	zoom *= ( cdim/view );

	if ( zoom < 0.001 ) zoom = 0.001;
	if ( zoom > 1000 ) zoom = 1000;
	SetZoom( zoom );
}

struct point_info
{
	size_t stage;
	size_t element;
	size_t index;
	
	static bool compare( const point_info &a, const point_info &b )
	{
		return a.stage <= b.stage && a.element < b.element;
	}
};

struct ray_segment {
	wxGLPoint3D point;
	wxColour color;
};

void IntersectionViewer::RebuildGeometry( ElementListBox *lb )
{
	wxBusyInfo info("Building 3D intersection point view...", this);

	SetCurrent( m_glContext );
	
	double Pos[3], Cos[3];
	int Elm, Stg, Ray;	
	
	RayData &R = m_prj.Results;
	size_t len = R.Length;
	if ( len == 0 ) return;
			
	// create new GL call list
	if ( !m_vertexListValid )
	{
		m_vertexListId = glGenLists( 1 );
		m_vertexListValid = true;
	}
	
	// recreate the vertex call list
	glNewList( m_vertexListId, GL_COMPILE );	
	
	glColor3f( 0, 0, 0 );
	glBegin( GL_POINTS );

	m_nplotted = 0;
	m_centroid[0] = m_centroid[1] = m_centroid[2] = 0.0;

	// store the last point for the ray line plots
	std::vector< std::vector<ray_segment> > rays;
	if ( m_rayNumbers.size() > 0 )
		rays.resize( m_rayNumbers.size(), std::vector<ray_segment>() ); 

	int last_stage = -1;
	int last_elem = -1;
	wxColour cur_color( *wxBLACK );
	wxColour last_color( *wxWHITE );
	m_min = wxGLPoint3D( 1e31f, 1e31f, 1e31f );
	m_max = wxGLPoint3D( -1e31f, -1e31f, -1e31f );
	
	// make list of all points to draw
	for( size_t i=0;i<len;i++ )
	{
		int istage = R.StageMap[i]-1;
		int ielem = abs(R.ElementMap[i])-1;
		bool absorbed = ( R.ElementMap[i] < 0 );

		bool selected = (lb && lb->IsSelected( istage, ielem ));
		
		if ( R.Transform( m_prj, m_coordSys, i, Pos, Cos, Elm, Stg, Ray ) )
		{
			if ( selected )
			{
				if ( m_pointColors == 1 )
					cur_color = m_colorList[ ((size_t)istage) % m_colorList.size() ];
				else if ( m_pointColors == 2 )
					cur_color = m_colorList[ ((size_t)ielem) % m_colorList.size() ];
				else
					cur_color = *wxBLACK;

				if ( cur_color != last_color )
					Color( cur_color );

				glVertex3f( (float)Pos[0], (float)Pos[1], (float)Pos[2] );

				if ( ((float)Pos[0]) < m_min.x ) m_min.x = (float)Pos[0];
				if ( ((float)Pos[0]) > m_max.x ) m_max.x = (float)Pos[0];
		
				if ( ((float)Pos[1]) < m_min.y ) m_min.y = (float)Pos[1];
				if ( ((float)Pos[1]) > m_max.y ) m_max.y = (float)Pos[1];
			
				if ( ((float)Pos[2]) < m_min.z ) m_min.z = (float)Pos[2];
				if ( ((float)Pos[2]) > m_max.z ) m_max.z = (float)Pos[2];
			
				m_centroid[0] += Pos[0];
				m_centroid[1] += Pos[1];
				m_centroid[2] += Pos[2];
				m_nplotted++;
				
				last_stage = istage;
				last_elem = ielem;
				last_color = cur_color;
			}

			if ( m_rayNumbers.size() > 0 
				&& ( Elm != 0
					|| (Elm == 0 && m_includeMissedRays) ) )
			{
				for( size_t k=0;k<m_rayNumbers.size();k++ )
				{
					if ( m_rayNumbers[k] == Ray ) 
					{
						if ( rays[k].size() == 0 )
						{
							ray_segment rs;
							rs.point.x = (float)(Pos[0] - Cos[0]*20.0);
							rs.point.y = (float)(Pos[1] - Cos[1]*20.0);
							rs.point.z = (float)(Pos[2] - Cos[2]*20.0);
							rs.color = Elm==0 ? *wxRED : wxColour(255,192,0);
							rays[k].push_back( rs );
						}

						wxGLPoint3D P( Pos[0], Pos[1], Pos[2] );
						if ( 0 == Elm )
						{
							P.x += (float)(Cos[0]*5.0);
							P.y += (float)(Cos[1]*5.0);
							P.z += (float)(Cos[2]*5.0);
						}

						ray_segment R;
						R.point = P;
						R.color = Elm==0 ? *wxRED : wxColour(255,192,0);
						rays[k].push_back( R );
					}
				}
			}
		}
	}

	glEnd(); // POINTS


	if ( rays.size() > 0 )
	{
		for( size_t k=0;k<rays.size();k++ )
		{
			const std::vector<ray_segment> &R = rays[k];
			if ( R.size() > 0 )
			{
				for( size_t j=0;j<R.size()-1;j++ )
				{
					Color( R[j+1].color );
					glBegin( GL_LINES );
						glVertex3f( R[j].point.x, R[j].point.y, R[j].point.z );
						glVertex3f( R[j+1].point.x, R[j+1].point.y, R[j+1].point.z );
					glEnd();
				}
			}
		}
	}

	
	glEndList();

	if ( m_nplotted > 0 )
	{
		m_centroid[0] /= m_nplotted;
		m_centroid[1] /= m_nplotted;
		m_centroid[2] /= m_nplotted;
	}
}

void IntersectionViewer::OnRender()
{
	if ( m_vertexListValid )
	{
		glEnable( GL_DEPTH_TEST );
		glCallList( m_vertexListId );

		if ( m_showAxes )
			Axes( m_min, m_max, m_showTicks ? 9.0f : 0.0f, true );

	}
	else
		Text( 0, 0, "No intersection data to display.", *wxWHITE, *wxBLACK_BRUSH );
}

enum { __idFirst = wxID_HIGHEST+31,
	ID_COORD_SYS, 
	ID_STAGE_SEL_ALL, ID_STAGE_SEL_NONE, ID_STAGE_SEL_SOME,
	ID_ELEMENT_SEL_ALL, ID_ELEMENT_SEL_NONE, ID_ELEMENT_SEL_SOME,
	ID_STAGE_LIST, ID_ELEMENT_LIST,
	ID_FINAL_ONLY, ID_POINT_COLOR, ID_DNI, ID_PLOT_RAYS, ID_RAY_NUMBERS, ID_INCL_MISSED,
	ID_SHOW_AXES, ID_SHOW_TICKS, ID_AXIS_TEXT_SIZE, ID_AXIS_ZOOM_RATE,
	ID_SCALE_X, ID_SCALE_Y, ID_SCALE_Z
};

BEGIN_EVENT_TABLE( IntersectionForm, wxWindow )
	EVT_CHOICE( ID_COORD_SYS, IntersectionForm::OnCommand )

	EVT_BUTTON( ID_STAGE_SEL_ALL, IntersectionForm::OnCommand )
	EVT_BUTTON( ID_STAGE_SEL_NONE, IntersectionForm::OnCommand )
	EVT_BUTTON( ID_STAGE_SEL_SOME, IntersectionForm::OnCommand )

	EVT_BUTTON( ID_ELEMENT_SEL_ALL, IntersectionForm::OnCommand )
	EVT_BUTTON( ID_ELEMENT_SEL_NONE, IntersectionForm::OnCommand )
	EVT_BUTTON( ID_ELEMENT_SEL_SOME, IntersectionForm::OnCommand )

	EVT_CHECKLISTBOX( ID_STAGE_LIST, IntersectionForm::OnCommand )
	EVT_LISTBOX( ID_ELEMENT_LIST, IntersectionForm::OnCommand )

	EVT_CHECKBOX( ID_FINAL_ONLY, IntersectionForm::OnCommand )
	EVT_CHOICE( ID_POINT_COLOR,  IntersectionForm::OnCommand )
	EVT_NUMERIC( ID_DNI, IntersectionForm::OnCommand )

	EVT_CHECKBOX( ID_PLOT_RAYS, IntersectionForm::OnCommand )
	EVT_TEXT_ENTER( ID_RAY_NUMBERS, IntersectionForm::OnCommand )
	EVT_CHECKBOX( ID_INCL_MISSED, IntersectionForm::OnCommand )
	
	EVT_CHECKBOX( ID_SHOW_AXES, IntersectionForm::OnCommand )
	EVT_CHECKBOX( ID_SHOW_TICKS, IntersectionForm::OnCommand )
	EVT_SLIDER( ID_AXIS_TEXT_SIZE, IntersectionForm::OnCommand )
	EVT_NUMERIC( ID_SCALE_X, IntersectionForm::OnCommand )
	EVT_NUMERIC( ID_SCALE_Y, IntersectionForm::OnCommand )
	EVT_NUMERIC( ID_SCALE_Z, IntersectionForm::OnCommand )
    EVT_TEXT_ENTER(ID_AXIS_ZOOM_RATE, IntersectionForm::OnCommand )
END_EVENT_TABLE()

IntersectionForm::IntersectionForm( wxWindow *parent, Project &prj )
	: wxWindow( parent, wxID_ANY ), m_prj(prj)
{
	wxPanel *panel = new wxPanel(this);

	wxString coordopts[] = { "Global", "Stage", "Element" };
	m_coordSys = new wxChoice( panel, ID_COORD_SYS, wxDefaultPosition, wxDefaultSize, 3, coordopts );

	wxBoxSizer *sizer1 = new wxBoxSizer( wxHORIZONTAL );
	sizer1->Add( new wxStaticText( panel, wxID_ANY, "Coordinate system:" ), 0, wxALIGN_CENTER_VERTICAL|wxLEFT|wxRIGHT, 5 );
	sizer1->Add( m_coordSys, 1, wxALL|wxALIGN_CENTER_VERTICAL, 5 );

	wxBoxSizer *sizer2 = new wxBoxSizer( wxHORIZONTAL );
	sizer2->Add( new wxStaticText( panel, wxID_ANY, "Stages:" ), 0, wxALIGN_CENTER_VERTICAL|wxLEFT|wxRIGHT, 1 );
	sizer2->AddStretchSpacer();
	sizer2->AddSpacer( 20 );
	sizer2->Add( new wxButton( panel, ID_STAGE_SEL_ALL, "All", wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT ), 0, wxALIGN_CENTER_VERTICAL|wxALL, 1 );
	sizer2->Add( new wxButton( panel, ID_STAGE_SEL_NONE, "None", wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT ), 0, wxALIGN_CENTER_VERTICAL|wxALL, 1 );
	sizer2->Add( new wxButton( panel, ID_STAGE_SEL_SOME, "Some", wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT ), 0, wxALIGN_CENTER_VERTICAL|wxALL, 1 );

	m_stageList = new wxCheckListBox( panel, ID_STAGE_LIST );

	
	wxBoxSizer *sizer3 = new wxBoxSizer( wxHORIZONTAL );
	sizer3->Add( new wxStaticText( panel, wxID_ANY, "Elements:" ), 0, wxALIGN_CENTER_VERTICAL|wxLEFT|wxRIGHT, 1 );
	sizer3->AddStretchSpacer();
	sizer3->AddSpacer( 20 );
	sizer3->Add( new wxButton( panel, ID_ELEMENT_SEL_ALL, "All", wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT ), 0, wxALIGN_CENTER_VERTICAL|wxALL, 1 );
	sizer3->Add( new wxButton( panel, ID_ELEMENT_SEL_NONE, "None", wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT ), 0, wxALIGN_CENTER_VERTICAL|wxALL, 1 );
	sizer3->Add( new wxButton( panel, ID_ELEMENT_SEL_SOME, "Some", wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT ), 0, wxALIGN_CENTER_VERTICAL|wxALL, 1 );

	m_elementList = new ElementListBox( m_prj, panel, ID_ELEMENT_LIST );

	m_finalOnly = new wxCheckBox( panel, ID_FINAL_ONLY, "Show final intersections only" );

	m_dni = new wxNumericCtrl( panel, ID_DNI, 1000.0 );
	wxBoxSizer *sizer4 = new wxBoxSizer( wxHORIZONTAL );
	sizer4->Add( new wxStaticText( panel, wxID_ANY, "DNI for calculations:" ), 0, wxALIGN_CENTER_VERTICAL|wxLEFT|wxRIGHT, 5 );
	sizer4->AddStretchSpacer();
	sizer4->Add( m_dni, 0, wxALL|wxALIGN_CENTER_VERTICAL, 5 );

	m_plotRays = new wxCheckBox( panel, ID_PLOT_RAYS, "Plot paths for rays: " );
	m_rayNumbers = new wxTextCtrl( panel, ID_RAY_NUMBERS, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxTE_PROCESS_ENTER );
	wxBoxSizer *sizer5 = new wxBoxSizer( wxHORIZONTAL );
	sizer5->Add( m_plotRays, 0, wxALIGN_CENTER_VERTICAL|wxLEFT|wxRIGHT, 5 );
	sizer5->AddStretchSpacer();
	sizer5->Add( m_rayNumbers, 0, wxLEFT|wxRIGHT|wxALIGN_CENTER_VERTICAL, 5 );

	m_inclMissedRays = new wxCheckBox( panel, ID_INCL_MISSED, "Include missed rays" );
	m_inclMissedRays->SetValue( true );

	m_summary = new wxTextCtrl( panel, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxTE_MULTILINE|wxTE_READONLY );

	
	wxString coloropts[] = { "uniform", "by stage", "by element" };
	m_pointColors = new wxChoice( panel, ID_POINT_COLOR, wxDefaultPosition, wxDefaultSize, 3, coloropts );
	wxBoxSizer *sizer6 = new wxBoxSizer( wxHORIZONTAL );
	sizer6->Add( new wxStaticText( panel, wxID_ANY, "Point colors:" ), 0, wxRIGHT|wxALIGN_CENTER_VERTICAL, 5 );
	sizer6->Add( m_pointColors, 1, wxLEFT|wxALIGN_CENTER_VERTICAL, 5 );
	
	m_showAxes = new wxCheckBox( panel, ID_SHOW_AXES, "Show axes" );
	m_showAxes->SetValue( true );
	m_showTicks = new wxCheckBox( panel, ID_SHOW_TICKS, "Show tick marks" );
	m_showTicks->SetValue( true );

	wxBoxSizer *sizer7 = new wxBoxSizer( wxHORIZONTAL );
	sizer7->Add( m_showAxes, 1, wxALL|wxEXPAND, 5 );
	sizer7->Add( m_showTicks, 1, wxALL|wxEXPAND, 5 );
	
	m_axisTextSize = new wxSlider( panel, ID_AXIS_TEXT_SIZE, wxNORMAL_FONT->GetPointSize(), 6, 20 );
	wxBoxSizer *sizer8 = new wxBoxSizer( wxHORIZONTAL );
	sizer8->Add( new wxStaticText( panel, wxID_ANY, "Text size:" ), 0, wxLEFT|wxRIGHT|wxALIGN_CENTER_VERTICAL, 5 );
	sizer8->Add( m_axisTextSize, 1, wxLEFT|wxRIGHT|wxALIGN_CENTER_VERTICAL, 5 );

	wxFlexGridSizer *sizer9 = new wxFlexGridSizer( 2, wxSize(5,5) );
	sizer9->AddGrowableCol( 0 );
	sizer9->Add( new wxStaticText( panel, wxID_ANY, "X axis scale factor:" ), wxALIGN_RIGHT|wxLEFT|wxRIGHT|wxALIGN_CENTER_VERTICAL, 5 );
	sizer9->Add( m_scaleX = new wxNumericCtrl( panel, ID_SCALE_X, 1.0 ), wxALIGN_CENTER_VERTICAL|wxLEFT|wxRIGHT, 5 );
	sizer9->Add( new wxStaticText( panel, wxID_ANY, "Y axis scale factor:" ), wxALIGN_RIGHT|wxLEFT|wxRIGHT|wxALIGN_CENTER_VERTICAL, 5 );
	sizer9->Add( m_scaleY = new wxNumericCtrl( panel, ID_SCALE_Y, 1.0 ), wxALIGN_CENTER_VERTICAL|wxLEFT|wxRIGHT, 5 );
	sizer9->Add( new wxStaticText( panel, wxID_ANY, "Z axis scale factor:" ), wxALIGN_RIGHT|wxLEFT|wxRIGHT|wxALIGN_CENTER_VERTICAL, 5 );
	sizer9->Add( m_scaleZ = new wxNumericCtrl( panel, ID_SCALE_Z, 1.0 ), wxALIGN_CENTER_VERTICAL|wxLEFT|wxRIGHT, 5 );
    sizer9->Add( new wxStaticText( panel, wxID_ANY, "Plot zoom rate (lower faster):"), wxALIGN_RIGHT|wxLEFT|wxRIGHT|wxALIGN_CENTER_VERTICAL, 5);
    sizer9->Add( m_scaleZoom = new wxNumericCtrl( panel, ID_AXIS_ZOOM_RATE, 10.), wxALIGN_CENTER_VERTICAL|wxLEFT|wxRIGHT, 5 );

	wxBoxSizer *panel_sizer = new wxBoxSizer( wxVERTICAL );
	panel_sizer->Add( sizer1, 0, wxALL|wxEXPAND, 0 );
	
	panel_sizer->Add( sizer2, 0, wxLEFT|wxRIGHT|wxEXPAND, 5 );
	panel_sizer->Add( m_stageList, 1, wxALL|wxEXPAND, 5 );

	panel_sizer->Add( sizer3, 0, wxLEFT|wxRIGHT|wxEXPAND, 5 );
	panel_sizer->Add( m_elementList, 1, wxALL|wxEXPAND, 5 );

	panel_sizer->Add( m_finalOnly, 0, wxALL|wxEXPAND, 5 );
		
	panel_sizer->Add( sizer4, 0, wxALL|wxEXPAND, 0 );
	panel_sizer->Add( sizer5, 0, wxALL|wxEXPAND, 0 );
	panel_sizer->Add( m_inclMissedRays, 0, wxALL|wxEXPAND, 5 );

	panel_sizer->Add( m_summary, 1, wxALL|wxEXPAND, 5 );

	panel_sizer->Add( new wxStaticLine( panel, wxID_ANY ), 0, wxALL|wxEXPAND, 1 );
	panel_sizer->Add( sizer6, 0, wxALL|wxEXPAND, 5 );
	panel_sizer->Add( sizer7, 0, wxALL|wxEXPAND, 0 );
	panel_sizer->Add( sizer9, 0, wxALL|wxEXPAND, 5 );
	panel_sizer->Add( sizer8, 0, wxALL|wxEXPAND, 0 );


	panel_sizer->AddSpacer( 50 );

	panel->SetSizer( panel_sizer );
	
	
	m_3d = new IntersectionViewer( this, prj );
	m_3d->SetFont( *wxNORMAL_FONT );

	wxBoxSizer *sizer_main = new wxBoxSizer( wxHORIZONTAL );
	sizer_main->Add( panel, 0, wxALL|wxEXPAND, 0 );
	sizer_main->Add( m_3d, 1, wxALL|wxEXPAND, 0 );
	SetSizer( sizer_main );

	m_coordSys->SetSelection( 0 ); // global
	m_pointColors->SetSelection( 2 ); // by element
}

IntersectionForm::~IntersectionForm()
{
}

void IntersectionForm::UpdateView()
{
	PopulateStages();
	PopulateElements();
	UpdatePlot();
}

void IntersectionForm::OnCommand( wxCommandEvent &evt )
{
	wxString text;

	switch( evt.GetId() )
	{
	case ID_COORD_SYS:
	case ID_FINAL_ONLY:
	case ID_PLOT_RAYS:
	case ID_INCL_MISSED:
	case ID_RAY_NUMBERS:
	case ID_POINT_COLOR:		
		UpdatePlot();
		break;

	case ID_STAGE_SEL_ALL:
		for( size_t i=0;i<m_stageList->GetCount();i++ )
			m_stageList->Check( i, true );

		PopulateElements();
		UpdatePlot();
		break;

	case ID_STAGE_SEL_NONE:
		for( size_t i=0;i<m_stageList->GetCount();i++ )
			m_stageList->Check( i, false );
		
		m_elementList->ClearSelections();

		PopulateElements();
		UpdatePlot();
		break;

	case ID_STAGE_SEL_SOME:
	{
		text = wxGetTextFromUser( "Select a range of stages (eg. 1, 4-6):", "Query", wxEmptyString, this );
		if ( !text.IsEmpty() )
		{
			std::vector<int> list = ::wxCommaDashListToIndices( text );			
			for( size_t i=0;i<m_stageList->GetCount();i++ )
			{
				bool sel = std::find( list.begin(), list.end(), (int)i+1 ) != list.end();
				m_stageList->Check( i, sel );
			}
		}
		CullElementSelections();

		PopulateElements();
		UpdatePlot();
	}
		break;
		
	case ID_ELEMENT_SEL_NONE:
		m_elementList->ClearSelections();

		UpdatePlot();
		break;

	case ID_ELEMENT_SEL_ALL:
		m_elementList->SelectAll();
		UpdatePlot();
		break;
		
	case ID_ELEMENT_SEL_SOME:
	{
		if ( m_elementList->Count() == 0 ) return;

		wxArrayString stages;
		for( size_t i=0;i<m_stageList->GetCount();i++ )
			if ( m_stageList->IsChecked( i ) )
				stages.Add( m_stageList->GetString(i) );

		text = wxGetTextFromUser( "Select a range of element numbers in '" +
			wxJoin( stages, ',' ) + "' (eg. 1, 4-6):", "Query", wxEmptyString, this );

		if ( !text.IsEmpty() )
		{
			m_elementList->ClearSelections();
			std::vector<int> list = ::wxCommaDashListToIndices( text );				
			for( size_t i=0;i<m_elementList->Count();i++ )
			{
				unsigned int s, e;
				if ( m_elementList->GetStageElementIndices( i, &s, &e )
					&& std::find( list.begin(), list.end(), e+1 ) != list.end() )
				{
					m_elementList->Select( s, e );
				}
			}

			m_elementList->Refresh();
			
		}

		UpdatePlot();
	}
		break;
		
		
	case ID_STAGE_LIST:
	{
		CullElementSelections();
		PopulateElements();
		UpdatePlot();
	}
		break;

	case ID_ELEMENT_LIST:
	{
		size_t row = (size_t)evt.GetSelection();
		unsigned int stage_num, element_num;
		if ( m_elementList->GetStageElementIndices( row, &stage_num, &element_num ) )
		{
			if ( m_elementList->IsSelected( stage_num, element_num ) )
				m_elementList->Select( stage_num, element_num, false );
			else
				m_elementList->Select( stage_num, element_num, true );
		}
		m_elementList->Refresh();
		UpdatePlot();
	}
		break;

	case ID_DNI:
		UpdateDetails();
		break;

	case ID_SHOW_AXES:
	case ID_SHOW_TICKS:
		m_3d->SetupAxes( m_showAxes->GetValue(), m_showTicks->GetValue() );
		m_3d->Refresh( false );
		break;

	case ID_AXIS_TEXT_SIZE:
		{
			wxFont font(*wxNORMAL_FONT);
			font.SetPointSize( m_axisTextSize->GetValue() );
			m_3d->SetFont( font );
			m_3d->Refresh( false );
		}
		break;

	case ID_SCALE_X:
		m_3d->SetScaleX( m_scaleX->Value() );
		m_3d->Refresh( false );
		break;

	case ID_SCALE_Y:
		m_3d->SetScaleY( m_scaleY->Value() );
		m_3d->Refresh( false );
		break;

	case ID_SCALE_Z:
		m_3d->SetScaleZ( m_scaleZ->Value() );
		m_3d->Refresh( false );
		break;
    case ID_AXIS_ZOOM_RATE:
        m_3d->SetZoomRate( m_scaleZoom->Value() );
        m_3d->Refresh( false );
        break;
	}

}

void IntersectionForm::UpdateDetails()
{
	RayData &rd = m_prj.Results;  // untranslated stage coordinates

	double DNI = m_dni->Value();
	double PowerPerRay = DNI
		* (rd.SunXMax - rd.SunXMin)
		* (rd.SunYMax - rd.SunYMin)
		/ rd.SunRayCount;

	wxString details;
	details += wxString::Format("Sun ray count: %d\n", rd.SunRayCount );
	details += wxString::Format("  over box of dimensions: %lg x %lg\n", 
				   rd.SunXMax-rd.SunXMin, rd.SunYMax-rd.SunYMin );
	details += wxString::Format("DNI: %lg\n", DNI);
	details += wxString::Format("Power per ray: %lg\n", PowerPerRay);
	details += wxString::Format("Num plotted rays: %d\n", (int)m_3d->GetNumPlotted() );
	details += wxString::Format("Power of plotted rays: %lg\n", (double)m_3d->GetNumPlotted()*PowerPerRay);

	double cx, cy, cz;
	m_3d->GetCentroid( &cx, &cy, &cz );
	details += wxString::Format("Centroid of Plotted Intersections (x,y,z):\n   ( %lg, %lg, %lg )", cx, cy, cz );

	m_summary->ChangeValue( details );
}

void IntersectionForm::UpdatePlot( )
{
	m_3d->SetZoomRate( 10 ); // higher number will be slower zooming in/out.  default was 20
	m_3d->SetViewZ( -10000, 10000 );
	m_3d->SetCoordSys( m_coordSys->GetSelection() );
	m_3d->SetPointColourMode( m_pointColors->GetSelection() );
	m_3d->ShowFinalOnly( m_finalOnly->GetValue() );
	m_3d->SetRayNumbers( m_plotRays->GetValue() ? wxCommaDashListToIndices( m_rayNumbers->GetValue() ) : std::vector<int>(), m_inclMissedRays->GetValue() );
	m_3d->Invalidate( m_elementList );

	UpdateDetails();
}

void IntersectionForm::PopulateStages()
{
	m_stageList->Freeze();
	m_stageList->Clear();

	if (system == 0) return;

	std::vector<size_t> cursel;
	for( size_t i=0;i<m_stageList->GetCount();i++)
		if (m_stageList->IsChecked(i))
			cursel.push_back( i );

	for (size_t i=0; i<m_prj.StageList.size();i++)
	{
		m_stageList->Append( m_prj.StageList[i]->Name );
		m_stageList->Check( i, 
			std::find( cursel.begin(), cursel.end(), i ) != cursel.end() );
	}
	m_stageList->Thaw();
}

void IntersectionForm::PopulateElements()
{
	m_elementList->Clear();
	
	for (size_t j=0;j< m_prj.StageList.size();j++)
	{
		Stage *stage = m_prj.StageList[j];

		if ( !m_stageList->IsChecked(j) )
			continue;

		for (size_t i=0;i<stage->ElementList.size();i++)
			m_elementList->Add( j, i );
	}

	m_elementList->Invalidate();
	m_elementList->Refresh();
}

void IntersectionForm::CullElementSelections()
{
	for (size_t j=0;j< m_prj.StageList.size();j++)
		if ( !m_stageList->IsChecked(j) )
			m_elementList->UnselectAllStage( j );

	m_elementList->Refresh();
}
