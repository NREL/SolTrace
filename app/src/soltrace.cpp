
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
#include <wx/frame.h>
#include <wx/simplebook.h>
#include <wx/filename.h>
#include <wx/log.h>
#include <wx/stdpaths.h>
#include <wx/busyinfo.h>

#ifdef __WXMSW__
#include <wex/mswfatal.h>
#endif

#include <wex/easycurl.h>
#include <wex/metro.h>
#include <wex/icons/cirplus.cpng>
#include <wex/icons/qmark.cpng>
#include <wex/utils.h>

#include "../resource/menu.cpng"
#include "../resource/notes_white.cpng"

#include "sunshape.h"
#include "optics.h"
#include "geometry.h"
#include "trace.h"
#include "script.h"
#include "intersections.h"
#include "fluxmap.h"
#include "raydata.h"

#include "soltrace.h"

int version_major = 3;
int version_minor = 3;
int version_micro = 0;

class CustomThemeProvider : public wxMetroThemeProvider
{
public:
	virtual ~CustomThemeProvider() { }
	virtual wxColour Colour( int id )
	{
		switch( id )
		{
		case wxMT_FOREGROUND: return wxColour( 130,186,0 );
		case wxMT_HOVER: return wxColour( 0, 138, 23 );
		case wxMT_DIMHOVER : return wxColour( 0, 102, 18 );
		default:
			return wxMetroThemeProvider::Colour( id );
		/*
		case wxMT_BACKGROUND:  return *wxWHITE;
		case wxMT_HOVER: return wxColour( 0, 88, 153 );
		case wxMT_DIMHOVER: return wxColour( 0, 107, 186 );
		case wxMT_LIGHTHOVER: return wxColour( 231, 232, 238 );
		case wxMT_ACCENT: return wxColour( 255, 143, 50 );
		case wxMT_TEXT: return wxColour( 135, 135, 135 ); 
		case wxMT_ACTIVE: return wxColour( 0, 114, 198 );
		case wxMT_SELECT:  return wxColour(193,210,238);
		case wxMT_HIGHLIGHT: return wxColour(224,232,246);
		*/
		}
	}
};

static wxLogWindow *g_logWindow = 0;
class MyLogWindow : public wxLogWindow
{
public:
	MyLogWindow( )	: wxLogWindow( 0, "soltrace-log" ) { 
		GetFrame()->SetPosition( wxPoint( 5, 5 ) );
		GetFrame()->SetClientSize( wxSize(1100,200) );
	}
	virtual bool OnFrameClose( wxFrame * ) {
		g_logWindow = 0; // clear the global pointer, then delete the frame
		return true;
	}
	
	static void Setup()
	{
		if ( g_logWindow != 0 )
			delete g_logWindow;

		g_logWindow = new MyLogWindow;
		wxLog::SetActiveTarget( g_logWindow );
		g_logWindow->Show();
	}
};



static wxArrayString g_appArgs;
class MyApp : public wxApp
{
public:
	virtual void OnFatalException()
	{
#ifdef __WXMSW__
		wxMSWHandleApplicationFatalException();
#endif
	}

	virtual bool OnInit()
	{
#ifdef _DEBUG
		MyLogWindow::Setup();
#endif


		bool is64 = (sizeof(void*) == 8);

#ifdef __WXMSW__
		wxMSWSetupExceptionHandler( "SolTrace", 
			wxString::Format("%d.%d.%d (%d bit)", version_major, version_minor, version_micro, is64 ? 64 : 32 ), 
			"soltrace.support@nrel.gov" );
#endif
		
		wxEasyCurl::Initialize();

		for( int i=0;i<argc;i++ )
			g_appArgs.Add( argv[i] );

		wxInitAllImageHandlers();

		wxMetroTheme::SetTheme( new CustomThemeProvider );
		wxLKScriptWindow::SetFactory( new SolTraceScriptWindowFactory );

		MainWindow *mw = new MainWindow;
		mw->Show();
		if ( g_appArgs.size() > 1 )
			mw->LoadProject( g_appArgs[1] );

		return true;
	}

	virtual int OnExit()
	{
		if ( g_logWindow ) g_logWindow->GetFrame()->Close();
		return wxApp::OnExit();
	}
};

#ifndef ST_CONSOLE_APP
    IMPLEMENT_APP( MyApp );
#endif

enum { ID_MAIN_MENU = wxID_HIGHEST+123, ID_TABS,
	ID_NEW_SCRIPT, ID_OPEN_SCRIPT };



BEGIN_EVENT_TABLE( MainWindow, wxFrame )
	EVT_BUTTON( ID_MAIN_MENU, MainWindow::OnCommand)
	EVT_LISTBOX( ID_TABS, MainWindow::OnCaseTabChange )
	EVT_BUTTON( ID_TABS, MainWindow::OnCaseTabButton )
	EVT_CLOSE( MainWindow::OnClose )
	EVT_MENU( wxID_ABOUT, MainWindow::OnCommand )
	EVT_MENU( wxID_HELP, MainWindow::OnCommand )
	EVT_MENU( wxID_NEW, MainWindow::OnCommand )
	EVT_MENU( wxID_OPEN, MainWindow::OnCommand)
	EVT_MENU( wxID_SAVE, MainWindow::OnCommand )
	EVT_MENU( wxID_SAVEAS, MainWindow::OnCommand )
	EVT_MENU( ID_NEW_SCRIPT, MainWindow::OnCommand )
	EVT_MENU( ID_OPEN_SCRIPT, MainWindow::OnCommand )
	EVT_MENU( wxID_CLOSE, MainWindow::OnCommand )
	EVT_MENU( wxID_EXIT, MainWindow::OnCommand )
	EVT_BUTTON( wxID_HELP, MainWindow::OnCommand )
END_EVENT_TABLE()

static MainWindow *g_mainWindow = 0;

MainWindow::MainWindow()
	: wxFrame( 0, wxID_ANY, "SolTrace", 
		wxDefaultPosition, wxSize( 1100, 700 ) )
{
#ifdef __WXMSW__
	SetIcon( wxICON( appicon ) );
#endif

	if ( g_mainWindow != 0 ) 
		wxMessageBox("internal error - only one main window can exist!");
	else
		g_mainWindow = this;

	wxBoxSizer *tools = new wxBoxSizer( wxHORIZONTAL );
	tools->Add( m_mainMenuButton = new wxMetroButton( this, ID_MAIN_MENU, wxEmptyString, wxBITMAP_PNG_FROM_DATA( menu ), wxDefaultPosition, wxDefaultSize /*, wxMB_DOWNARROW */), 0, wxALL|wxEXPAND, 0 );
	m_tabList = new wxMetroTabList( this, ID_TABS, wxDefaultPosition, wxDefaultSize );
	tools->Add( m_tabList, 1, wxALL|wxEXPAND, 0 );		
	tools->Add( new wxMetroButton( this, wxID_HELP, wxEmptyString, wxBITMAP_PNG_FROM_DATA(qmark), wxDefaultPosition, wxDefaultSize), 0, wxALL|wxEXPAND, 0 );
	
	m_notebook = new wxSimplebook( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxBORDER_NONE );
		
	wxBoxSizer *sizer = new wxBoxSizer( wxVERTICAL );
	sizer->Add( tools, 0, wxALL|wxEXPAND, 0 );
	sizer->Add( m_notebook, 1, wxALL|wxEXPAND, 0 );
	SetSizer(sizer);
		
	std::vector<wxAcceleratorEntry> entries;
	entries.push_back( wxAcceleratorEntry( wxACCEL_CTRL, 'o', wxID_OPEN ) );
	entries.push_back( wxAcceleratorEntry( wxACCEL_CTRL, 's', wxID_SAVE ) );
	entries.push_back( wxAcceleratorEntry( wxACCEL_CTRL | wxACCEL_SHIFT, 'n', ID_NEW_SCRIPT ) );
	entries.push_back( wxAcceleratorEntry( wxACCEL_CTRL | wxACCEL_SHIFT, 'o', ID_OPEN_SCRIPT ) );
	entries.push_back( wxAcceleratorEntry( wxACCEL_CTRL, 'w', wxID_CLOSE ) );
	entries.push_back( wxAcceleratorEntry( wxACCEL_NORMAL, WXK_F1, wxID_HELP ) );
	SetAcceleratorTable( wxAcceleratorTable( entries.size(), &entries[0] ) );

	m_tabList->Append( "Sun" );
	m_sunShapeForm = new SunShapeForm( m_notebook, m_project );
	m_notebook->AddPage( m_sunShapeForm, "Sun" );
	
	m_tabList->Append( "Optics" );
	m_opticsForm = new OpticsForm( m_notebook, m_project );
	m_notebook->AddPage( m_opticsForm, "Optics" );
	
	m_tabList->Append( "Geometry" );
	m_geometryForm = new GeometryForm( m_notebook, m_project );
	m_notebook->AddPage( m_geometryForm, "Geometry" );
	
	m_tabList->Append( "Trace" );
	m_traceForm = new TraceForm( m_notebook, m_project );
	m_notebook->AddPage( m_traceForm, "Trace" );
	
	m_tabList->Append( "Intersections" );
	m_intersectionForm = new IntersectionForm( m_notebook, m_project );
	m_notebook->AddPage( m_intersectionForm, "Intersections (3D)" );

	m_tabList->Append( "Flux maps" );
	m_fluxMapForm = new FluxMapForm( m_notebook, m_project );
	m_notebook->AddPage( m_fluxMapForm, "Flux maps" );

	m_tabList->Append( "Data" );
	m_rayDataForm = new RayDataForm( m_notebook, m_project );
	m_notebook->AddPage( m_rayDataForm, "Data" );

	m_modified = false;
	m_currentpage = "Sun";

	UpdateFrameTitle();
}

MainWindow::~MainWindow() {
	g_mainWindow = 0;
}

MainWindow &MainWindow::Instance()
{
	if ( !g_mainWindow )
	{
		wxMessageBox("internal error - no mainwindow instance defined" );
		::exit( -1 );
	}

	return *g_mainWindow;
}

void MainWindow::SetModified( bool b )
{
	m_modified = b;
	UpdateFrameTitle();
}

wxString MainWindow::GetWorkDir()
{
	return m_traceForm->GetWorkDir();
}

wxString MainWindow::GetAppDataDir()
{
	static wxString path;
	if ( path.IsEmpty() )
	{
		wxFileName dir( wxPathOnly(wxStandardPaths::Get().GetExecutablePath()) + "/.." );
		dir.Normalize();
		path = dir.GetFullPath();
	}
	return path;		
}
		
bool MainWindow::LoadProject( const wxString &file, bool quiet )
{	
	FILE *fp = fopen( file.c_str(), "r" );
	if (!fp)
	{
		if (!quiet) wxMessageBox( "Failed to open file:\n" + file );
		m_fileName.clear();
		return false;
	}

	wxBusyInfo info( "Loading project " + wxFileNameFromPath(file) + "..." );

	bool ok = m_project.Read(fp);
	if (ok) m_fileName = file;
	else if (!quiet) wxMessageBox( "Invalid file format for loading:\n" + file);

	fclose(fp);

	UpdateAllInputForms();
	UpdateResults();

	SetModified(false);
	return ok;
}
		
bool MainWindow::SaveProject( const wxString &file, bool quiet )
{
	
	if (!quiet &&
		wxFileName(file).GetPath()
		== wxFileName( GetAppDataDir() + "/samples/").GetPath())
	{
		if ( wxID_NO == wxMessageBox("You are attempting to write to a standard sample file included with SolTrace.  "
			"Are you sure you want to overwrite it?\n\n" + wxFileName(file).GetName(), "Notice", wxYES_NO ) )
		{
			return false;
		}
	}
	
	FILE *fp = fopen( file.c_str(), "w" );
	if (!fp) return false;
	
	wxBusyInfo info( "Saving project " + wxFileNameFromPath(file) + "..." );
	m_project.Write( fp );
	fclose(fp);

	m_fileName = file;
	SetModified( false );

	return true;
}
		
void MainWindow::Save()
{
	if ( m_fileName.IsEmpty() )
	{
		SaveAs();
		return;
	}

	if ( !SaveProject( m_fileName ) )
		wxMessageBox("Error writing project to disk:\n\n" + m_fileName, "Notice", wxOK, this );		
}
		
void MainWindow::SaveAs()
{
	wxFileDialog dlg( this, "Save SolTrace input file as", wxPathOnly(m_fileName), 
		m_fileName, "SolTrace Project File (*.stinput)|*.stinput", 
		wxFD_SAVE|wxFD_OVERWRITE_PROMPT );
	if ( dlg.ShowModal() == wxID_OK )
	{
		m_fileName = dlg.GetPath();
		Save();
	}
}

bool MainWindow::CloseProject( bool force )
{	
	if ( IsModified() && !force )
	{
		wxString name( m_fileName );
		if ( name.IsEmpty() ) name = "untitled";
		int ret = wxMessageBox("The project '" + wxFileName(name).GetName() + "' has been modified.  Save changes?", "Query", wxICON_EXCLAMATION|wxYES_NO|wxCANCEL, this );
		if (ret == wxYES)
		{
			Save( );
			if ( IsModified() ) // if failed to save, cancel
				return false;
		}
		else if (ret == wxCANCEL)
			return false;
	}


	m_project.Sun.ResetToDefaults();
	m_project.ClearOptics();
	m_project.ClearStages();	
	m_project.Results.FreeMemory();

	m_fileName.Clear();
	SetModified( false );

	UpdateAllInputForms();
	UpdateResults();
	UpdateFrameTitle();
	return true;

}

void MainWindow::OnClose( wxCloseEvent &evt )
{
	if ( !SolTraceScriptWindow::CloseAll() )
	{
		evt.Veto();
		return;
	}

	Raise();
	if ( !CloseProject() )
	{
		evt.Veto();
		return;
	}
	
	// destroy the window
#ifndef ST_CONSOLE_APP
	wxGetApp().ScheduleForDestruction( this );
#endif
}

void MainWindow::OnCommand( wxCommandEvent &evt )
{
	switch( evt.GetId() )
	{
	case wxID_NEW:
	case wxID_CLOSE:
		CloseProject();
		break;
	case wxID_SAVEAS:
		SaveAs();
		break;
	case wxID_SAVE:
		Save();
		break;
	case wxID_OPEN:
		{
			if ( !CloseProject() ) return;			
			wxFileDialog dlg(this, "Open SolTrace input file", wxEmptyString, wxEmptyString, "SolTrace Input Files (*.stinput)|*.stinput", wxFD_OPEN );
			if (dlg.ShowModal() == wxID_OK)
				if( !LoadProject( dlg.GetPath(), false ) )
					wxMessageBox("Error loading input file:\n\n" 
						+ dlg.GetPath() + "\n\n" /*+ m_project.GetLastError()*/, "Notice", wxOK, this );
		}
		break;

	case ID_NEW_SCRIPT:
		SolTraceScriptWindow::CreateNewWindow();
		break;

	case ID_OPEN_SCRIPT:
		SolTraceScriptWindow::OpenFiles();
		break;

	case ID_MAIN_MENU:
		{
			wxPoint p = m_mainMenuButton->ClientToScreen( wxPoint( 0, m_mainMenuButton->GetClientSize().y ) );
			wxMetroPopupMenu menu;
			menu.Append( wxID_NEW, "New project\tCtrl-N" );
			menu.Append( ID_NEW_SCRIPT, "New script" );
			menu.AppendSeparator();
			menu.Append( wxID_OPEN, "Open project\tCtrl-O" );
			menu.Append( ID_OPEN_SCRIPT, "Open script" );
			menu.AppendSeparator();
			menu.Append( wxID_SAVE, "Save\tCtrl-S" );
			menu.Append( wxID_SAVEAS, "Save as..." );
			menu.AppendSeparator();
			menu.Append( wxID_CLOSE, "Close\tCtrl-W" );
			menu.Append( wxID_EXIT, "Quit" );
			menu.Popup( this, p );
		}
		break;
	case wxID_EXIT:
		Close();
		break;
	case wxID_HELP:
		ShowHelpTopic( m_currentpage );
		break;

	};
}

void MainWindow::OnCaseTabChange( wxCommandEvent &evt )
{
	m_notebook->SetSelection( evt.GetSelection() );
	m_currentpage = m_notebook->GetPageText(evt.GetSelection());
}

void MainWindow::OnCaseTabButton( wxCommandEvent & )
{
	// not used currently
}

void MainWindow::UpdateFrameTitle()
{
	bool is64 = (sizeof(void*)==8);
	wxString title = "SolTrace " + wxString::Format("%d.%d.%d (%d bit)", version_major, version_minor, version_micro, is64 ? 64 : 32 );
	if ( !m_fileName.IsEmpty() )	title += ": " + m_fileName;
	else title += ": untitled";

	if ( m_modified ) title += " *";

	if ( GetTitle() != title )
		SetTitle( title );
}

void MainWindow::UpdateAllInputForms()
{
	m_sunShapeForm->UpdateFromData();
	m_opticsForm->UpdateList( 0 );
	m_geometryForm->UpdateForm();
}

void MainWindow::UpdateResults()
{
	m_intersectionForm->UpdateView();
	m_fluxMapForm->UpdateList();
	m_fluxMapForm->UpdatePlot();
	m_rayDataForm->UpdateDataDisplay();
}

void MainWindow::ShowHelpTopic( const wxString &topic )
{
	wxString topic_map;
	if (topic == "Sun")
		topic_map = "Sun";
	else if (topic == "Optics")
		topic_map = "OpticalProperties";
	else if (topic == "Geometry")
		topic_map == "Geometry";
	else if (topic == "Trace")
		topic_map = "Tracing";
	else if ((topic == "Intersections (3D)") || (topic == "Flux maps"))
		topic_map = "Visualization";
	else if (topic == "Data")
		topic_map = "ExportingData";
	else
		topic_map = "Introduction";
	
	try
	{
		wxFileName help_dir;
		help_dir.SetPath(wxPathOnly(wxStandardPaths::Get().GetExecutablePath()) + "/..");
		help_dir.AppendDir("help");
		try
		{
			wxString help = wxString::Format("HH.EXE \"ms-its:%s/SolTrace.chm::/%s.htm\"", help_dir.GetPath().ToStdString(), topic_map);
			wxExecute(help);
		}
		catch (...)
		{
			//wxFileName fn(MainWindow::Instance().GetAppDataDir() + "/help/59163.pdf");
			wxFileName fn( MainWindow::Instance().GetAppDataDir() + "SolTrace.pdf" );
			fn.Normalize( );
			wxLaunchDefaultBrowser( "file:///" + fn.GetFullPath( ) );
		}
	}
	catch (...)
	{
		wxMessageBox("Sorry, help is not available right now.");
	}
}
