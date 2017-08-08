#ifndef __soltrace_h
#define __soltrace_h

#include <wx/frame.h>
#include "project.h"

#define ST_CONSOLE_APP      //define this if we want to compile as a console application (no gui)

class wxSimplebook;
class wxPanel;
class wxMetroButton;
class wxMetroTabList;
class SunShapeForm;
class OpticsForm;
class GeometryForm;
class TraceForm;
class RayDataForm;
class FluxMapForm;
class IntersectionForm;

class MainWindow : public wxFrame
{
public:
	MainWindow();
	virtual ~MainWindow();

	static MainWindow &Instance();

	static void ShowHelpTopic( const wxString &topic );
		
	bool LoadProject( const wxString &file, bool quiet = false );
	bool SaveProject( const wxString &file, bool quiet = false );

	void Save();
	void SaveAs();
	bool CloseProject( bool force = false );

	bool IsModified() { return m_modified; }
	void SetModified( bool b = true );

	wxString GetFileName() { return m_fileName; }
	wxString GetWorkDir();
	wxString GetAppDataDir();
	
	void UpdateFrameTitle();
	void UpdateAllInputForms();
	void UpdateResults();

	Project &GetProject() { return m_project; }
	SunShapeForm *GetSunShape() { return m_sunShapeForm; }
	OpticsForm *GetOptics() { return m_opticsForm; }
	GeometryForm *GetGeometry() { return m_geometryForm; }
	TraceForm *GetTrace() { return m_traceForm; }
	FluxMapForm *GetFluxMaps() { return m_fluxMapForm; }
	RayDataForm *GetRayData() { return m_rayDataForm; }
	IntersectionForm *GetIntersectionForm() { return m_intersectionForm; }
	
protected:
	void OnClose( wxCloseEvent & );
	void OnCommand( wxCommandEvent & );
	void OnCaseTabChange( wxCommandEvent & );
	void OnCaseTabButton( wxCommandEvent & );
	
private:
	wxMetroButton *m_mainMenuButton;
	wxMetroTabList *m_tabList;
	wxSimplebook *m_notebook;

	SunShapeForm *m_sunShapeForm;
	OpticsForm *m_opticsForm;
	GeometryForm *m_geometryForm;
	TraceForm *m_traceForm;
	IntersectionForm *m_intersectionForm;
	FluxMapForm *m_fluxMapForm;
	RayDataForm *m_rayDataForm;

	bool m_modified;
	Project m_project;
	wxString m_fileName;
	
	DECLARE_EVENT_TABLE();
};


#endif
