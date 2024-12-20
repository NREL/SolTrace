
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


#ifndef __soltrace_h
#define __soltrace_h

#include <wx/frame.h>
#include "project.h"

//#define ST_CONSOLE_APP      //define this if we want to compile as a console application (no gui)

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
	wxString m_currentpage;
	
	DECLARE_EVENT_TABLE();
};


#endif
