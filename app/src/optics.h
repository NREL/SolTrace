#ifndef __optics_h
#define __optics_h

#include <wx/wx.h>
#include <wx/panel.h>
#include <wx/grid.h>

#include "project.h"

class wxExtTextCtrl;
class wxNotebook;

class OpticalPropertyForm;
class OpticsForm : public wxPanel
{
public:
	OpticsForm( wxWindow *parent, Project &prj );

	void UpdateList(int sel=-1);
	void UpdateOptForms(int idx=-1);

	void Modified();

	void AddOptic(const wxString &name = "New optic");
	void DeleteOptic( int sel );
	void ClearOptics();

private:
	Project &m_prj;

	void OnButton(wxCommandEvent &evt);
	void OnListSelect(wxCommandEvent &evt);
	void OnNameChange(wxCommandEvent &evt);

	wxNotebook *m_optTabs;
	wxListBox *m_opticList;
	OpticalPropertyForm *m_frontOpt;
	OpticalPropertyForm *m_backOpt;
	wxExtTextCtrl *m_optName;

	DECLARE_EVENT_TABLE()
};

#endif