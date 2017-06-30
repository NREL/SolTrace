#ifndef __raydata_h
#define __raydata_h

#include <wx/wx.h>
#include "project.h"

class wxExtGridCtrl;
class wxChoice;

class RayDataForm : public wxPanel
{
public:
	RayDataForm( wxWindow *parent, Project &prj );

	void UpdateDataDisplay();

private:
	Project &m_prj;
	wxExtGridCtrl *m_grid;
	wxChoice *m_coords, *m_stage;

	void CopyDataToClipboard();
	void SaveDataAsCSV();

	void OnCommand( wxCommandEvent & );

	DECLARE_EVENT_TABLE();
};

#endif
