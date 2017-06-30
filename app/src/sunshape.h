#ifndef __sunshape_h
#define __sunshape_h

#include <wx/panel.h>
#include <wx/grid.h>

#include "project.h"

class wxNumericCtrl;
class wxPLPlotCtrl;
class wxRadioChoice;
class wxCheckBox;
class wxExtGridCtrl;
class wxStaticText;
class wxButton;
class wxStaticBox;

class SunShapeForm : public wxPanel
{
public:
	SunShapeForm( wxWindow *parent, Project &prj );	
	void UpdateFromData();

private:
	Project &m_prj;

	void OnCommand( wxCommandEvent & );
	void OnGridCellChange( wxGridEvent & );
	void UpdateCoordinates();
	void UpdatePointSourceWidgets();
	void UpdateWidgetsEnabled();
	void UpdatePlot();
	void UpdateTable();
	
	wxRadioChoice *m_xyzOrLdh;
	wxStaticText *m_lblX, *m_lblY, *m_lblZ;
	wxNumericCtrl *m_x, *m_y, *m_z;
	wxCheckBox *m_pointSrc;
	
	wxStaticBox *m_grpSunShape;
	wxRadioChoice *m_shape;
	wxNumericCtrl *m_sigma, *m_halfWidth;
	wxNumericCtrl *m_userPoints;
	wxExtGridCtrl *m_userData;
	wxStaticText *m_dirLabel;
	wxPLPlotCtrl *m_plot;
	wxButton *m_predefined, *m_load, *m_save, *m_copy, *m_paste;

	void Modified();

	DECLARE_EVENT_TABLE();
};

#endif
