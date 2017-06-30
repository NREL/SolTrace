#ifndef __fluxmap_h
#define __fluxmap_h

#include <wx/wx.h>
#include "project.h"

class wxPLPlotCtrl;
class ElementListBox;
class wxNumericCtrl;

class FluxMapForm : public wxPanel
{
public:
	FluxMapForm( wxWindow *parent, Project &prj );
	virtual ~FluxMapForm();

	void UpdateList();
	void UpdatePlot();
private:
	Project &m_prj;
	ElementStatistics m_es;
	
	ElementListBox *m_elementList;
	wxPLPlotCtrl *m_plot;
	wxNumericCtrl *m_numXBins, *m_numYBins;
	wxNumericCtrl *m_minX, *m_minY, *m_maxX, *m_maxY;
	wxCheckBox *m_autoExtents;
	wxCheckBox *m_finalOnly;
	wxNumericCtrl *m_dni;
	wxStaticText *m_summary;
	wxNumericCtrl *m_contourLevels;
	wxChoice *m_colorScheme;

	void OnCommand( wxCommandEvent & );

	void WriteTecFlx();
	void GetCurrentSelection( int *stage, int *element, char *surfidx = 0 );

	DECLARE_EVENT_TABLE();
};


#endif
