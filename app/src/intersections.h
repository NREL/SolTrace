#ifndef __intersection_h
#define __intersection_h

#include <wx/window.h>
#include <wex/gleasy.h>

#include "project.h"

class wxNumericCtrl;
class wxCheckBox;
class wxCheckListBox;
class wxChoice;
class ElementListBox;


class IntersectionViewer : public wxGLEasyCanvas
{
public:
	IntersectionViewer( wxWindow *parent, Project &prj );
	virtual ~IntersectionViewer();

	void SetupAxes( bool show = true, bool ticks = true );
	void Invalidate( ElementListBox *lb = 0);

	void Clear();
	void SetCoordSys( int coord ) { m_coordSys = coord; }
	void SetPointColourMode( int mode ) { m_pointColors = mode; }
	void ShowFinalOnly( bool b ) { m_finalOnly = b; }
	void SetRayNumbers( const std::vector<int> &rnums, bool incl_missed ) { m_rayNumbers = rnums; m_includeMissedRays = incl_missed; }

	void FitView();

	size_t GetNumPlotted() { return m_nplotted; }
	void GetCentroid( double *x, double *y, double *z ) {
		*x = m_centroid[0]; *y = m_centroid[1]; *z = m_centroid[2];
	}

	void RebuildGeometry( ElementListBox *lb );

private:
	virtual void OnRender();
	
	Project &m_prj;

	GLuint m_vertexListId;
	bool m_vertexListValid;
	wxGLPoint3D m_min, m_max;

	int m_coordSys, m_pointColors;
	bool m_finalOnly;
	std::vector<int> m_rayNumbers;
	bool m_includeMissedRays;

	std::vector<wxColour> m_colorList;
	
	size_t m_nplotted;
	double m_centroid[3];

	bool m_showAxes, m_showTicks;
};


class IntersectionForm : public wxWindow
{
public:
	IntersectionForm( wxWindow *parent, Project &prj );
	virtual ~IntersectionForm();

	void UpdateView();

private:
	Project &m_prj;

	wxChoice *m_coordSys;
	wxCheckListBox *m_stageList;
	ElementListBox *m_elementList;
	wxCheckBox *m_finalOnly;
	wxNumericCtrl *m_dni;
	wxCheckBox *m_plotRays, *m_inclMissedRays;
	wxTextCtrl *m_rayNumbers;
	wxTextCtrl *m_summary;
	
	wxChoice *m_pointColors;
	wxSlider *m_axisTextSize;
	wxCheckBox *m_showAxes, *m_showTicks;
	wxNumericCtrl *m_scaleX, *m_scaleY, *m_scaleZ;
	
	std::vector<wxColour> m_colourList;
	
	IntersectionViewer *m_3d;

	void OnCommand( wxCommandEvent & );

	
	void UpdateDetails();
	void UpdatePlot();
	void PopulateStages();
	void PopulateElements();
	bool GetStageElementIndices( size_t i, size_t *stage, size_t *elem );
	bool IsSelected( size_t stage_num, size_t element_num );
	void CullElementSelections();

	DECLARE_EVENT_TABLE();
};

#endif
