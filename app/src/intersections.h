
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
    wxNumericCtrl *m_scaleX, *m_scaleY, *m_scaleZ, *m_scaleZoom;
	
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
