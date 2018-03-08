
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


#ifndef __geometry_h
#define __geometry_h

#include <wx/wx.h>
#include <wx/panel.h>
#include <wx/grid.h>

#include "project.h"

class wxNotebook;

class wxExtGridCtrl;
class wxExtTextCtrl;
class wxNumericCtrl;

class ElementTable;
class GeometryForm;

class StageForm : public wxPanel
{
public:
	StageForm( wxWindow *parent, GeometryForm *geo, Project &prj, Stage *s );

	void UpdateFromData();

	void Import();
	void Export();
	void Clear();
	void Delete(int idx=-1);
	void Append(int n=-1);
	void Insert();

	void EditStageZRot();
	void EditZRot(int idx=-1);
	void EditAperture(int idx=-1);
	void EditSurface(int idx=-1);
	void EditOptics(int idx=-1);

	void Modified();
	void NameModified();

private:

	void UpdateProperties();
	Project &m_prj;
	Stage *m_stage;

	GeometryForm *m_geoForm;

	wxExtTextCtrl *m_stageName;
	wxCheckBox *m_virtualStage, *m_multipleHits, *m_traceThrough;
	wxNumericCtrl *m_x, *m_y, *m_z, *m_ax, *m_ay, *m_az, *m_zrot;

	ElementTable *m_gridTable;
	wxExtGridCtrl *m_grid;


	void OnGridCellChange(wxGridEvent &evt);
	void OnGridCellSelect(wxGridEvent &evt);
	void OnGridRangeSelect(wxGridRangeSelectEvent &evt);
	void OnGridCellDClick(wxGridEvent &evt);
	void OnGridCellRightClick(wxGridEvent &evt);
	
	void OnCommand( wxCommandEvent &evt );

	DECLARE_EVENT_TABLE()
};

class GeometryForm : public wxPanel
{
public:
	GeometryForm( wxWindow *parent, Project &prj );

	void NewStage(const wxString &name = "New stage", int pos = -1 );
	void DeleteStage(int sel, bool quiet=false);
	void ClearStages(bool quiet=false);
	StageForm *GetStageForm(int idx);
	StageForm *GetStageForm( Stage *stg );

	void UpdateForm();
	void UpdateStageNames();
	void UpdateStageFlagsAndGeom(int stage_num);
	
	void Modified();
private:
	Project &m_prj;
	wxNotebook *m_stageTabs;

	void OnButton(wxCommandEvent &evt);

	DECLARE_EVENT_TABLE()
};

#endif
