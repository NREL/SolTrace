
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
#include <wx/renderer.h>
#include <wx/dcbuffer.h>
#include <wex/utils.h>

#include "elementlist.h"

wxString ElementItem::Label( Project &prj )
{
	if ( Element *e = prj.GetElement( stage, element ) )
	{
		wxString name( wxString::Format("%d,%d: ", (int)(stage+1), (int)(element+1) ) );
		if (!e->Comment.IsEmpty()) name += e->Comment;
		else name += wxString("Surf.'") + e->SurfaceIndex + wxString("'");

		name += wxString::Format(" (%d hits)", e->RayHits);
		return name;
	}
	else
		return "<index error>";
}


BEGIN_EVENT_TABLE(ElementListBox, wxScrolledWindow)
	EVT_LEFT_DOWN( ElementListBox::OnLeftDown )
	EVT_PAINT( ElementListBox::OnPaint )
	EVT_ERASE_BACKGROUND( ElementListBox::OnErase )
	EVT_SIZE( ElementListBox::OnResize )	
END_EVENT_TABLE()

ElementListBox::ElementListBox( Project &prj, wxWindow *parent, int id, 
	const wxPoint &pos, const wxSize &size)
	: wxScrolledWindow( parent, id, pos, size, wxCLIP_CHILDREN|wxBORDER_SIMPLE ), m_prj(prj)
{
	SetBackgroundStyle( wxBG_STYLE_CUSTOM );
	m_bestSize.Set(100,100);
	m_scrollRate = 25;
	m_space = (int)( 4.0 * wxGetScreenHDScale() );
	m_itemHeight = 20;  // default. gets set in Invalidate() properly.
	m_chkBoxSize = wxRendererNative::Get().GetCheckBoxSize( this );

	SetBackgroundColour( *wxWHITE );
	SetFont( *wxNORMAL_FONT );
}

bool ElementListBox::IsSelected( size_t idx )
{
	if ( idx < m_items.size() )
		return m_sel[ Code( m_items[idx].stage, m_items[idx].element ) ];
	else return false;
}

bool ElementListBox::IsSelected( unsigned int s, unsigned int e )
{
	return m_sel[ Code(s,e) ];
}

bool ElementListBox::GetStageElementIndices( size_t i, unsigned int *stage, unsigned int *elem )
{
	if ( i < m_items.size() )
	{
		*stage = m_items[i].stage;
		*elem = m_items[i].element;
		return true;
	}
	else return false;

}


unsigned int ElementListBox::Code( unsigned int s, unsigned int e )
{
	return ( s << 24 ) | ( 0x00FFFFFF & e );
}

void ElementListBox::Select( unsigned int s, unsigned int e, bool tf )
{
	m_sel[ Code(s,e) ] = tf;
}

void ElementListBox::SelectAll()
{
	m_sel.clear();

	for( size_t i=0;i<m_items.size();i++ )
		Select( m_items[i].stage, m_items[i].element );

	Refresh();
}

void ElementListBox::UnselectAllStage( unsigned int stage )
{
	for( std::unordered_map<unsigned int,bool>::iterator it = m_sel.begin();
		it != m_sel.end();
		++it )
	{
		unsigned int S = it->first;
		S = S >> 24;
		if ( stage == S )
			it->second = false;
	}
}

void ElementListBox::ClearSelections()
{
	m_sel.clear();
	Refresh();
}

void ElementListBox::Clear() { m_items.clear(); Invalidate(); }
void ElementListBox::Reserve( size_t n ) { m_items.reserve( n ); }
void ElementListBox::Add( unsigned int stage, unsigned int element ) { m_items.push_back( ElementItem(stage,element) ); }
size_t ElementListBox::Count() { return m_items.size(); }

void ElementListBox::RecalculateBestSize()
{
	wxClientDC dc(this);
	dc.SetFont( GetFont() );
	m_itemHeight = dc.GetCharHeight();
	if ( m_itemHeight < m_chkBoxSize.y )
		m_itemHeight = m_chkBoxSize.y;
	m_itemHeight += m_space;
		
	m_bestSize = GetClientSize();
	m_bestSize.y = m_items.size() * m_itemHeight; 
}

void ElementListBox::OnPaint( wxPaintEvent &evt )
{
	wxAutoBufferedPaintDC dc(this);
	DoPrepareDC(dc);

	int cwidth = 0, cheight = 0;
	GetClientSize( &cwidth, &cheight );

	wxSize client( GetClientSize()  );
	
	// paint the background.
	wxColour bg = GetBackgroundColour();
	dc.SetBrush(wxBrush(bg));
	dc.SetPen(wxPen(bg,1));
	wxRect windowRect( wxPoint(0,0), client );
	CalcUnscrolledPosition(windowRect.x, windowRect.y,
		&windowRect.x, &windowRect.y);
	dc.DrawRectangle(windowRect);


	dc.SetFont( GetFont() );
	dc.SetTextForeground( *wxBLACK );
	int text_height = dc.GetCharHeight();
	wxRendererNative &rndr( wxRendererNative::Get() );
	for (int i=0;i<m_items.size();i++)	
	{
		int py_start = i*m_itemHeight;

		if ( py_start >= windowRect.y - m_itemHeight
			&& py_start <= windowRect.y + client.y+m_itemHeight )
		{
			wxRect chkrct( m_space/2, py_start + m_itemHeight/2 - m_chkBoxSize.y/2, m_chkBoxSize.x, m_chkBoxSize.y );
			rndr.DrawCheckBox( this, dc, chkrct, IsSelected(i) ? wxCONTROL_CHECKED : 0 );
			dc.DrawText( m_items[i].Label(m_prj), m_space + m_chkBoxSize.x, 
				py_start + m_itemHeight/2 - text_height/2 );
		}
	}			
}
		
void ElementListBox::OnLeftDown(wxMouseEvent &evt)
{
	int vsx, vsy;
	GetViewStart( &vsx, &vsy );
	vsx *= m_scrollRate;
	vsy *= m_scrollRate;
	int px =evt.GetY()+vsy;

	SetFocus();

	int clickindex = px/m_itemHeight;
	if ( clickindex >= 0 && clickindex < (int)m_items.size() )
	{
		wxCommandEvent selevt(wxEVT_COMMAND_LISTBOX_SELECTED, this->GetId() );
		selevt.SetEventObject(this);
		selevt.SetInt(clickindex);
		selevt.SetString( m_items[clickindex].Label(m_prj) );
		GetEventHandler()->ProcessEvent(selevt);
		return;
	}

}
	
void ElementListBox::OnErase(wxEraseEvent &evt)
{
	/* nothing to do */
}

	
void ElementListBox::OnResize(wxSizeEvent &evt)
{
	ResetScrollbars();
}
		
void ElementListBox::Invalidate()
{
	RecalculateBestSize();
	ResetScrollbars();
	InvalidateBestSize();
}

void ElementListBox::ResetScrollbars()
{
	int hpos, vpos;
	GetViewStart( &hpos, &vpos );
	SetScrollbars( m_scrollRate, m_scrollRate, m_bestSize.GetWidth()/m_scrollRate,m_bestSize.GetHeight()/m_scrollRate,hpos,vpos );	
	InvalidateBestSize();
}

wxSize ElementListBox::DoGetBestSize() const
{
	return m_bestSize;
}
