#include <math.h>
#include "types.h"
#include "procs.h"

#define sqr(x) ((x)*(x))
				
bool AperturePlane(
			TElement *Element )
{
/*{Calculates the aperture plane of the element in element coord system.  Applicable
 to rotationally symmetric apertures surfaces with small curvature: g, s, p, o, c, v, m, e, r, i.
  input - Element = Element record containing geometry of element
  output -
         - Element.ZAperture  where ZAperture is the distance from the origin to the plane.
}*/
	double radius2 = 0.0;
	double FocalLength = 0.0;
	double SumOfAlphaTerms = 0.0;
	double lastz = 0.0;
	double alphas = 0.0;
	double rho2i = 0.0;
	double MaximumRadius=0.0, RadiusTemp=0.0;
	double zm = 0.0;

  //This routine is poorly written (I don't know what I was thinking).  But it works, so I haven't reduced the redundant code,
  //but right here I am doing a special thing for the case of irregular triangular or quadrilateral apertures so I don't have
  //to add MORE redundant code.
	switch ( Element->ShapeIndex )
	{
	case 'i':
	case 'I':
			MaximumRadius = sqrt(sqr(Element->ParameterA)+sqr(Element->ParameterB));
			RadiusTemp = sqrt(sqr(Element->ParameterC)+sqr(Element->ParameterD));
			if (RadiusTemp > MaximumRadius ) MaximumRadius = RadiusTemp;
			RadiusTemp = sqrt(sqr(Element->ParameterE)+sqr(Element->ParameterF));
			if (RadiusTemp > MaximumRadius ) MaximumRadius = RadiusTemp;
		break;
	
	case 'q':
	case 'Q':
			MaximumRadius = sqrt(sqr(Element->ParameterA)+sqr(Element->ParameterB));
			RadiusTemp = sqrt(sqr(Element->ParameterC)+sqr(Element->ParameterD));
			if (RadiusTemp > MaximumRadius ) MaximumRadius = RadiusTemp;
			RadiusTemp = sqrt(sqr(Element->ParameterE)+sqr(Element->ParameterF));
			if (RadiusTemp > MaximumRadius ) MaximumRadius = RadiusTemp;
			RadiusTemp = sqrt(sqr(Element->ParameterG)+sqr(Element->ParameterH));
			if (RadiusTemp > MaximumRadius ) MaximumRadius = RadiusTemp;
		break;
	}


	switch ( Element->SurfaceIndex )
	{
// *************************************************************************
    //   {Surface described by general Spencer & Murty formulation}
	case 'g':
	case 'G':
		
			switch ( Element->ShapeIndex )
			{
			case 'c':
			case 'C':
			case 'h':
			case 'H':
			case 't':
			case 'T': // {circle, hexagon, triangle}
					radius2 = sqr(Element->ParameterA/2.0);
					SumOfAlphaTerms = 0.0;
					rho2i = 1.0;
					for (st_uint_t i=0;i<5;i++)
					{
						alphas = Element->Alpha[i];
						rho2i = rho2i*radius2;
						SumOfAlphaTerms += alphas*rho2i;
					}
					
					Element->ZAperture = Element->VertexCurvX*radius2/(1+sqrt(1-Element->Kappa*Element->VertexCurvX*Element->VertexCurvX*radius2))+SumOfAlphaTerms;
				break;
				
			case 'r':
			case 'R': // {Rectangle}
					radius2 = sqr(Element->ParameterA/2.0) + sqr(Element->ParameterB/2.0);
					SumOfAlphaTerms = 0.0;
					for(st_uint_t i=0;i<5;i++)
						SumOfAlphaTerms += Element->Alpha[i]*pow(radius2, (int)(i+1));
					
					Element->ZAperture = Element->VertexCurvX*radius2/(1+sqrt(1-Element->Kappa*sqr(Element->VertexCurvX)*radius2))+SumOfAlphaTerms;
				break;
	
			case 'a':
			case 'A':
			case 'l':
			case 'L':    //annulus or single axis curvature
					if (fabs(Element->ParameterB) >= fabs(Element->ParameterA))
						radius2 = sqr(Element->ParameterB);
					else
						radius2 = sqr(Element->ParameterA);
						
					SumOfAlphaTerms = 0.0;
					for (st_uint_t i=0;i<5;i++)
						SumOfAlphaTerms =  SumOfAlphaTerms + Element->Alpha[i]*pow(radius2, (int)(i+1));
					
					Element->ZAperture = Element->VertexCurvX*radius2/(1+sqrt(1-Element->Kappa*sqr(Element->VertexCurvX)*radius2))+SumOfAlphaTerms;
				break;
			
			case 'i':
			case 'I':
			case 'q':
			case 'Q': //irregular triangle or quadrilateral
					radius2 = sqr(MaximumRadius);
					
					///// *** NEXT line added by Aron Dobos 4 Nov 09, b/c otherwise SumOfAlphaTerms would be uninitialized...
					SumOfAlphaTerms = 0.0;
					Element->ZAperture = Element->VertexCurvX*radius2/(1+sqrt(1-Element->Kappa*sqr(Element->VertexCurvX)*radius2))+SumOfAlphaTerms;
				break;
			} // end nested switch (Element->ShapeIndex)
		break;
		
// *************************************************************************	
	case 's':
	case 'S': // {Surface described by sphere}
			FocalLength = 1.0/(2.0*Element->VertexCurvX);  //VertexY not used
			switch( Element->ShapeIndex )
			{
			case 'c':
			case 'C':
			case 'h':
			case 'H':
			case 't':
			case 'T': // {circle, hexagon, triangle}
					Element->ZAperture = 2.0*FocalLength - sqrt(sqr(2.0*FocalLength) - sqr(Element->ParameterA/2.0));
				break;
				
			case 'r':
			case 'R': // {Rectangle}
					radius2 = sqr(Element->ParameterA/2.0) + sqr(Element->ParameterB/2.0);
					Element->ZAperture = 2.0*FocalLength - sqrt(sqr(2.0*FocalLength) - radius2);
				break;
			
			case 'a':
			case 'A':
			case 'l':
			case 'L':	   //annulus or single axis curvature				
					if (fabs(Element->ParameterB) >= fabs(Element->ParameterA))
						radius2 = sqr(Element->ParameterB);
					else
						radius2 = sqr(Element->ParameterA);
						
					Element->ZAperture = 2.0*FocalLength - sqrt(sqr(2.0*FocalLength) - radius2);
				break;
			
			case 'i':
			case 'I':
			case 'q':
			case 'Q': //irregular triangle or quadrilateral
					radius2 = sqr(MaximumRadius);
					Element->ZAperture = 2.0*FocalLength - sqrt(sqr(2.0*FocalLength) - radius2);
				break;			
			}
			
		break;
		
// *************************************************************************
	case 'p':
	case 'P': //{Surface described by Parabola}
			FocalLength = 1.0/(2.0*Element->VertexCurvX);
			switch (Element->ShapeIndex)
			{
			case 'c':
			case 'C':
			case 'h':
			case 'H':
			case 't':
			case 'T'://	{circle, hexagon, triangle}
					Element->ZAperture = sqr(Element->ParameterA/2.0)/(4.0*FocalLength);
				break;
				
			case 'r':
			case 'R': // {Rectangle}				
					radius2 = sqr(Element->ParameterA/2.0) + sqr(Element->ParameterB/2.0);
					Element->ZAperture = radius2/(4.0*FocalLength);
				break;
				
			case 'a':
			case 'A':
			case 'l':
			case 'L':    //annulus or single axis curvature
					if (fabs(Element->ParameterB) >= fabs(Element->ParameterA))
						radius2 = sqr(Element->ParameterB);
					else
						radius2 = sqr(Element->ParameterA);
					Element->ZAperture =radius2/(4.0*FocalLength);
				break;
			
			case 'i':
			case 'I':
			case 'q':
			case 'Q': //irregular triangle or quadrilateral
					radius2 = sqr(MaximumRadius);
					Element->ZAperture =radius2/(4.0*FocalLength);
				break;			
			}
			
		break;
// *************************************************************************                 
	case 'f':
	case 'F': //      {Surface described by flat plane}
			Element->ZAperture = 0.0;
		break;
		
// *************************************************************************                 
	case 'o':
	case 'O': //       {Surface described by other - hyperboloids and ellipsoids}
	
			switch( Element->ShapeIndex )
			{
			case 'c':
			case 'C':
			case 'h':
			case 'H':
			case 't':
			case 'T': // {circle, hexagon, triangle}
					Element->ZAperture = (1-sqrt(1-Element->Kappa*sqr(Element->VertexCurvX)*sqr(Element->ParameterA/2.0)))/(Element->Kappa*Element->VertexCurvX);
				break;
				
			case 'r':
			case 'R': // 			{Rectangle}	
					radius2 = sqr(Element->ParameterA/2.0) + sqr(Element->ParameterB/2.0);
					Element->ZAperture = (1-sqrt(1-Element->Kappa*sqr(Element->VertexCurvX)*radius2))/(Element->Kappa*Element->VertexCurvX);
				break;
			
			case 'a':
			case 'A':
			case 'l':
			case 'L': //annulus or single axis curvature			
					if (fabs(Element->ParameterB) >= fabs(Element->ParameterA) )
						radius2 = sqr(Element->ParameterB);
					else
						radius2 = sqr(Element->ParameterA);
					Element->ZAperture = (1-sqrt(1-Element->Kappa*sqr(Element->VertexCurvX)*radius2))/(Element->Kappa*Element->VertexCurvX);
						 
				break;
			
			case 'i':
			case 'I':
			case 'q':
			case 'Q': //irregular triangle or quadrilateral
					radius2 = sqr(MaximumRadius);
					Element->ZAperture = (1-sqrt(1-Element->Kappa*sqr(Element->VertexCurvX)*radius2))/(Element->Kappa*Element->VertexCurvX);
				break;
			}
			
		break;
// *************************************************************************
	case 'c':
	case 'C': //       {Surface described by cone}	
			switch ( Element->ShapeIndex )
			{
			case 'c':
			case 'C':
			case 'h':
			case 'H':
			case 't':
			case 'T': // {circle, hexagon, triangle}
					Element->ZAperture = Element->ParameterA/(2.0*tan(Element->ConeHalfAngle*(ACOSM1O180)));
				break;
				
			case 'r':
			case 'R': // {Rectangle}
					radius2 = sqr(Element->ParameterA/2.0) + sqr(Element->ParameterB/2.0);
					Element->ZAperture = sqrt(radius2)/(2.0*tan(Element->ConeHalfAngle*(ACOSM1O180)));
				break;
			
			case 'a':
			case 'A': //annulus
					Element->ZAperture = Element->ParameterB/(tan(Element->ConeHalfAngle*(ACOSM1O180)));
				break;
				
			case 'i':
			case 'I':
			case 'q':
			case 'Q': //irregular triangle or quadrilateral				
					radius2 = sqr(MaximumRadius);
					Element->ZAperture = sqrt(radius2)/(2.0*tan(Element->ConeHalfAngle*(ACOSM1O180)));
				break;
			}
			
		break;
		
// *************************************************************************
	case 't':
	case 'T': //       {Surface described by cylinder}
			Element->ZAperture = 2.0/Element->CurvOfRev;
		break;
		
// *************************************************************************
	case 'd':
	case 'D': //        {Surface described by torus}        
			Element->ZAperture = 2.0*Element->CrossSectionRadius;
		break;
		
// *************************************************************************
	case 'm':
	case 'M': //        {Surface described by Zernike monomial}
	
			//Assume that parabola is intended shape.  Use average of B(2,0) and B(2,2) terms
			//for 1/4F and from F and facet shape determine ZAperture
			//FocalLength = (1.0/(4.0*Coefficients[4]) + 1.0/(4.0*Coefficients[6]))/2.0;
			FocalLength = 1.0;

			if (Element->BCoefficients.nrows() > 2
				&& Element->BCoefficients.ncols() > 2)
			{
				FocalLength = (1.0/(4.0 * Element->BCoefficients.at(2,0) ) 
					+ 1.0/(4.0 * Element->BCoefficients.at(2,2)))/2.0;
			}
			else
				return false;
			
			switch( Element->ShapeIndex )
			{
			case 'c':
			case 'C':
			case 'h':
			case 'H':
			case 't':
			case 'T': // 		{circle, hexagon, triangle}
					Element->ZAperture =sqr(Element->ParameterA/2.0)/(4.0*FocalLength);
				break;
			
			case 'r':
			case 'R': // 			{Rectangle}		
					radius2 = sqr(Element->ParameterA/2.0) + sqr(Element->ParameterB/2.0);
					Element->ZAperture = radius2/(4.0*FocalLength);
				break;
				
			case 'a':
			case 'A':
			case 'l':
			case 'L': //annulus or single axis curvature			
					if (fabs(Element->ParameterB) >= fabs(Element->ParameterA) )
						radius2 = sqr(Element->ParameterB);
					else
						radius2 = sqr(Element->ParameterA);
					Element->ZAperture =radius2/(4.0*FocalLength);
				break;
				
			case 'i':
			case 'I':
			case 'q':
			case 'Q': //irregular triangle or quadrilateral
					radius2 = sqr(MaximumRadius);
					Element->ZAperture =radius2/(4.0*FocalLength);
				break;
			}
		break;
		
// *************************************************************************
	case 'r':
	case 'R': //       {Surface described by Polynomial surface of revolution}
			switch( Element->ShapeIndex)
			{
			case 'c':
			case 'C':
			case 'h':
			case 'H':
			case 't':
			case 'T': // 			{circle, hexagon, triangle}
					EvalPoly(Element->ParameterA/2.0, 0, Element->PolyCoeffs, Element->FitOrder, &Element->ZAperture);
				break;
			
			case 'r':
			case 'R': // 		{Rectangle}	
					radius2 = sqr(Element->ParameterA/2.0) + sqr(Element->ParameterB/2.0);
					EvalPoly(sqrt(radius2), 0, Element->PolyCoeffs, Element->FitOrder, &Element->ZAperture);
				break;
			
			case 'a':
			case 'A':
			case 'l':
			case 'L': // annulus or single axis curvature
					if (fabs(Element->ParameterB) >= fabs(Element->ParameterA) )
						radius2 = sqr(Element->ParameterB);
					else
						radius2 = sqr(Element->ParameterA);
						
					EvalPoly(sqrt(radius2), 0, Element->PolyCoeffs, Element->FitOrder, &Element->ZAperture);
				break;
			
			case 'i':
			case 'I':
			case 'q':
			case 'Q': //irregular triangle or quadrilateral
					radius2 = sqr(MaximumRadius);
					EvalPoly(sqrt(radius2), 0, Element->PolyCoeffs, Element->FitOrder, &Element->ZAperture);
				break;
			
			}
		break;
		
// *************************************************************************
	case 'i':
	case 'I': //     {Surface described by cubic spline interpolation data}
			Element->ZAperture = Element->CubicSplineYData[Element->CubicSplineYData.size()-1];  //= to y value of last interpolation data point; can't go beyond this.
		break;
		
// *************************************************************************
	case 'e':
	case 'E': //        {Surface described by finite element data}
                    //assumes that finite element surface has vertex at element origin and that
                    //finite element coordinate coincides with element coordinate system
			lastz = 0.0;
			for (st_uint_t i=0;i<Element->FEData.nrows();i++)
				if (Element->FEData.at(i,2) > lastz) 
					lastz = Element->FEData.at(i,2);
				
			Element->ZAperture = lastz;
		break;
		
// *************************************************************************
	case 'v':
	case 'V':
       //{Surface described by VSHOT data}   {Don't believe this is the correct way to evaluate ZAperture.  It should be the largest Z value from the input dataset, not the terms described below
       //                                    which have no meaning for higher order fits or for single axis curvature sections where the average of B(20) and B(22) is not a good estimate.}
			lastz = 0.0;
			for (st_uint_t i=0;i<Element->VSHOTData.nrows();i++)
			{
				EvalMono(Element->VSHOTData.at(i,0), Element->VSHOTData.at(i,1), 
					Element->BCoefficients, Element->FitOrder, 0.0, 0.0, &zm);

				if (zm > lastz) lastz = zm;
			}
			Element->ZAperture = lastz;
		break;
	}

	return true;
}
//End of Procedure--------------------------------------------------------------
