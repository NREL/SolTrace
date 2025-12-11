
#include "pt_optimizations.hpp"

#include <algorithm>

#include "constants.hpp"
#include "simulation_data_export.hpp"

namespace SolTrace::NativeRunner
{

	void SetupPTOptimizations(
		// system info
		TSystem *System, const bool AsPowerTower,
		// outputs
		st_hash_tree &sun_hash,
		st_hash_tree &rec_hash,
		double (&reccm_helio)[3])
	{
		// Calculate the center of mass of the receiver stage (StageList[1]) in
		// heliostat stage coordinates.
		double reccm[] = {0., 0., 0.};
		int nelrec = 0;
		if (AsPowerTower)
		{
			for (uint_fast64_t j = 0; j < System->StageList[1]->ElementList.size(); j++)
			{
				TElement *el = System->StageList[1]->ElementList.at(j).get();

				// if (!el->Enabled)
				// 	continue;

				nelrec++;

				for (int jj = 0; jj < 3; jj++)
					reccm[jj] += el->Origin[jj];
			}
			for (int jj = 0; jj < 3; jj++)
				reccm[jj] /= (double)nelrec; // average

			// Transform to reference
			double dum1[] = {0., 0., 1.};
			double dum2[3];
			double reccm_global[3];
			TransformToReference(reccm, dum1, System->StageList[1]->Origin,
								 System->StageList[1]->RLocToRef, reccm_global, dum2);

			// Transform to local (heliostat). reccm_helio is the x,y,z position
			// of the receiver centroid in heliostat stage coordinates.
			TransformToLocal(reccm_global, dum1, System->StageList[0]->Origin,
							 System->StageList[0]->RRefToLoc, reccm_helio, dum2);
		}
		// Create an array that stores the element address and the projected size
		// in polar coordinates
		std::vector<eprojdat> el_proj_dat;
		el_proj_dat.reserve(System->StageList[0]->ElementList.size());

		// calculate the smallest zone size. This should be on the order of the
		// largest element in the stage. load stage 0 elements into the mesh
		double d_elm_max = -9.e9;
		double d_elm;

		for (uint_fast64_t i = 0; i < System->StageList[0]->ElementList.size(); i++)
		{
			TElement *el = System->StageList[0]->ElementList.at(i).get();

			// el->element_number = i + 1; // use index for element number

			d_elm = el->aperture->diameter_circumscribed_circle();

			// double d_elm;

			// switch (el->ShapeIndex)
			// {
			// // circular aperture
			// case 'c':
			// case 'C':
			// 	// hexagonal aperture
			// case 'h':
			// case 'H':
			// 	// triangular aperture
			// case 't':
			// case 'T':
			// 	d_elm = el->ParameterA;
			// 	break;
			// 	// rectangular aperture
			// case 'r':
			// case 'R':
			// 	d_elm = sqrt(el->ParameterA * el->ParameterA + el->ParameterB * el->ParameterB);
			// 	break;
			// 	// annular aperture
			// case 'a':
			// case 'A':
			// 	d_elm = el->ParameterB;
			// 	break;
			// case 'l':
			// case 'L':
			// 	// off axis aperture section of line focus trough  or cylinder
			// 	d_elm = sqrt(el->ParameterB * el->ParameterB * 4. + el->ParameterC * el->ParameterC);
			// 	break;
			// 	// Irregular triangle
			// case 'i':
			// case 'I':
			// 	// irregular quadrilateral
			// case 'q':
			// case 'Q':
			// {
			// 	double xmax = fmax(el->ParameterA, fmax(el->ParameterC, el->ParameterE));
			// 	double xmin = fmin(el->ParameterA, fmin(el->ParameterC, el->ParameterE));
			// 	double ymax = fmax(el->ParameterB, fmax(el->ParameterD, el->ParameterF));
			// 	double ymin = fmin(el->ParameterB, fmin(el->ParameterD, el->ParameterF));

			// 	if (el->ShapeIndex == 'q' || el->ShapeIndex == 'Q')
			// 	{
			// 		xmax = fmax(xmax, el->ParameterG);
			// 		xmin = fmin(xmin, el->ParameterG);
			// 		ymax = fmax(ymax, el->ParameterH);
			// 		ymin = fmin(ymin, el->ParameterH);
			// 	}

			// 	double dx = xmax - xmin;
			// 	double dy = ymax - ymin;

			// 	d_elm = sqrt(dx * dx + dy * dy);

			// 	break;
			// }
			// default:
			// 	break;
			// }

			d_elm_max = fmax(d_elm_max, d_elm);

			if (AsPowerTower)
			{
				// Calculate the distance from the receiver to the element and the max projected size
				double dX[3];
				for (int jj = 0; jj < 3; jj++)
					dX[jj] = el->Origin[jj] - reccm_helio[jj]; // vector from receiver to heliostat (not unitized)
				double r_elm = 0.;
				for (int jj = 0; jj < 3; jj++)
					r_elm += dX[jj] * dX[jj];
				r_elm = sqrt(r_elm);			   // vector length
				double d_elm_proj = d_elm / r_elm; // Projected size of the element from the view of the receiver (radians)

				// calculate az,zen coordinate
				double az, zen;
				az = atan2(dX[0] / r_elm, dX[1] / r_elm); // Az coordinate of the heliostat from the receiver's perspective
				zen = asin(dX[2] / r_elm);				  // Zen coordinate """"

				el_proj_dat.push_back(eprojdat(el, d_elm_proj, az, zen));
			}
		}

		if (AsPowerTower)
		{
			// Sort the polar projections by size, largest to smallest
			std::sort(el_proj_dat.begin(), el_proj_dat.end(), eprojdat_compare_refactored);
		}

		// set up the layout data object that provides configuration details for
		// the hash tree
		KDLayoutData sun_ld;
		sun_ld.xlim[0] = System->Sun.MinXSun;
		sun_ld.xlim[1] = System->Sun.MaxXSun;
		sun_ld.ylim[0] = System->Sun.MinYSun;
		sun_ld.ylim[1] = System->Sun.MaxYSun;
		sun_ld.min_unit_dx = d_elm_max;
		sun_ld.min_unit_dy = d_elm_max;

		sun_hash.create_mesh(sun_ld);

		// load stage 0 elements into the mesh
		for (uint_fast64_t i = 0; i < System->StageList[0]->ElementList.size(); i++)
		{
			TElement *el = System->StageList[0]->ElementList.at(i).get();
			sun_hash.add_object((void *)el, el->PosSunCoords[0], el->PosSunCoords[1]);
		}

		// calculate and associate neighbors with each zone
		sun_hash.add_neighborhood_data();

		if (AsPowerTower)
		{
			// Set things up for the polar coordinate tree
			KDLayoutData rec_ld;
			rec_ld.xlim[0] = -PI;
			rec_ld.xlim[1] = PI;
			rec_ld.ylim[0] = -PI / 2.;
			rec_ld.ylim[1] = PI / 2.;
			// use smallest element to set the minimum size
			rec_ld.min_unit_dx = rec_ld.min_unit_dy = el_proj_dat.back().d_proj; // radians at equator

			rec_hash.create_mesh(rec_ld);

			// load stage 0 elements into the receiver mesh in the order of largest projection to smallest
			for (int i = 0; i < el_proj_dat.size(); i++)
			{
				eprojdat *D = &el_proj_dat.at(i);

				// Calculate the angular span of the element
				double angspan[2];
				double adjmult = 1.5;
				angspan[0] = D->d_proj / cos(fabs(D->zen)) * adjmult; // azimuthal span
				angspan[0] = fmin(angspan[0], 2. * PI);				  // limit to circumference
				angspan[1] = D->d_proj / PI * adjmult;				  // zenithal span
				rec_hash.add_object((void *)D->el_addr, D->az, D->zen, angspan);
			}

			// associate neighbors with each zone
			rec_hash.add_neighborhood_data();
		}
	}

	uint_fast64_t GetPTElements(
		// system info
		const bool AsPowerTower,

		// Stage info
		// const TStage *Stage
		const tstage_ptr Stage,
		const int i,

		// Ray info
		const bool in_multi_hit_loop,
		const double (&PosRayStage)[3],
		const double (&reccm_helio)[3],
		st_hash_tree *rec_hash,

		const std::vector<void *> &sunint_elements,

		// Outputs
		std::vector<void *> &reflint_elements,
		bool &has_elements)
	{
		uint_fast64_t nintelements = 0;

		if (i == 0)
		{
			if (in_multi_hit_loop)
			{
				if (AsPowerTower)
				{
					//>=Second time through - checking for first stage multiple element interactions

					// get ray position in receiver polar coordinates
					double raypvec[3];
					for (int jj = 0; jj < 3; jj++)
						raypvec[jj] = PosRayStage[jj] - reccm_helio[jj];
					double raypvecmag = sqrt(raypvec[0] * raypvec[0] + raypvec[1] * raypvec[1] + raypvec[2] * raypvec[2]);
					double raypol[2];
					raypol[0] = atan2(raypvec[0], raypvec[1]);
					raypol[1] = asin(raypvec[2] / raypvecmag);
					// get elements in the vicinity of the ray's polar coordinates
					reflint_elements.clear();
					rec_hash->get_all_data_at_loc(reflint_elements, raypol[0], raypol[1]);
					nintelements = reflint_elements.size();
					has_elements = nintelements > 0;
				}
				else
				{
					nintelements = Stage->ElementList.size();
				}
			}
			else
			{
				// First time through - checking for sun ray intersections
				if (has_elements)
					nintelements = sunint_elements.size();
				else
					nintelements = 0;
			}
		}
		else
			nintelements = Stage->ElementList.size();

		return nintelements;
	}

} // namespace SolTrace::NativeRunner
