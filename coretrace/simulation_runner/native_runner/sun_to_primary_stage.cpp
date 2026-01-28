
#include "sun_to_primary_stage.hpp"

#include "native_runner_types.hpp"
#include "simulation_data_export.hpp"
#include "trace_logger.hpp"

namespace SolTrace::NativeRunner
{

    bool SunToPrimaryStage(
        trace_logger_ptr logger,
        TSystem *System,
        TStage *Stage,
        TSun *Sun,
        double PosSunStage[3])
    {

        /*{Purpose: To compute the sun position within primary sage and the maximum radius of a cicle seen from sun which encircles
        all elements  within the primary stage.  Used for genenating rays from sun.   modified on 09/26/05 to establish rectangular
        region of interest - more efficient.
             Input - Sys     = Primary stage
                     Sun       = Sun description block
             Output - PosSunStage  = position of sun in primary stage coordinate system
                      Sun.MaxRad =    maximum radius of a cicle seen from sun which encircles
                                  all elements within the primary stage  relative to center of mass of all elements in that stage
                      Sun.Xcm
                      Sun.Ycm    =  center of mass of all elements in primary stage as seen in sun coordinate system
        }*/

        double dx = 0, dy = 0, dz = 0, dtot = 0;
        double CosSunGlob[3] = {0.0, 0.0, 0.0};
        double PosSunGlob[3] = {0.0, 0.0, 0.0};
        double CosSunStage[3] = {0.0, 0.0, 0.0};

        uint_fast64_t i = 0;
        double x = 0, y = 0, radius = 0;
        double Origin[3] = {0.0, 0.0, 0.0};
        double CosDum[3] = {0.0, 0.0, 0.0};
        double PosLoc[3] = {0.0, 0.0, 0.0};
        double CosLoc[3] = {0.0, 0.0, 0.0};
        double RRefToLoc[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        double radius1 = 0.0, radius2 = 0.0, radius3 = 0.0, radius4 = 0.0, radiustemp = 0.0;
        double Xsum = 0.0, Ysum = 0.0, xminsun = 0.0, yminsun = 0.0, xmaxsun = 0.0, ymaxsun = 0.0;
        double XLegofRadius = 0.0;

        double element_radius;

        // PosSunGlob[0] = Sun.Origin[0];//Position of sun coord. system origin in global system
        // PosSunGlob[1] = Sun.Origin[1]; //changed 5/1/00 to place sun at primary stage origin; direction vector
        // PosSunGlob[2] = Sun.Origin[2]; //calculated below from difference between entered sun position and global
        PosSunGlob[0] = Stage->Origin[0]; // origin
        PosSunGlob[1] = Stage->Origin[1];
        PosSunGlob[2] = Stage->Origin[2];

        // First calculate direction cosines of sun z-axis in global coord. system
        dx = 0.0 - Sun->Origin[0]; // changed 5/1/00 to tie the sun direction to global coordinate system origin
        dy = 0.0 - Sun->Origin[1]; // for any stage; not different for each stage.  this represents reality because
        dz = 0.0 - Sun->Origin[2]; // sun is essentially inifinitely far away.
        dtot = sqrt(dx * dx + dy * dy + dz * dz);

        if (dtot == 0.0)
        {
            logger->error_log("error calculating sun position in primary stage, dtot = 0.0\n");
            return false;
        }

        dx = dx / dtot; // unit vector in global coord.system of sun coord. system z axis.
        dy = dy / dtot;
        dz = dz / dtot;

        CosSunGlob[0] = dx; // direction cosines of sun Z-axis in global system.
        CosSunGlob[1] = dy;
        CosSunGlob[2] = dz;

        // Transform sun direction vector to Stage system; CosSunStage is dir cosines of sun ray in Stage coord. system
        // PosSunStage is position of sun coord. system origin in Stage system
        TransformToLocal(PosSunGlob, CosSunGlob, Stage->Origin, Stage->RRefToLoc, PosSunStage, CosSunStage);

        Sun->Euler[0] = atan2(CosSunStage[0], CosSunStage[2]); // Euler angles relating sun to Stage system
        Sun->Euler[1] = asin(CosSunStage[1]);
        Sun->Euler[2] = 0.0;

        /*     {Now we have the Euler angles from Stage to the sun coordinate system.  We have to now transform the
              element locations in the stage system to the sun coordinate system and find the smallest circle in the
              xy plane of the sun system that completely encompasses the projected images of the elements onto that plane}*/

        Origin[0] = 0.0; // Origin of transformed system and stage system the same
        Origin[1] = 0.0;
        Origin[2] = 0.0;

        CosDum[0] = 0.0; // direction cosines not important; only interested in point locations
        CosDum[1] = 0.0;
        CosDum[2] = 1.0;

        Sun->MaxRad = 0.0;
        Sun->Xcm = 0.0;
        Sun->Ycm = 0.0;
        Sun->MaxXSun = -1.0e20;
        Sun->MinXSun = 1.0e20;
        Sun->MaxYSun = -1.0e20;
        Sun->MinYSun = 1.0e20;

        CalculateTransformMatrices(Sun->Euler, RRefToLoc, Sun->RLocToRef);

        //{Now calculate center of mass of projected distribution. Added 09/26/05}
        Xsum = 0.0;
        Ysum = 0.0;
        // for (i = 0; i < Stage->ElementList.size(); i++)
        // {
        // 	// if (!Stage->ElementList[i]->Enabled)
        // 	// 	continue;
        // 	TransformToLocal(Stage->ElementList[i]->Origin, CosDum, Origin, RRefToLoc, PosLoc, CosLoc);
        // 	// Now have PosLoc which is the projected position of element[i] in xy plane of sun coord. system
        // 	Xsum = Xsum + PosLoc[0];
        // 	Ysum = Ysum + PosLoc[1];
        // }

        telement_ptr elem;
        for (auto iter = Stage->ElementList.cbegin();
             iter != Stage->ElementList.cend();
             ++iter)
        {
            // TransformToLocal(iter->second->Origin, CosDum, Origin,
            // 				 RRefToLoc, PosLoc, CosLoc);
            elem = *iter;
            TransformToLocal(elem->Origin, CosDum, Origin,
                             RRefToLoc, PosLoc, CosLoc);
            Xsum += PosLoc[0];
            Ysum += PosLoc[1];
        }

        // center of mass of distribution of element locations as projected
        // in sun coord system.   Added 09/26/05
        Sun->Xcm = Xsum / Stage->ElementList.size();
        Sun->Ycm = Ysum / Stage->ElementList.size();

        // std::cout << "Xcm = " << Sun->Xcm
        // 		  << "\nYcm = " << Sun->Ycm
        // 		  << "\nnelement = " << Stage->ElementList.size()
        // 		  << std::endl;

        size_t nelements = 0;
        elem = nullptr;

        for (auto iter = Stage->ElementList.cbegin();
             iter != Stage->ElementList.cend();
             ++iter)
        {
            elem = *iter;

            // TransformToLocal(Stage->ElementList[i]->Origin, CosDum, Origin, RRefToLoc, PosLoc, CosLoc);
            TransformToLocal(elem->Origin, CosDum, Origin,
                             RRefToLoc, PosLoc, CosLoc);
            // Now have PosLoc which is the projected position of element[i] in xy plane of sun coord. system
            x = PosLoc[0] - Sun->Xcm; // changes origin to center of mass of all elements  09/26/05
            y = PosLoc[1] - Sun->Ycm;
            radius = sqrt(x * x + y * y);

            xminsun = PosLoc[0];
            xmaxsun = PosLoc[0];
            yminsun = PosLoc[1];
            ymaxsun = PosLoc[1];

            // save the projected position of the element on the sun coordinate plane
            elem->PosSunCoords[0] = PosLoc[0];
            elem->PosSunCoords[1] = PosLoc[1];
            elem->PosSunCoords[2] = PosLoc[2];

            // element_radius = 0.5 * Stage->ElementList[i]->aperture->diameter_circumscribed_circle();
            element_radius = elem->aperture->radius_circumscribed_circle();
            radius += element_radius;
            xminsun -= element_radius;
            yminsun -= element_radius;
            xmaxsun += element_radius;
            ymaxsun += element_radius;

            if (radius > Sun->MaxRad)
                Sun->MaxRad = radius; // establishes a circular region

            // restablishes a rectangular region instead of a circular region
            // Added 09/26/05
            if (xminsun < Sun->MinXSun)
                Sun->MinXSun = xminsun;
            if (xmaxsun > Sun->MaxXSun)
                Sun->MaxXSun = xmaxsun;
            if (yminsun < Sun->MinYSun)
                Sun->MinYSun = yminsun;
            if (ymaxsun > Sun->MaxYSun)
                Sun->MaxYSun = ymaxsun;

            nelements++;
        }

        if (nelements == 0)
            logger->error_log("error calculating sun position in primary stage because there are no elements");

        return (nelements > 0);
    }

} // namespace SolTrace::NativeRunner
