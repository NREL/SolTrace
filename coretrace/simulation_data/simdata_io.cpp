
#include "simdata_io.hpp"

#include <array>
#include <cctype>
#include <cstring>
#include <exception>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>

#include <nlohmann/json.hpp>

#include "basic_sun_position.hpp"
#include "constants.hpp"
#include "json_helpers.hpp"
#include "ray_source.hpp"
#include "simulation_data.hpp"
#include "single_element.hpp"
#include "stage_element.hpp"
#include "sun.hpp"
#include "surface.hpp"
#include "utilities.hpp"
#include "virtual_element.hpp"

namespace SolTrace::Data {

DistributionType char_to_distribution(const char dist_char) {
    switch (dist_char) {
    case ('g'): {
        return DistributionType::GAUSSIAN;
    }
    case ('p'): {
        return DistributionType::PILLBOX;
    }
    case ('f'): {
        return DistributionType::DIFFUSE;
    }
    case ('d'): {
        return DistributionType::USER_DEFINED;
    }
    default: {
        return DistributionType::GAUSSIAN;
    }
    }
}

SunShape char_to_sunshape(const char dist_char) {
    switch (dist_char) {
    case ('g'): {
        return SunShape::GAUSSIAN;
    }
    case ('p'): {
        return SunShape::PILLBOX;
    }
    case ('d'): {
        return SunShape::USER_DEFINED;
    }
    default: {
        return SunShape::GAUSSIAN;
    }
    }
}

InteractionType int_to_interaction(const int interaction_int) {
    switch (interaction_int) {
    case (1): {
        return InteractionType::REFRACTION;
    }
    case (2): {
        return InteractionType::REFLECTION;
    }
    default: {
        return InteractionType::REFLECTION;
    }
    }
}

ApertureType char_to_aperture(const char aperture_char) {
    switch (aperture_char) {
    case ('c'): {
        return ApertureType::CIRCLE;
    }
    case ('h'): {
        return ApertureType::HEXAGON;
    }
    case ('t'): {
        return ApertureType::EQUILATERAL_TRIANGLE;
    }
    case ('r'): {
        return ApertureType::RECTANGLE;
    }
    case ('a'): {
        return ApertureType::ANNULUS;
    }
    case ('l'): {
        return ApertureType::SINGLE_AXIS_CURVATURE_SECTION;
    }
    case ('i'): {
        return ApertureType::IRREGULAR_TRIANGLE;
    }
    case ('q'): {
        return ApertureType::IRREGULAR_QUADRILATERAL;
    }
    default: {
        return ApertureType::APERTURE_UNKNOWN;
    }
    }
}

SurfaceType char_to_surface(const char surface_char) {
    switch (surface_char) {
    case ('s'): return SurfaceType::SPHERE;
    case ('p'): return SurfaceType::PARABOLA;
    case ('o'): return SurfaceType::HYPER;
    case ('g'): return SurfaceType::GENERAL_SPENCER_MURTY;
    case ('f'): return SurfaceType::FLAT;
    case ('c'): return SurfaceType::CONE;
    case ('t'): return SurfaceType::CYLINDER;
    case ('d'): return SurfaceType::TORUS;
    default: return SurfaceType::SURFACE_UNKNOWN;
    }
}

static void read_line(char* buf, int len, FILE* fp) {
    fgets(buf, len, fp);
    int nch = strlen(buf);
    if (nch > 0 && buf[nch - 1] == '\n') buf[nch - 1] = 0;
    if (nch - 1 > 0 && buf[nch - 2] == '\r') buf[nch - 2] = 0;
}

std::vector<std::string> split(const std::string& str,
                               const std::string& delim,
                               bool               ret_empty,
                               bool               ret_delim) {
    std::vector<std::string> list;

    char                   cur_delim[2] = { 0, 0 };
    std::string::size_type m_pos        = 0;
    std::string            token;

    while (m_pos < str.length()) {
        std::string::size_type pos = str.find_first_of(delim, m_pos);
        if (pos == std::string::npos) {
            cur_delim[0] = 0;
            token.assign(str, m_pos, std::string::npos);
            m_pos = str.length();
        } else {
            cur_delim[0]               = str[pos];
            std::string::size_type len = pos - m_pos;
            token.assign(str, m_pos, len);
            m_pos = pos + 1;
        }

        if (token.empty() && !ret_empty) continue;

        list.push_back(token);

        if (ret_delim && cur_delim[0] != 0 && m_pos < str.length())
            list.push_back(std::string(cur_delim));
    }

    return list;
}

bool process_sun(FILE* fp, SimulationData& sd) {
    char buf[1024];

    // Read Sun info
    int    bi = 0, count = 0;
    char   cshape = 'g';
    double Sigma, HalfWidth;
    bool   PointSource;
    double X, Y, Z, Latitude, Day, Hour;
    bool   UseLDHSpec;

    read_line(buf, 1023, fp);
    sscanf(buf,
           "SUN\tPTSRC\t%d\tSHAPE\t%c\tSIGMA\t%lg\tHALFWIDTH\t%lg",
           &bi,
           &cshape,
           &Sigma,
           &HalfWidth);
    PointSource = (bi != 0);
    cshape      = static_cast<char>(std::tolower(static_cast<unsigned char>(cshape)));

    read_line(buf, 1023, fp);

    sscanf(buf,
           "XYZ\t%lg\t%lg\t%lg\tUSELDH\t%d\tLDH\t%lg\t%lg\t%lg",
           &X,
           &Y,
           &Z,
           &bi,
           &Latitude,
           &Day,
           &Hour);
    UseLDHSpec = (bi != 0);

    read_line(buf, 1023, fp);
    sscanf(buf, "USER SHAPE DATA\t%d", &count);
    std::vector<double> angle_vec;
    std::vector<double> intensity_vec;
    if (count > 0) {
        for (int i = 0; i < count; i++) {
            double x, y;
            read_line(buf, 1023, fp);
            sscanf(buf, "%lg\t%lg", &x, &y);
            angle_vec.push_back(x);
            intensity_vec.push_back(y);
        }
    }

    // Make sun
    auto sun = make_ray_source<Sun>();

    // Define sun position
    if (UseLDHSpec) {
        // sun->set_position(Latitude, Day, Hour);
        st_sun_position(Latitude, Day, Hour, &X, &Y, &Z);
    }
    sun->set_position(X, Y, Z);
    // else
    // {
    //  sun->set_position(X, Y, Z);
    // }

    // Define sun shape
    SunShape sun_shape = char_to_sunshape(cshape);
    sun->set_shape(sun_shape, Sigma, HalfWidth, 0.0, angle_vec, intensity_vec);
    // TOD: Buie sun shape not implemented here

    // TODO set point source

    // Attach sun to simulation data
    sd.add_ray_source(sun);
    return true;
}

bool read_optic_surface(FILE*               fp,
                        bool                is_front,
                        OpticalPropertySet& optics,
                        int&                OpticalSurfaceNumber,
                        double&             refraction) {
    if (!fp) return false;
    char buf[1024];
    read_line(buf, 1023, fp);
    std::vector<std::string> parts = split(std::string(buf), "\t", true, false);
    if (parts.size() < 15) {
        printf("too few tokens for optical surface: %zu\n", parts.size());
        printf("\t>> %s\n", buf);
        return false;
    }

    char ErrorDistribution = 'g';
    if (parts[1].length() > 0) ErrorDistribution = parts[1][0];

    int ApertureStopOrGratingType = atoi(parts[2].c_str());
    OpticalSurfaceNumber          = atoi(parts[3].c_str());
    int    DiffractionOrder       = atoi(parts[4].c_str());
    double Reflectivity           = atof(parts[5].c_str());
    double Transmissivity         = atof(parts[6].c_str());
    double RMSSlope               = atof(parts[7].c_str());
    double RMSSpecularity         = atof(parts[8].c_str());
    refraction                    = atof(parts[9].c_str());
    // double RefractionIndexImag = atof(parts[10].c_str());   // This is not
    // used
    double GratingCoeffs[4];
    GratingCoeffs[0] = atof(parts[11].c_str());
    GratingCoeffs[1] = atof(parts[12].c_str());
    GratingCoeffs[2] = atof(parts[13].c_str());
    GratingCoeffs[3] = atof(parts[14].c_str());

    bool    UseReflectivityTable = false;
    int     refl_npoints         = 0;
    double* refl_angles          = 0;
    double* refls                = 0;

    bool    UseTransmissivityTable = false;
    int     trans_npoints          = 0;
    double* trans_angles           = 0;
    double* transs                 = 0;

    if (parts.size() >= 17) {
        UseReflectivityTable = (atoi(parts[15].c_str()) > 0);
        refl_npoints         = atoi(parts[16].c_str());
        if (parts.size() >= 19) {
            UseTransmissivityTable = (atoi(parts[17].c_str()) > 0);
            trans_npoints          = atoi(parts[18].c_str());
        }
    }

    if (UseReflectivityTable) {
        refl_angles = new double[refl_npoints];
        refls       = new double[refl_npoints];

        for (int i = 0; i < refl_npoints; i++) {
            read_line(buf, 1023, fp);
            sscanf(buf, "%lg %lg", &refl_angles[i], &refls[i]);
        }
    }
    if (UseTransmissivityTable) {
        trans_angles = new double[trans_npoints];
        transs       = new double[trans_npoints];

        for (int i = 0; i < trans_npoints; i++) {
            read_line(buf, 1023, fp);
            sscanf(buf, "%lg %lg", &trans_angles[i], &transs[i]);
        }
    }

    // Define optical properties
    DistributionType dist = char_to_distribution(ErrorDistribution);
    OpticalSide      side = is_front ? OpticalSide::Front : OpticalSide::Back;
    optics.set_properties(
        side, dist, Transmissivity, Reflectivity, RMSSlope, RMSSpecularity);

    if (refl_angles != 0) delete[] refl_angles;
    if (refls != 0) delete[] refls;
    if (trans_angles != 0) delete[] trans_angles;
    if (transs != 0) delete[] transs;
    return true;
}

bool process_optics(FILE*                                      fp,
                    std::map<std::string, OpticalPropertySet>& optics_map) {
    char buf[1024];

    // Read number of optics
    int count = 0;
    read_line(buf, 1023, fp);
    sscanf(buf, "OPTICS LIST COUNT\t%d", &count);

    // Define each optics
    for (int i = 0; i < count; i++) {
        // Read optical pair info line
        read_line(buf, 1023, fp);

        if (strncmp(buf, "OPTICAL PAIR", 12) == 0) {
            // int iopt = st_add_optic(cxt, (const char*)(buf + 13));
            std::string optics_name = std::string(buf + 13);
            double      refrac_front, refrac_back;
            int         OpticalSurfaceNumber = 0;

            OpticalPropertySet optics_set(InteractionType::UNKNOWN,
                                          optics_name);
            if(!read_optic_surface(fp, true, optics_set, OpticalSurfaceNumber, refrac_front))
                return false;
            if(!read_optic_surface(fp, false, optics_set, OpticalSurfaceNumber, refrac_back))
                return false;
            optics_set.set_refraction_indices(refrac_front, refrac_back);

            optics_map[optics_name] = optics_set;
        } else
            return false;
    }

    return true;
}

bool read_element(FILE*                                      fp,
                  std::map<std::string, OpticalPropertySet>& optics_map,
                  element_ptr&                               el,
                  SimulationData&                            sd,
                  bool                                       virt) {
    char buf[1024];
    read_line(buf, 1023, fp);

    std::vector<std::string> tok = split(buf, "\t", true, false);
    if (tok.size() < 29) {
        printf("too few tokens for element: %zu\n", tok.size());
        printf("\t>> %s\n", buf);
        return false;
    }

    bool       enabled = atoi(tok[0].c_str()) ? 1 : 0;
    glm::dvec3 xyz     = { atof(tok[1].c_str()),
                           atof(tok[2].c_str()),
                           atof(tok[3].c_str()) };
    glm::dvec3 aim     = { atof(tok[4].c_str()),
                           atof(tok[5].c_str()),
                           atof(tok[6].c_str()) };
    double     zrot    = atof(tok[7].c_str());

    char ShapeIndex = ' ';
    if (tok[8].length() > 0) {
        ShapeIndex = tok[8][0];
    } else {
        printf("no aperture index specified for element\n");
        return false;
    }

    std::vector<double> aperture_params;
    for (int i = 0; i < 8; i++) {
        aperture_params.push_back(atof(tok[i + 9].c_str()));
    }

    char SurfaceIndex = ' ';
    if (tok[17].length() > 0) {
        SurfaceIndex = tok[17][0];
    } else {
        printf("no surface index specified for element\n");
        return false;
    }

    std::vector<double> surface_params;
    for (int i = 0; i < 8; i++) {
        surface_params.push_back(atof(tok[i + 18].c_str()));
    }

    // Skipping surface file for now
    std::string     SurfaceFile = tok[26];
    std::string     optics_name = tok[27].c_str();
    InteractionType interaction = int_to_interaction(atoi(tok[28].c_str()));

    // Create element
    if (virt) el = make_element<VirtualElement>();
    else
        el = make_element<SingleElement>();

    if (!enabled) el->disable();

    // Make aperture
    ApertureType aperture_type = char_to_aperture(ShapeIndex);
    if (aperture_type == ApertureType::APERTURE_UNKNOWN) {
        std::stringstream ss;
        ss << "Aperture character " << ShapeIndex
           << " returned unknown aperture type " << aperture_type;
        throw std::invalid_argument(ss.str());
    }

    aperture_ptr ap_ptr =
        Aperture::make_aperture_from_type(aperture_type, aperture_params);
    if (ap_ptr == nullptr) {
        std::stringstream ss;
        ss << "Unable to make aperture pointer -- "
           << "\nChar: " << ShapeIndex << "\nType: " << aperture_type
           << "\nParams: [";
        for (auto cit = aperture_params.cbegin(); cit != aperture_params.cend();
             ++cit) {
            ss << *cit << ", ";
        }
        ss << "]" << std::endl;
        throw std::runtime_error(ss.str());
    }
    el->set_aperture(ap_ptr);

    // Make surface
    SurfaceType surface_type = char_to_surface(SurfaceIndex);
    if (surface_type == SurfaceType::SURFACE_UNKNOWN) {
        std::stringstream ss;
        ss << "Unknown surface type " << surface_type;
        throw std::invalid_argument(ss.str());
    }
    surface_ptr surf_ptr = make_surface_from_type(surface_type, surface_params);
    el->set_surface(surf_ptr);
    if (surface_type == SurfaceType::CYLINDER) {
        ap_ptr    = el->get_aperture();
        auto rect = std::dynamic_pointer_cast<Rectangle>(ap_ptr);
        auto cyl  = std::dynamic_pointer_cast<Cylinder>(surf_ptr);
        if (rect == nullptr || cyl == nullptr) {
            throw std::invalid_argument("This should not happen!");
        }
        double r = cyl->radius;
        // Force aperture to have correct dimensions
        rect->set_x_length(2.0 * r);
        rect->set_x_coord(-1.0 * r);
        // Adjust for difference between legacy cylinder coordinates
        // and new cylinder coordinates
        glm::dvec3 s_loc = glm::dvec3(0, 0, r);
        glm::dvec3 s_ref;
        el->set_reference_frame_geometry(
            glm::dvec3(xyz), glm::dvec3(aim), zrot);
        el->convert_vector_local_to_reference(s_ref, s_loc);
        xyz = el->get_origin_ref() + s_ref;
        // Because this is actually a point, need to adjust it as well...
        aim = el->get_aim_vector_ref() + s_ref;
    }

    // Set element position and orientation
    el->set_reference_frame_geometry(glm::dvec3(xyz), glm::dvec3(aim), zrot);

    // Set optical properties
    if (!virt) {
        auto optics_iter = optics_map.find(optics_name);
        if (optics_iter == optics_map.end()) {
            std::stringstream ss;
            ss << "Element references unknown optical property set: "
               << optics_name;
            throw std::runtime_error(ss.str());
        }

        OpticalPropertySet optics_set = optics_iter->second;
        optics_set.set_interaction_type(interaction);

        auto optics_ref = sd.find_or_add_optical_property_set(optics_set);
        el->set_optical_property_set(optics_ref);
    }

    return true;
}

bool process_stages(FILE*                                      fp,
                    SimulationData&                            sd,
                    std::map<std::string, OpticalPropertySet>& optics_map) {
    char buf[1024];

    // Loop through stages
    int count_stage = 0;
    read_line(buf, 1023, fp);
    sscanf(buf, "STAGE LIST COUNT\t%d", &count_stage);

    for (int i_stage = 0; i_stage < count_stage; i_stage++) {
        int    virt = 0, multi = 1, count_element = 0, tr = 0;
        double X, Y, Z, AX, AY, AZ, ZRot;

        read_line(buf, 1023, fp);
        sscanf(buf,
               "STAGE\tXYZ\t%lg\t%lg\t%lg\tAIM\t%lg\t%lg\t%lg\tZROT\t%"
               "lg\tVIRTUAL\t%d\tMULTIHIT\t%d\tELEMENTS\t%d\tTRACETHROUGH\t%d",
               &X,
               &Y,
               &Z,
               &AX,
               &AY,
               &AZ,
               &ZRot,
               &virt,
               &multi,
               &count_element,
               &tr);

        read_line(buf, 1023, fp); // read name

        // Make stage
        stage_ptr stage = make_stage(i_stage);
        stage->set_origin(X, Y, Z);
        stage->set_aim_vector(AX, AY, AZ);
        stage->set_zrot(ZRot);
        stage->compute_coordinate_rotations();
        if (virt) {
            stage->mark_virtual();
        } else {
            stage->unmark_virtual();
        }

        // Loop through elements
        for (int i_element = 0; i_element < count_element; i_element++) {
            element_ptr el;
            if (!read_element(fp, optics_map, el, sd, virt) || el == nullptr)
            {
                return false;
            }
            el->set_name(std::to_string(i_element));
            // TODO make virtual if stage is virtual?

            if (!Element::is_success(stage->add_element(el))) {
                std::cout << "Failed to add element to stage" << std::endl;
            }
        }

        if (!Element::is_success(sd.add_stage(stage))) {
            std::cout << "Failed to add stage to SimulationData" << std::endl;
        }
    }

    return true;
}

bool process_sim_par(FILE* fp, SimulationData& sd) {
    char buf[1024];

    // Check if end of file
    if (feof(fp)) return false;

    // Check for simulation parameters
    read_line(buf, 1023, fp);
    if (strncmp(buf, "TRACE", 5) != 0) return false;

    // Get simulation parameters
    int n_rays, n_rays_sun, n_cpu, seed, ss, err, pf;
    sscanf(buf,
           "TRACE\tNRAY\t%d\tNSUN\t%d\tCPU\t%d\tSEED\t%d\tSUNSHAPE\t%"
           "d\tERRORS\t%d\tPTFOCUS\t%d",
           &n_rays,
           &n_rays_sun,
           &n_cpu,
           &seed,
           &ss,
           &err,
           &pf);

    // Assign simulation parameters
    SimulationParameters& par    = sd.get_simulation_parameters();
    par.number_of_rays           = n_rays;
    par.max_number_of_rays       = n_rays_sun;
    par.seed                     = seed;
    par.include_sun_shape_errors = ss;
    par.include_optical_errors   = err;

    // TODO Assign number CPUs, point focus?
    return true;
}

bool load_stinput_file(SimulationData& sd, std::string filename) {
    // Read in file
    FILE* fp = fopen(filename.data(), "r");
    if (!fp) {
        printf("failed to open system input file: %s\n", filename.data());
        return false;
    }

    // Clear simulation data after the file has successfully opened
    sd.clear();

    // Buffer to store read line
    char buf[1024];

    // Get version info (if first line starts with '#')
    int  vmaj = 0, vmin = 0, vmic = 0;
    char c = fgetc(fp);
    if (c == '#') {
        read_line(buf, 1023, fp);
        sscanf(
            buf, " SOLTRACE VERSION %d.%d.%d INPUT FILE", &vmaj, &vmin, &vmic);

        // unsigned int file_version = vmaj*10000 + vmin*100 + vmic;

        printf("loading input file version %d.%d.%d\n", vmaj, vmin, vmic);
    } else {
        ungetc(c, fp);
    }

    // Read in Sun
    if(!process_sun(fp, sd))
        return false;

    // Read in Optics
    std::map<std::string, OpticalPropertySet> optics_map;
    if (!process_optics(fp, optics_map))
        return false;

    // Read in Stages
    if (!process_stages(fp, sd, optics_map))
        return false;

    // Read in simulation parameters (if any)
    process_sim_par(fp, sd);

    fclose(fp);

    return true;
}

} // namespace SolTrace::Data
