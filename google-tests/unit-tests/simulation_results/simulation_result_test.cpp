#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <sstream>

#include <simulation_result.hpp>

#include "common.hpp"

TEST(InteractionRecord, Constructors)
{
    element_id elid = 5;
    RayEvent it = RayEvent::ABSORB;
    glm::dvec3 loc(-3.2, -2.5, 1.2);
    glm::dvec3 dir(-1.0, 1.5, -2.5);
    InteractionRecord ir1(elid, it, loc, dir);
    InteractionRecord ir2(elid, it,
                          loc[0], loc[1], loc[2],
                          dir[0], dir[1], dir[2]);

    // EXPECT_TRUE(ir1.index < 0);
    // EXPECT_TRUE(ir2.index < 0);
    EXPECT_EQ(ir1.element, ir2.element);
    EXPECT_EQ(ir1.event, ir2.event);
    EXPECT_EQ(ir1.location[0], ir2.location[0]);
    EXPECT_EQ(ir1.location[1], ir2.location[1]);
    EXPECT_EQ(ir1.location[2], ir2.location[2]);
}

TEST(InteractionRecord, OutputOperator)
{
    element_id elid = 5;
    RayEvent it = RayEvent::ABSORB;
    glm::dvec3 loc(-3.2, -2.5, 1.2);
    glm::dvec3 dir(-1.0, 1.5, -2.5);
    InteractionRecord ir1(elid, it, loc, dir);

    // Basic test to see that the output operator is there
    // and functioning.
    std::stringstream ss;
    ss << ir1;
    EXPECT_TRUE(ss.str().size() > 0);
}

TEST(InteractionRecord, RayEventToString)
{
    EXPECT_EQ(SolTrace::Result::ray_event_string(RayEvent::CREATE), "CREATE");
    EXPECT_EQ(SolTrace::Result::ray_event_string(RayEvent::ABSORB), "ABSORB");
    EXPECT_EQ(SolTrace::Result::ray_event_string(RayEvent::REFLECT), "REFLECT");
    EXPECT_EQ(SolTrace::Result::ray_event_string(RayEvent::TRANSMIT), "TRANSMIT");
    EXPECT_EQ(SolTrace::Result::ray_event_string(RayEvent::VIRTUAL), "VIRTUAL");
    EXPECT_EQ(SolTrace::Result::ray_event_string(RayEvent::EXIT), "EXIT");
    EXPECT_EQ(SolTrace::Result::ray_event_string(RayEvent::UNKNOWN), "UNKNOWN");
}

TEST(RayRecord, Accessors)
{
    // Test constants
    const double TOL = 1e-12;
    const uint_fast32_t NINTER = 5;
    const int_fast64_t ID = 5;
    RayEvent my_types[NINTER] = {
        RayEvent::CREATE,
        RayEvent::TRANSMIT,
        RayEvent::REFLECT,
        RayEvent::ABSORB,
        RayEvent::EXIT};

    // Adding records and sizing
    ray_record_ptr rr = make_ray_record(ID);
    for (uint_fast32_t ell = 0; ell < NINTER; ++ell)
    {
        glm::dvec3 loc(1.0 * ell * ell, 2.0 * ell * ell, 3.0 * ell * ell);
        glm::dvec3 dir(2.0 * ell + 1.0, 4.0 * ell + 2.0, 6.0 * ell + 3.0);
        RayEvent it = my_types[ell];
        interaction_ptr ir = make_interaction_record(ell, it, loc, dir);
        rr->add_interaction_record(ir);
        EXPECT_EQ(rr->get_number_of_interactions(), ell + 1);
    }

    // Getting direction cosines and indexing
    int_fast64_t idx = 0;
    glm::dvec3 ret1;
    glm::dvec3 ret2;
    glm::dvec3 lcos;
    for (auto iter = rr->interactions.cbegin();
         iter != rr->interactions.cend();
         ++iter)
    {
        EXPECT_EQ(*iter, (*rr)[idx]);

        rr->get_direction(idx, ret1);
        rr->get_direction(*iter, ret2);
        EXPECT_TRUE(is_identical(ret1, ret2));

        lcos = {1.0, 2.0, 3.0};
        lcos *= 2 * idx + 1;

        SolTrace::Data::normalize_inplace(lcos);
        EXPECT_NEAR(ret1[0], lcos[0], TOL);
        EXPECT_NEAR(ret1[1], lcos[1], TOL);
        EXPECT_NEAR(ret1[2], lcos[2], TOL);

        ret1 = glm::dvec3(0.0);
        rr->get_position(idx, ret1);
        rr->get_position(*iter, ret2);
        EXPECT_TRUE(is_identical(ret1, ret2));

        EXPECT_EQ(rr->get_element(idx), idx);
        EXPECT_EQ(rr->get_element(*iter), idx);

        EXPECT_EQ(rr->get_event(idx), my_types[idx]);
        EXPECT_EQ(rr->get_event(*iter), my_types[idx]);

        ++idx;
    }
}

TEST(RayRecord, OutputOperator)
{
    // Test constants
    const uint_fast32_t NINTER = 5;
    const int_fast64_t ID = 5;
    RayEvent my_types[NINTER] = {
        RayEvent::CREATE,
        RayEvent::TRANSMIT,
        RayEvent::REFLECT,
        RayEvent::ABSORB,
        RayEvent::EXIT};

    // Adding records and sizing
    RayRecord rr(ID);
    for (uint_fast32_t ell = 0; ell < NINTER; ++ell)
    {
        glm::dvec3 loc(1.0 * ell * ell, 2.0 * ell * ell, 3.0 * ell * ell);
        glm::dvec3 dir(2.0 * ell + 1.0, 4.0 * ell + 2.0, 6.0 * ell + 3.0);
        RayEvent it = my_types[ell];
        interaction_ptr ir = make_interaction_record(ell, it, loc, dir);
        rr.add_interaction_record(ir);
    }

    std::stringstream ss;
    ss << rr;
    EXPECT_TRUE(ss.str().size() > 0);
}

TEST(ElementRecord, OutputOperator)
{
    // Test constants
    const uint_fast32_t NINTER = 5;
    const int_fast64_t ID = 5;
    RayEvent my_types[NINTER] = {
        RayEvent::CREATE,
        RayEvent::TRANSMIT,
        RayEvent::REFLECT,
        RayEvent::ABSORB,
        RayEvent::EXIT};

    // Adding records and sizing
    ElementRecord erec(ID);
    for (uint_fast32_t ell = 0; ell < NINTER; ++ell)
    {
        glm::dvec3 loc(1.0 * ell * ell, 2.0 * ell * ell, 3.0 * ell * ell);
        glm::dvec3 dir(2.0 * ell + 1.0, 4.0 * ell + 2.0, 6.0 * ell + 3.0);
        RayEvent it = my_types[ell];
        interaction_ptr ir = make_interaction_record(ell, it, loc, dir);
        erec.add_interaction_record(ir);
    }

    std::stringstream ss;
    ss << erec;
    EXPECT_TRUE(ss.str().size() > 0);

}

TEST(SimulationResult, Accessors)
{
    // Test constants
    const uint_fast32_t NRAYS = 3;
    const uint_fast32_t NINTER = 5;
    RayEvent my_types[NINTER] = {
        RayEvent::REFLECT,
        RayEvent::TRANSMIT,
        RayEvent::REFLECT,
        RayEvent::ABSORB,
        RayEvent::EXIT};

    SimulationResult sr;

    for (uint_fast32_t k = 0; k < NRAYS; ++k)
    {
        ray_record_ptr rr = make_ray_record(k);
        for (uint_fast32_t ell = 0; ell < NINTER; ++ell)
        {
            glm::dvec3 loc(1.0 * ell, 2.0 * ell, 3.0 * ell);
            glm::dvec3 dir(1.0, 2.0, 3.0);
            RayEvent it = my_types[ell];
            interaction_ptr ir = make_interaction_record(ell, it, loc, dir);
            rr->add_interaction_record(ir);
        }
        sr.add_ray_record(rr);
        EXPECT_EQ(sr.get_number_of_records(), k + 1);
        EXPECT_EQ(sr[k], rr);
    }

    int_fast64_t idx = 0;
    for (auto iter = sr.get_ray_record_iterator(); !sr.is_at_end(iter); ++iter)
    {
        EXPECT_EQ(*iter, sr[idx]);
        ++idx;
    }
    EXPECT_EQ(sr.get_number_of_records(), idx);
}

TEST(SimulationResult, SunSamplingStats)
{
    const double TOL = 1e-12;

    // --- Three-argument overload (width, height, sun_ray_count) ---
    {
        SimulationResult sr;
        sr.set_sun_sampling_stats(4.0, 3.0, 12);

        double w, h;
        sr.get_sun_dimensions(w, h);
        EXPECT_NEAR(w, 4.0, TOL);
        EXPECT_NEAR(h, 3.0, TOL);
        EXPECT_NEAR(sr.get_sun_A_box(), 12.0, TOL);
        EXPECT_EQ(sr.get_sun_ray_count(), 12u);
        EXPECT_NEAR(sr.get_ray_area_weight(), 12.0 / 12.0, TOL);
    }

    // --- Two-argument overload (sun_A_box, sun_ray_count) ---
    {
        SimulationResult sr;
        sr.set_sun_sampling_stats(25.0, 5);

        double w, h;
        sr.get_sun_dimensions(w, h);
        EXPECT_TRUE(std::isnan(w));
        EXPECT_TRUE(std::isnan(h));
        EXPECT_NEAR(sr.get_sun_A_box(), 25.0, TOL);
        EXPECT_EQ(sr.get_sun_ray_count(), 5u);
        EXPECT_NEAR(sr.get_ray_area_weight(), 25.0 / 5.0, TOL);
    }

    // --- Three-argument overload: invalid inputs ---
    {
        SimulationResult sr;
        EXPECT_THROW(sr.set_sun_sampling_stats(4.0, 3.0, 0),  std::invalid_argument); // zero rays
        EXPECT_THROW(sr.set_sun_sampling_stats(0.0, 3.0, 10), std::invalid_argument); // zero width
        EXPECT_THROW(sr.set_sun_sampling_stats(-1.0, 3.0, 10), std::invalid_argument); // negative width
        EXPECT_THROW(sr.set_sun_sampling_stats(4.0, 0.0, 10), std::invalid_argument); // zero height
        EXPECT_THROW(sr.set_sun_sampling_stats(4.0, -1.0, 10), std::invalid_argument); // negative height
    }

    // --- Two-argument overload: invalid inputs ---
    {
        SimulationResult sr;
        EXPECT_THROW(sr.set_sun_sampling_stats(25.0, 0),   std::invalid_argument); // zero rays
        EXPECT_THROW(sr.set_sun_sampling_stats(0.0, 10),   std::invalid_argument); // zero area
        EXPECT_THROW(sr.set_sun_sampling_stats(-1.0, 10),  std::invalid_argument); // negative area
    }
}

TEST(SimulationResult, OstreamOperator)
{

    // Test constants
    const uint_fast32_t NRAYS = 3;
    const uint_fast32_t NINTER = 5;
    RayEvent my_types[NINTER] = {
        RayEvent::REFLECT,
        RayEvent::TRANSMIT,
        RayEvent::REFLECT,
        RayEvent::ABSORB,
        RayEvent::EXIT};

    SimulationResult sr;

    for (uint_fast32_t k = 0; k < NRAYS; ++k)
    {
        ray_record_ptr rr = make_ray_record(k);
        for (uint_fast32_t ell = 0; ell < NINTER; ++ell)
        {
            glm::dvec3 loc(1.0 * ell, 2.0 * ell, 3.0 * ell);
            glm::dvec3 dir(1.0, 2.0, 3.0);
            RayEvent it = my_types[ell];
            interaction_ptr ir = make_interaction_record(ell, it, loc, dir);
            rr->add_interaction_record(ir);
        }
        sr.add_ray_record(rr);
    }

    std::stringstream ss;
    ss << sr;

    EXPECT_TRUE(ss.str().size() > 0);
}

TEST(SimulationResult, IndexOperator)
{
    // Test constants
    const uint_fast32_t NRAYS = 3;
    const uint_fast32_t NINTER = 5;
    const int_fast32_t INDEX = 0;
    RayEvent my_types[NINTER] = {
        RayEvent::REFLECT,
        RayEvent::TRANSMIT,
        RayEvent::REFLECT,
        RayEvent::ABSORB,
        RayEvent::EXIT};

    SimulationResult sr;

    for (uint_fast32_t k = 0; k < NRAYS; ++k)
    {
        ray_record_ptr rr = make_ray_record(k);
        for (uint_fast32_t ell = 0; ell < NINTER; ++ell)
        {
            glm::dvec3 loc(1.0 * ell, 2.0 * ell, 3.0 * ell);
            glm::dvec3 dir(1.0, 2.0, 3.0);
            RayEvent it = my_types[ell];
            interaction_ptr ir = make_interaction_record(ell, it, loc, dir);
            rr->add_interaction_record(ir);
        }
        sr.add_ray_record(rr);
    }

    EXPECT_THROW(sr[-1], std::invalid_argument);
    EXPECT_THROW(sr[NRAYS + 1], std::invalid_argument);

    ray_record_ptr rr = sr[INDEX];
    EXPECT_NE(rr, nullptr);
    EXPECT_THROW((*rr)[-1], std::invalid_argument);
    EXPECT_THROW((*rr)[NINTER + 1], std::invalid_argument);

    interaction_ptr ir = (*rr)[INDEX];
    EXPECT_EQ(ir->event, my_types[INDEX]);
}

TEST(SimulationResult, WriteCSV)
{
    // Test constants
    const uint_fast32_t NRAYS = 3;
    const uint_fast32_t NINTER = 5;
    RayEvent my_types[NINTER] = {
        RayEvent::REFLECT,
        RayEvent::TRANSMIT,
        RayEvent::REFLECT,
        RayEvent::ABSORB,
        RayEvent::EXIT};

    std::string csv_file("temp_string.csv");

    SimulationResult sr;

    for (uint_fast32_t k = 0; k < NRAYS; ++k)
    {
        ray_record_ptr rr = make_ray_record(k);
        for (uint_fast32_t ell = 0; ell < NINTER; ++ell)
        {
            glm::dvec3 loc(1.0 * ell, 2.0 * ell, 3.0 * ell);
            glm::dvec3 dir(1.0, 2.0, 3.0);
            RayEvent it = my_types[ell];
            interaction_ptr ir = make_interaction_record(ell, it, loc, dir);
            rr->add_interaction_record(ir);
        }
        sr.add_ray_record(rr);
    }

    // Write the CSV files
    sr.write_csv_file(csv_file);
    sr.write_csv_file("temp_char.csv");

    // Cleanup
    std::filesystem::remove(csv_file);
    std::filesystem::remove("temp_char.csv");
}
