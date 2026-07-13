#include <gtest/gtest.h>

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

TEST(RayHistoryResult, Accessors)
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

    RayHistoryResult sr;

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

TEST(RayHistoryResult, OstreamOperator)
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

    RayHistoryResult sr;

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

TEST(RayHistoryResult, IndexOperator)
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

    RayHistoryResult sr;

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

TEST(RayHistoryResult, WriteCSV)
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

    RayHistoryResult sr;

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
