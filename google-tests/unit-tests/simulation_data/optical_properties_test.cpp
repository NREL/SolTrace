#include <gtest/gtest.h>

#include <composite_element.hpp>
#include <limits>
#include <single_element.hpp>
#include <sun.hpp>
#include <simulation_data.hpp>

#include "common.hpp"

TEST(OpticalProperties, DangleOpticalPointer)
{
    SimulationData sd;

    // Make optical property set
    auto opt_set = OpticalPropertySet();
    opt_set.set_ideal_reflection(OpticalSide::Both);

    // Add to simulation data
    auto opt_ref = sd.add_optical_property_set(opt_set);

    // Make element
    auto test_element = make_element<SingleElement>();
    test_element->set_aperture(make_aperture<Circle>(1));
    test_element->set_surface(make_surface<Flat>());

    // Attach optical property reference to element
    test_element->set_optical_property_set(opt_ref);

    // Add element to simulation data
    ASSERT_NO_THROW(sd.add_element(test_element));

    sd.clear();

    EXPECT_EQ(sd.get_optical_property_set(*test_element), nullptr);
}

TEST(OpticalProperties, MissingOptical)
{
    SimulationData sd;

    // Make element
    auto test_element = make_element<SingleElement>();
    test_element->set_aperture(make_aperture<Circle>(1));
    test_element->set_surface(make_surface<Flat>());

    // Add element to simulation data
    EXPECT_THROW(sd.add_element(test_element), std::invalid_argument);
}

TEST(OpticalProperties, InvalidErrorSettingsThrow)
{
    OpticalPropertySet optics(InteractionType::REFLECTION, "invalid_errors");
    optics.set_ideal_reflection(OpticalSide::Both);

    EXPECT_THROW(
        optics.set_errors(OpticalSide::Both, DistributionType::UNKNOWN, 1.0, 1.0),
        std::invalid_argument);
    EXPECT_THROW(
        optics.set_errors(OpticalSide::Both, DistributionType::GAUSSIAN, -1.0, 1.0),
        std::invalid_argument);
    EXPECT_THROW(optics.set_errors(OpticalSide::Both,
                                   DistributionType::GAUSSIAN,
                                   1.0,
                                   std::numeric_limits<double>::quiet_NaN()),
                 std::invalid_argument);
}

TEST(OpticalProperties, MutateOptics)
{
    SimulationData sd;

    // Make optical property set
    auto opt_set = OpticalPropertySet();
    opt_set.set_ideal_reflection(OpticalSide::Both);

    // Add to simulation data
    auto opt_ref = sd.add_optical_property_set(opt_set);

    // Make element 1
    auto test_element_1 = make_element<SingleElement>();
    test_element_1->set_aperture(make_aperture<Circle>(1));
    test_element_1->set_surface(make_surface<Flat>());

    // Attach optical property reference to element 1
    test_element_1->set_optical_property_set(opt_ref);

    // Make element 2
    auto test_element_2 = make_element<SingleElement>();
    test_element_2->set_aperture(make_aperture<Circle>(1));
    test_element_2->set_surface(make_surface<Flat>());

    // Attach same optical properties to elemetn 2
    test_element_2->set_optical_property_set(opt_ref);

    // Modify optical properties directly from simulation data
    auto* optics_ptr =  sd.get_mutable_optical_property_set(*test_element_1);
    optics_ptr->set_reflectivity(OpticalSide::Front, 0.5);

    // Confirm reflectivity is modified for both
    EXPECT_EQ(test_element_1->get_optical_property_set()->get_reflectivity(OpticalSide::Front),
        0.5);
    EXPECT_EQ(test_element_2->get_optical_property_set()->get_reflectivity(OpticalSide::Front),
        0.5);

    EXPECT_EQ(test_element_1->get_optical_property_set(),
        test_element_2->get_optical_property_set());
}

TEST(OpticalProperties, NullPointer)
{
    // Make optics
    OpticalPropertySet opt_set = OpticalPropertySet();
    opt_set.set_ideal_reflection(OpticalSide::Both);
    auto opt_ptr = std::make_shared<OpticalPropertySet>(opt_set);

    // Make optics reference
    OpticalPropertySetReference ref = { 0, opt_ptr };

    // Make element
    auto test_element = make_element<SingleElement>();
    test_element->set_aperture(make_aperture<Circle>(1));
    test_element->set_surface(make_surface<Flat>());
    test_element->set_optical_property_set(ref);

    // Add element to simulation data
    SimulationData sd;
    sd.add_element(test_element);

    // Delete optics pointer
    opt_ptr.reset();

    // Try to get pointer
    EXPECT_EQ(test_element->get_optical_property_set(), nullptr);

}

TEST(OpticalProperties, ElementCanAccessSimulationOwnedOptics)
{
    SimulationData sd;

    OpticalPropertySet opt_set;
    opt_set.set_ideal_reflection(OpticalSide::Both);

    auto opt_ref = sd.add_optical_property_set(opt_set);

    auto test_element = make_element<SingleElement>();
    test_element->set_aperture(make_aperture<Circle>(1));
    test_element->set_surface(make_surface<Flat>());
    test_element->set_optical_property_set(opt_ref);

    EXPECT_NE(test_element->get_optical_property_set(), nullptr);
    EXPECT_EQ(test_element->get_optical_property_set()->get_reflectivity(OpticalSide::Front), 1.0);
}

TEST(OpticalProperties, RemovingElementDoesNotRemoveOptics)
{
    SimulationData sd;

    OpticalPropertySet opt_set;
    opt_set.set_ideal_reflection(OpticalSide::Both);

    auto opt_ref = sd.add_optical_property_set(opt_set);

    auto test_element = make_element<SingleElement>();
    test_element->set_aperture(make_aperture<Circle>(1));
    test_element->set_surface(make_surface<Flat>());
    test_element->set_optical_property_set(opt_ref);

    auto id = sd.add_element(test_element);
    ASSERT_TRUE(SolTrace::Data::Element::is_success(id));

    sd.remove_element(id);

    EXPECT_NE(sd.get_optical_property_set(*test_element), nullptr);
}

TEST(OpticalProperties, FindOrAddDeduplicatesEquivalentOptics)
{
    SimulationData sd;

    OpticalPropertySet opt_set_1;
    opt_set_1.set_ideal_reflection(OpticalSide::Both);

    OpticalPropertySet opt_set_2;
    opt_set_2.set_ideal_reflection(OpticalSide::Both);

    auto ref_1 = sd.find_or_add_optical_property_set(opt_set_1);
    auto ref_2 = sd.find_or_add_optical_property_set(opt_set_2);

    EXPECT_EQ(ref_1.id, ref_2.id);
    EXPECT_EQ(ref_1.optical_property_set.lock(), ref_2.optical_property_set.lock());
}

TEST(OpticalProperties, FindOrAddKeepsDifferentOpticsSeparate)
{
    SimulationData sd;

    OpticalPropertySet opt_set_1;
    opt_set_1.set_ideal_reflection(OpticalSide::Both);

    OpticalPropertySet opt_set_2;
    opt_set_2.set_ideal_reflection(OpticalSide::Both);
    opt_set_2.set_reflectivity(OpticalSide::Front, 0.5);

    auto ref_1 = sd.find_or_add_optical_property_set(opt_set_1);
    auto ref_2 = sd.find_or_add_optical_property_set(opt_set_2);

    EXPECT_NE(ref_1.id, ref_2.id);
    EXPECT_NE(ref_1.optical_property_set.lock(), ref_2.optical_property_set.lock());
}

TEST(OpticalProperties, AddElementWithExpiredOpticsThrows)
{
    auto test_element = make_element<SingleElement>();
    test_element->set_aperture(make_aperture<Circle>(1));
    test_element->set_surface(make_surface<Flat>());

    {
        SimulationData temporary_sd;

        OpticalPropertySet opt_set;
        opt_set.set_ideal_reflection(OpticalSide::Both);

        auto opt_ref = temporary_sd.add_optical_property_set(opt_set);
        test_element->set_optical_property_set(opt_ref);
    }

    SimulationData sd;

    EXPECT_THROW(sd.add_element(test_element), std::invalid_argument);
}

TEST(OpticalProperties, ReplaceElementDoesNotInvalidateOptics)
{
    SimulationData sd;

    OpticalPropertySet opt_set;
    opt_set.set_ideal_reflection(OpticalSide::Both);

    auto opt_ref = sd.add_optical_property_set(opt_set);

    auto old_element = make_element<SingleElement>();
    old_element->set_aperture(make_aperture<Circle>(1));
    old_element->set_surface(make_surface<Flat>());
    old_element->set_optical_property_set(opt_ref);

    auto id = sd.add_element(old_element);
    ASSERT_TRUE(SolTrace::Data::Element::is_success(id));

    auto new_element = make_element<SingleElement>();
    new_element->set_aperture(make_aperture<Circle>(2));
    new_element->set_surface(make_surface<Flat>());
    new_element->set_optical_property_set(opt_ref);

    EXPECT_TRUE(sd.replace_element(id, new_element));

    EXPECT_NE(new_element->get_optical_property_set(), nullptr);
    EXPECT_NE(sd.get_optical_property_set(*new_element), nullptr);
}