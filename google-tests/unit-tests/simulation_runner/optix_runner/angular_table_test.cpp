#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include <CspElement.h>
#include <angular_table_registry.h>
#include <optical_properties.hpp>

using OptixCSP::AngularTableRange;
using OptixCSP::AngularTableRegistry;
using OptixCSP::CspElement;
using OptixCSP::MaterialData;
using SolTrace::Data::AngularTablePoint;
using SolTrace::Data::DistributionType;
using SolTrace::Data::InteractionType;
using SolTrace::Data::OpticalPropertySet;
using SolTrace::Data::OpticalSide;

namespace
{

std::shared_ptr<OpticalPropertySet> make_optics(InteractionType type)
{
    auto optics =
        std::make_shared<OpticalPropertySet>(type, 1.0, 1.0, "TestOptics");
    optics->set_properties(
        OpticalSide::Both, DistributionType::NONE, 0.0, 0.0, 0.0, 0.0);
    return optics;
}

const std::vector<AngularTablePoint> kSampleTable = {
    { 0.0, 1.0 }, { 100.0, 0.8 }, { 300.0, 0.4 }
};

} // namespace

// --- AngularTableRegistry --------------------------------------------------

TEST(AngularTableRegistry, ReturnsEmptyRangeWhenNoTableEnabled)
{
    AngularTableRegistry registry;
    auto                 optics = make_optics(InteractionType::REFLECTION);

    const AngularTableRange range = registry.intern(
        1, *optics, OpticalSide::Front, /*use_refraction=*/false);

    EXPECT_EQ(range.offset, 0u);
    EXPECT_EQ(range.count, 0u);
    EXPECT_TRUE(registry.cos_pool().empty());
    EXPECT_TRUE(registry.value_pool().empty());
}

TEST(AngularTableRegistry, InternsReflectivityTableAndMatchesCache)
{
    AngularTableRegistry registry;
    auto                 optics = make_optics(InteractionType::REFLECTION);
    optics->set_reflectivity_table(OpticalSide::Front, kSampleTable);

    const AngularTableRange range = registry.intern(
        1, *optics, OpticalSide::Front, /*use_refraction=*/false);

    ASSERT_EQ(range.count, kSampleTable.size());
    EXPECT_EQ(range.offset, 0u);

    const std::vector<double>& expected_cos =
        optics->get_reflectivity_cosine_cache(OpticalSide::Front);
    const std::vector<double> expected_values =
        optics->get_reflectivity_values(OpticalSide::Front);

    for (std::size_t i = 0; i < range.count; ++i)
    {
        EXPECT_FLOAT_EQ(registry.cos_pool()[range.offset + i],
                        static_cast<float>(expected_cos[i]));
        EXPECT_FLOAT_EQ(registry.value_pool()[range.offset + i],
                        static_cast<float>(expected_values[i]));
    }
}

TEST(AngularTableRegistry, InternsTransmissivityTableOnlyWhenUseRefraction)
{
    AngularTableRegistry registry;
    auto                 optics = make_optics(InteractionType::REFRACTION);
    optics->set_transmissivity_table(OpticalSide::Back, kSampleTable);

    // A reflectivity lookup on a side/optics that only has a transmissivity
    // table enabled must not pick that table up.
    const AngularTableRange refl_range = registry.intern(
        2, *optics, OpticalSide::Back, /*use_refraction=*/false);
    EXPECT_EQ(refl_range.count, 0u);

    const AngularTableRange trans_range = registry.intern(
        3, *optics, OpticalSide::Back, /*use_refraction=*/true);
    ASSERT_EQ(trans_range.count, kSampleTable.size());
}

TEST(AngularTableRegistry, DeduplicatesRepeatedInternCalls)
{
    AngularTableRegistry registry;
    auto                 optics = make_optics(InteractionType::REFLECTION);
    optics->set_reflectivity_table(OpticalSide::Front, kSampleTable);

    const AngularTableRange first = registry.intern(
        7, *optics, OpticalSide::Front, /*use_refraction=*/false);
    const AngularTableRange second = registry.intern(
        7, *optics, OpticalSide::Front, /*use_refraction=*/false);

    EXPECT_EQ(first.offset, second.offset);
    EXPECT_EQ(first.count, second.count);
    EXPECT_EQ(registry.cos_pool().size(), kSampleTable.size());
}

TEST(AngularTableRegistry, TracksSidesIndependently)
{
    AngularTableRegistry registry;
    auto                 optics = make_optics(InteractionType::REFLECTION);
    optics->set_reflectivity_table(OpticalSide::Front, kSampleTable);
    // Back side left without a table.

    const AngularTableRange front_range = registry.intern(
        4, *optics, OpticalSide::Front, /*use_refraction=*/false);
    const AngularTableRange back_range = registry.intern(
        4, *optics, OpticalSide::Back, /*use_refraction=*/false);

    EXPECT_EQ(front_range.count, kSampleTable.size());
    EXPECT_EQ(back_range.count, 0u);
}

TEST(AngularTableRegistry, DistinctIdsAppendSeparatePoolEntries)
{
    AngularTableRegistry registry;
    auto                 optics_a = make_optics(InteractionType::REFLECTION);
    optics_a->set_reflectivity_table(OpticalSide::Front, kSampleTable);

    const std::vector<AngularTablePoint> other_table = { { 50.0, 0.9 },
                                                          { 500.0, 0.1 } };
    auto optics_b = make_optics(InteractionType::REFLECTION);
    optics_b->set_reflectivity_table(OpticalSide::Front, other_table);

    const AngularTableRange range_a = registry.intern(
        10, *optics_a, OpticalSide::Front, /*use_refraction=*/false);
    const AngularTableRange range_b = registry.intern(
        11, *optics_b, OpticalSide::Front, /*use_refraction=*/false);

    EXPECT_EQ(range_a.offset, 0u);
    EXPECT_EQ(range_a.count, kSampleTable.size());
    EXPECT_EQ(range_b.offset, kSampleTable.size());
    EXPECT_EQ(range_b.count, other_table.size());
    EXPECT_EQ(registry.cos_pool().size(),
             kSampleTable.size() + other_table.size());
}

// --- CspElement::toDeviceMaterialDataFront/Back ----------------------------

TEST(CspElementAngularTable, PopulatesMaterialDataForEnabledSideOnly)
{
    auto optics = make_optics(InteractionType::REFLECTION);
    optics->set_reflectivity_table(OpticalSide::Front, kSampleTable);

    CspElement element;
    element.set_optics(optics, /*id=*/1);

    AngularTableRegistry registry;
    const MaterialData   front = element.toDeviceMaterialDataFront(registry);
    const MaterialData   back  = element.toDeviceMaterialDataBack(registry);

    EXPECT_TRUE(front.use_angular_table);
    EXPECT_EQ(front.angular_table_count, kSampleTable.size());

    EXPECT_FALSE(back.use_angular_table);
    EXPECT_EQ(back.angular_table_count, 0u);
}

TEST(CspElementAngularTable, RefractionSelectsTransmissivityTable)
{
    auto optics = make_optics(InteractionType::REFRACTION);
    optics->set_transmissivity_table(OpticalSide::Both, kSampleTable);

    CspElement element;
    element.set_optics(optics, /*id=*/2);

    AngularTableRegistry registry;
    const MaterialData   front = element.toDeviceMaterialDataFront(registry);

    EXPECT_TRUE(front.use_refraction);
    EXPECT_TRUE(front.use_angular_table);
    EXPECT_EQ(front.angular_table_count, kSampleTable.size());
}

TEST(CspElementAngularTable, ElementsSharingOpticsShareRegistryEntry)
{
    auto optics = make_optics(InteractionType::REFLECTION);
    optics->set_reflectivity_table(OpticalSide::Front, kSampleTable);

    CspElement element_a;
    element_a.set_optics(optics, /*id=*/5);
    CspElement element_b;
    element_b.set_optics(optics, /*id=*/5);

    AngularTableRegistry registry;
    const MaterialData   md_a = element_a.toDeviceMaterialDataFront(registry);
    const MaterialData   md_b = element_b.toDeviceMaterialDataFront(registry);

    EXPECT_EQ(md_a.angular_table_offset, md_b.angular_table_offset);
    EXPECT_EQ(md_a.angular_table_count, md_b.angular_table_count);
    EXPECT_EQ(registry.cos_pool().size(), kSampleTable.size());
}
