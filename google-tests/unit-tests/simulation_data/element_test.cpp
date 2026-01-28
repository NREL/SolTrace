#include <gtest/gtest.h>

#include <sstream>

#include <simulation_data_export.hpp>

#include "common.hpp"

TEST(OpticalProperties, OutputOperator)
{
    std::stringstream ss;

    OpticalProperties op;
    op.set_ideal_reflection();
    ss << op;
    EXPECT_GE(ss.str().length(), 0);
}

TEST(Element, ConstructionSmokeTest)
{
    // Smoke tests
    auto elem = SolTrace::Data::make_element<SingleElement>();
    auto comp = SolTrace::Data::make_element<CompositeElement>();
    // auto ve = SolTrace::Data::make_element<VirtualElement>();
    // auto vp = SolTrace::Data::make_element<VirtualPlane>(10.0, 10.0);
    auto st = SolTrace::Data::make_stage(1);
}

TEST(Element, SingleElementAccessors)
{
    SingleElement ref;
    Vector3d zero(0.0, 0.0, 0.0);
    Vector3d khat(0.0, 0.0, 1.0);
    EXPECT_TRUE(is_identical(ref.get_origin_ref(), zero));
    EXPECT_TRUE(is_identical(ref.get_aim_vector_ref(), khat));
    EXPECT_TRUE(is_identical(ref.get_euler_angles(), zero));
    EXPECT_EQ(ref.get_aperture(), nullptr);
    EXPECT_EQ(ref.get_surface(), nullptr);

    EXPECT_TRUE(ref.is_enabled());
    ref.disable();
    EXPECT_FALSE(ref.is_enabled());
    ref.enable();
    EXPECT_TRUE(ref.is_enabled());

    EXPECT_FALSE(ref.is_composite());
    EXPECT_FALSE(ref.is_virtual());
    EXPECT_TRUE(ref.is_single());

    EXPECT_FALSE(ref.is_virtual());
    ref.mark_virtual();
    EXPECT_TRUE(ref.is_virtual());
    ref.unmark_virtual();
    EXPECT_FALSE(ref.is_virtual());

    auto pos = Vector3d(2.0, 1.0, -3.0);
    ref.set_origin(pos);
    EXPECT_TRUE(is_identical(ref.get_origin_ref(), pos));

    auto aim = Vector3d(-1.0, 0.0, 1.0);
    ref.set_aim_vector(aim);
    EXPECT_TRUE(is_identical(ref.get_aim_vector_ref(), aim));

    // auto eulers = Vector3d(0.1, 0.2, -0.3);
    // // ref.set_euler_angles(eulers);
    // EXPECT_TRUE(is_identical(ref.get_euler_angles(), eulers));

    const double ZROT = 0.58;
    ref.set_zrot(ZROT);
    EXPECT_EQ(ref.get_zrot(), ZROT);

    EXPECT_EQ(ref.set_reference_frame_geometry(pos, aim, ZROT), 0);
    EXPECT_FALSE(is_identical(ref.get_euler_angles(), zero));

    const double D = 1.0;
    auto ap = SolTrace::Data::make_aperture<SolTrace::Data::Circle>(D);
    ref.set_aperture(ap);
    auto rap = std::dynamic_pointer_cast<SolTrace::Data::Circle>(ref.get_aperture());
    EXPECT_NE(rap, nullptr);
    EXPECT_EQ(rap->diameter, D);

    const double HA = 0.25;
    auto sp = SolTrace::Data::make_surface<SolTrace::Data::Cone>(HA);
    ref.set_surface(sp);
    auto rsp = std::dynamic_pointer_cast<SolTrace::Data::Cone>(ref.get_surface());
    EXPECT_NE(rsp, nullptr);
    EXPECT_EQ(rsp->half_angle, HA);

    auto opf = ref.get_front_optical_properties();
    auto opb = ref.get_back_optical_properties();
    EXPECT_EQ(opf->transmitivity, opb->transmitivity);
    EXPECT_EQ(opf->reflectivity, opb->reflectivity);
    EXPECT_EQ(opf->slope_error, opb->slope_error);
    EXPECT_EQ(opf->specularity_error, opb->specularity_error);

    OpticalProperties op(InteractionType::REFLECTION,
                         DistributionType::GAUSSIAN,
                         0.75, 0.25, 0.1, 0.001, 1.0, 1.0);
    ref.set_front_optical_properties(op);
    // EXPECT_EQ(*opf, op);
    EXPECT_EQ(opf->transmitivity, op.transmitivity);
    EXPECT_EQ(opf->reflectivity, op.reflectivity);
    EXPECT_EQ(opf->slope_error, op.slope_error);
    EXPECT_EQ(opf->specularity_error, op.specularity_error);

    ref.set_back_optical_properties(op);
    EXPECT_EQ(opb->transmitivity, op.transmitivity);
    EXPECT_EQ(opb->reflectivity, op.reflectivity);
    EXPECT_EQ(opb->slope_error, op.slope_error);
    EXPECT_EQ(opb->specularity_error, op.specularity_error);
}

TEST(Element, VirtualElement)
{
    VirtualElement ve;

    EXPECT_TRUE(ve.is_virtual());

    const double LX = 1.0;
    const double LY = 2.0;
    ve.set_aperture(make_aperture<Rectangle>(LX, LY));
    auto aptr = ve.get_aperture();
    EXPECT_EQ(aptr->get_type(), ApertureType::RECTANGLE);
    auto rptr = std::dynamic_pointer_cast<Rectangle>(aptr);
    EXPECT_NE(rptr, nullptr);
    if (rptr != nullptr)
    {
        EXPECT_EQ(rptr->x_length, LX);
        EXPECT_EQ(rptr->y_length, LY);
    }

    ve.set_surface(make_surface<Flat>());
    auto sptr = ve.get_surface();
    EXPECT_EQ(sptr->get_type(), SurfaceType::FLAT);
    auto fptr = std::dynamic_pointer_cast<Flat>(sptr);
    EXPECT_NE(fptr, nullptr);

    EXPECT_TRUE(ve.is_virtual());
    EXPECT_TRUE(ve.is_single());
    EXPECT_FALSE(ve.is_composite());

    // These functions should have no effects
    OpticalProperties op(InteractionType::REFLECTION,
                         DistributionType::GAUSSIAN,
                         0.75, 0.25, 0.1, 0.001, 1.0, 1.0);
    ve.set_front_optical_properties(op);
    auto opf = ve.get_front_optical_properties();
    EXPECT_EQ(opf->reflectivity, 0.0);
    EXPECT_EQ(opf->slope_error, 0.0);
    EXPECT_EQ(opf->specularity_error, 0.0);
    EXPECT_EQ(opf->transmitivity, 1.0);
    ve.set_back_optical_properties(op);
    auto opb = ve.get_back_optical_properties();
    EXPECT_EQ(opb->reflectivity, 0.0);
    EXPECT_EQ(opb->slope_error, 0.0);
    EXPECT_EQ(opb->specularity_error, 0.0);
    EXPECT_EQ(opb->transmitivity, 1.0);

    return;
}

TEST(Element, VirtualPlane)
{
    const double LX = 5.0;
    const double LY = 2.5;
    VirtualPlane vp(LX, LY);

    EXPECT_TRUE(vp.is_virtual());
    EXPECT_EQ(vp.get_aperture()->get_type(), ApertureType::RECTANGLE);
    EXPECT_EQ(vp.get_surface()->get_type(), SurfaceType::FLAT);

    auto rptr = std::dynamic_pointer_cast<Rectangle>(vp.get_aperture());
    EXPECT_NE(rptr, nullptr);
    if (rptr != nullptr)
    {
        EXPECT_EQ(rptr->x_length, LX);
        EXPECT_EQ(rptr->y_length, LY);
    }

    // These functions should have no effects
    vp.set_aperture(make_aperture<Circle>(2.0));
    EXPECT_EQ(vp.get_aperture()->get_type(), ApertureType::RECTANGLE);
    vp.set_surface(make_surface<Parabola>(1.0, 1.0));
    EXPECT_EQ(vp.get_surface()->get_type(), SurfaceType::FLAT);

    return;
}

TEST(Element, CompositeElementAccessors)
{
    auto cmp = SolTrace::Data::make_element<CompositeElement>();
    EXPECT_TRUE(cmp->is_composite());

    // Things that should be empty...
    EXPECT_EQ(cmp->get_aperture(), nullptr);
    EXPECT_EQ(cmp->get_surface(), nullptr);
    EXPECT_EQ(cmp->get_front_optical_properties(), nullptr);
    EXPECT_EQ(cmp->get_back_optical_properties(), nullptr);
    const aperture_ptr ap = cmp->get_aperture();
    EXPECT_EQ(ap, nullptr);
    const surface_ptr sp = cmp->get_surface();
    EXPECT_EQ(sp, nullptr);
    const OpticalProperties *op = cmp->get_back_optical_properties();
    EXPECT_EQ(op, nullptr);
    op = cmp->get_front_optical_properties();
    EXPECT_EQ(op, nullptr);

    // These should do nothing...
    cmp->set_front_optical_properties(OpticalProperties());
    cmp->set_back_optical_properties(OpticalProperties());

    // Add/remove/change elements
    const int NUM_ELEMENTS = 4;
    for (int i = 0; i < NUM_ELEMENTS; ++i)
    {
        auto elem = make_configured_element();
        cmp->add_element(elem);
        // EXPECT_EQ(elem->get_stage(), STAGE);
    }
    EXPECT_EQ(cmp->get_number_of_elements(), NUM_ELEMENTS);

    cmp->remove_element(0);
    EXPECT_EQ(cmp->get_number_of_elements(), NUM_ELEMENTS - 1);
    EXPECT_EQ(cmp->get_element(0), nullptr);
    auto elem = make_configured_element();
    auto id = cmp->add_element(elem);
    EXPECT_EQ(cmp->get_element(id).get(), elem.get());

    auto elem2 = make_configured_element();
    EXPECT_NE(elem.get(), elem2.get());
    cmp->replace_element(id, elem2);
    EXPECT_EQ(cmp->get_element(id).get(), elem2.get());

    // Recursive CompositeElements -- this is disallowed!
    auto cmp2 = SolTrace::Data::make_element<CompositeElement>();
    auto cid = cmp->add_element(cmp2);
    EXPECT_FALSE(SolTrace::Data::Element::is_success(cid));
    EXPECT_EQ(cid, SolTrace::Data::ELEMENT_INVALID_SETUP);
    EXPECT_EQ(cmp->get_number_of_elements(), NUM_ELEMENTS);

    elem2->enable();
    cmp->disable();
    EXPECT_FALSE(cmp->is_enabled());
    EXPECT_FALSE(elem2->is_enabled());
    cmp->enable();
    EXPECT_TRUE(cmp->is_enabled());
    EXPECT_TRUE(elem2->is_enabled());

    EXPECT_FALSE(elem2->is_virtual());
    EXPECT_FALSE(cmp->is_virtual());
    cmp->mark_virtual();
    EXPECT_TRUE(cmp->is_virtual());
    EXPECT_TRUE(elem2->is_virtual());
    cmp->unmark_virtual();
    EXPECT_FALSE(cmp->is_virtual());
    EXPECT_FALSE(elem2->is_virtual());

    cmp->mark_virtual();
    auto elem3 = make_configured_element();
    EXPECT_FALSE(elem3->is_virtual());
    cmp->add_element(elem3);
    EXPECT_TRUE(elem3->is_virtual());
    cmp->unmark_virtual();
    EXPECT_FALSE(elem3->is_virtual());

    // Check that pass through functions are hooked up correctly
    auto iter = cmp->get_iterator();
    while (!cmp->is_at_end(iter))
    {
        ++iter;
    }
    EXPECT_TRUE(cmp->is_at_end(iter));

    auto citer = cmp->get_const_iterator();
    while (!cmp->is_at_end(citer))
    {
        ++citer;
    }
    EXPECT_TRUE(cmp->is_at_end(citer));
}

TEST(Element, StageElementAccessors)
{
    const int_fast64_t STAGE = 10;
    const int_fast64_t RESET_STAGE = 20;
    auto st1 = SolTrace::Data::make_stage(10);

    auto el1 = make_configured_element();
    auto cmp1 = SolTrace::Data::make_element<CompositeElement>();
    auto sub1 = make_configured_element();
    auto sub2 = make_configured_element();
    auto sub3 = make_configured_element();
    EXPECT_TRUE(SolTrace::Data::Element::is_success(cmp1->add_element(sub1)));
    EXPECT_TRUE(SolTrace::Data::Element::is_success(cmp1->add_element(sub2)));
    EXPECT_TRUE(SolTrace::Data::Element::is_success(cmp1->add_element(sub3)));

    EXPECT_TRUE(SolTrace::Data::Element::is_success(st1->add_element(el1)));
    EXPECT_TRUE(SolTrace::Data::Element::is_success(st1->add_element(cmp1)));

    EXPECT_EQ(st1->get_number_of_elements(), 4);
    EXPECT_EQ(st1->get_stage(), STAGE);
    EXPECT_EQ(el1->get_stage(), STAGE);
    EXPECT_EQ(cmp1->get_stage(), STAGE);
    EXPECT_EQ(sub1->get_stage(), STAGE);

    st1->set_stage(RESET_STAGE);
    EXPECT_EQ(st1->get_stage(), RESET_STAGE);
    EXPECT_EQ(el1->get_stage(), RESET_STAGE);
    EXPECT_EQ(cmp1->get_stage(), RESET_STAGE);
    EXPECT_EQ(sub1->get_stage(), RESET_STAGE);
}

TEST(Element, CoordinateComputationsIdentity)
{
    auto el = make_configured_element();
    auto st = SolTrace::Data::make_stage(0);
    Vector3d origin(0.0, 0.0, 0.0);
    Vector3d aim(0.0, 0.0, 1.0);
    Vector3d local(1.0, 1.0, 1.0);
    Vector3d ref(1.0, 1.0, 1.0);
    Vector3d result;
    result.zero();
    double zrot = 0.0;

    el->set_reference_frame_geometry(origin, aim, zrot);
    st->set_reference_frame_geometry(origin, aim, zrot);

    // Conversion test
    el->convert_local_to_global(result, local);
    EXPECT_TRUE(is_identical(ref, result));
    result.zero();
    el->convert_local_to_stage(result, local);
    EXPECT_TRUE(is_identical(ref, result));
    result.zero();
    el->convert_local_to_reference(result, local);
    EXPECT_TRUE(is_identical(ref, result));
    result.zero();

    st->convert_local_to_global(result, local);
    EXPECT_TRUE(is_identical(ref, result));
    result.zero();
    st->convert_local_to_stage(result, local);
    EXPECT_TRUE(is_identical(ref, result));
    result.zero();
    st->convert_local_to_reference(result, local);
    EXPECT_TRUE(is_identical(ref, result));
    result.zero();

    st->add_element(el);
    el->convert_local_to_global(result, local);
    EXPECT_TRUE(is_identical(ref, result));
    result.zero();
    el->convert_local_to_stage(result, local);
    EXPECT_TRUE(is_identical(ref, result));
    result.zero();
    el->convert_local_to_reference(result, local);
    EXPECT_TRUE(is_identical(ref, result));
    result.zero();

    // Origin location tests
    EXPECT_TRUE(is_identical(origin, st->get_origin_ref()));
    EXPECT_TRUE(is_identical(origin, st->get_origin_stage()));
    EXPECT_TRUE(is_identical(origin, st->get_origin_global()));

    EXPECT_TRUE(is_identical(origin, el->get_origin_ref()));
    EXPECT_TRUE(is_identical(origin, el->get_origin_stage()));
    EXPECT_TRUE(is_identical(origin, el->get_origin_global()));

    // Aim vector conversion tests
    EXPECT_TRUE(is_identical(aim, st->get_aim_vector_ref()));
    EXPECT_TRUE(is_identical(aim, st->get_aim_vector_stage()));
    EXPECT_TRUE(is_identical(aim, st->get_aim_vector_global()));

    EXPECT_TRUE(is_identical(aim, el->get_aim_vector_ref()));
    EXPECT_TRUE(is_identical(aim, el->get_aim_vector_stage()));
    EXPECT_TRUE(is_identical(aim, el->get_aim_vector_global()));

    // Operators
    Matrix3d Q;
    Q.identity();
    Matrix3d RtoL = el->get_reference_to_local();
    EXPECT_TRUE(is_identical(RtoL, Q));
    Matrix3d LtoR = el->get_local_to_reference();
    EXPECT_TRUE(is_identical(LtoR, Q));
    Matrix3d StoL = el->get_stage_to_local();
    EXPECT_TRUE(is_identical(StoL, Q));
    Matrix3d LtoS = el->get_local_to_stage();
    EXPECT_TRUE(is_identical(LtoS, Q));
    Matrix3d GtoL = el->get_global_to_local();
    EXPECT_TRUE(is_identical(GtoL, Q));
    Matrix3d LtoG = el->get_local_to_global();
    EXPECT_TRUE(is_identical(LtoG, Q));
}

TEST(Element, CoordinateComputationsRotations)
{
    using SolTrace::Data::MatrixTranspose;
    using SolTrace::Data::PI;

    // **** Setup Answers **** //
    // Origin
    Vector3d Origin;
    Origin.zero();
    // Coordinate transform matrix
    Matrix3d Q1;
    Q1.set_value(0, 0, 0.5);
    Q1.set_value(1, 0, sqrt(3.0) / 2.0);
    Q1.set_value(2, 0, 0.0);
    Q1.set_value(0, 1, -1.0 / sqrt(2.0));
    Q1.set_value(1, 1, 1.0 / sqrt(6.0));
    Q1.set_value(2, 1, -1.0 / sqrt(3.0));
    Q1.set_value(0, 2, -0.5);
    Q1.set_value(1, 2, sqrt(3.0) / 6.0);
    Q1.set_value(2, 2, sqrt(6.0) / 3.0);
    Matrix3d Q1t;
    MatrixTranspose(Q1.data, 3, Q1t.data);
    // Corresponding Euler angles in radians
    const double a1 = 0.0;
    const double b1 = asin(-1.0 / sqrt(3.0));
    const double g1 = acos(1.0 / cos(b1) * 1.0 / sqrt(6.0)); // approximately 0.615
    // Corresponding aim vector (local z-axis in reference coordinates)
    // Vector3d aim1(0.0, -sqrt(3.0) / 2.0, sqrt(2.0 / 3.0));
    Vector3d aim1(0.0, -1.0 / sqrt(3.0), sqrt(2.0 / 3.0));

    // Z-Rotation is the last of the Euler angles but in degrees
    const double zrot1 = g1 * 180.0 / PI;

    Matrix3d Q2;
    Q2.set_value(0, 0, (sqrt(8.0) + sqrt(6.0)) / 8.0);
    Q2.set_value(0, 1, -0.75);
    Q2.set_value(0, 2, (sqrt(6.0) - sqrt(8.0)) / 8.0);
    Q2.set_value(1, 0, (2.0 * sqrt(6.0) - sqrt(2.0)) / 8.0);
    Q2.set_value(1, 1, sqrt(3.0) / 4.0);
    Q2.set_value(1, 2, (-2.0 * sqrt(6.0) - sqrt(2.0)) / 8.0);
    Q2.set_value(2, 0, sqrt(6.0) / 4.0);
    Q2.set_value(2, 1, 0.5);
    Q2.set_value(2, 2, sqrt(6.0) / 4.0);
    Matrix3d Q2t;
    MatrixTranspose(Q2.data, 3, Q2t.data);
    // Corresponding Euler angles in radians
    const double a2 = PI / 4.0;
    const double b2 = PI / 6.0;
    const double g2 = PI / 3.0;
    // Corresponding aim vector (local z-axis in reference coordinates)
    Vector3d aim2(sqrt(3.0 / 8.0), 0.5, sqrt(3.0 / 8.0));

    // Z-Rotation is the last of the Euler angles but in degrees
    const double zrot2 = 60.0;

    // **** Setup Elements **** //
    auto el = make_configured_element();
    el->set_reference_frame_geometry(Origin, aim1, zrot1);

    auto st = SolTrace::Data::make_stage(0);
    st->set_reference_frame_geometry(Origin, aim2, zrot2);
    st->add_element(el);

    // **** Tests **** //
    const double TOL = 1e-12;
    Vector3d scratch;
    Vector3d result_vec;
    Matrix3d result_mat;
    // Origin tests
    EXPECT_TRUE(is_identical(el->get_origin_stage(), el->get_origin_ref()));
    EXPECT_TRUE(is_identical(el->get_origin_ref(), Origin));
    EXPECT_TRUE(is_identical(el->get_origin_global(), Origin));

    // Euler angles tests
    result_vec = el->get_euler_angles();
    EXPECT_NEAR(result_vec[0], a1, TOL);
    EXPECT_NEAR(result_vec[1], b1, TOL);
    EXPECT_NEAR(result_vec[2], g1, TOL);
    result_vec = st->get_euler_angles();
    EXPECT_NEAR(result_vec[0], a2, TOL);
    EXPECT_NEAR(result_vec[1], b2, TOL);
    EXPECT_NEAR(result_vec[2], g2, TOL);

    // Aim vector tests
    EXPECT_TRUE(is_identical(el->get_aim_vector_ref(), aim1, TOL));
    EXPECT_TRUE(is_identical(el->get_aim_vector_stage(), el->get_aim_vector_ref(), TOL));
    matrix_vector_product(Q2t, aim1, result_vec);
    EXPECT_TRUE(is_identical(el->get_aim_vector_global(), result_vec, TOL));

    Vector3d v_local(-1.0, 2.0, 4.0);
    Vector3d v_stage;
    Vector3d v_global;
    matrix_vector_product(Q1t, v_local, v_stage);
    matrix_vector_product(Q2t, v_stage, v_global);

    Vector3d result;

    el->convert_local_to_stage(result, v_local);
    EXPECT_TRUE(is_identical(result, v_stage, TOL));
    result.zero();

    el->convert_stage_to_local(result, v_stage);
    EXPECT_TRUE(is_identical(result, v_local, TOL));
    result.zero();

    el->convert_local_to_global(result, v_local);
    EXPECT_TRUE(is_identical(result, v_global, TOL));
    result.zero();

    el->convert_global_to_local(result, v_global);
    EXPECT_TRUE(is_identical(result, v_local, TOL));
}

TEST(Element, CoordinateComputationsTranslations)
{
    // **** Setup Answers **** //
    // Origin
    Vector3d Origin1(1.0, 2.0, 3.0);
    // Corresponding Euler angles in radians
    const double a = 0.0;
    const double b = 0.0;
    const double g = 0.0;
    // Corresponding aim vector (local z-axis in reference coordinates)
    // Vector3d aim1(0.0, -sqrt(3.0) / 2.0, sqrt(2.0 / 3.0));
    Vector3d aim1(0.0, 0.0, 1.0);
    vector_add(1.0, Origin1, 1.0, aim1);

    // Z-Rotation is the last of the Euler angles but in degrees
    const double zrot = 0.0;

    // Origin
    Vector3d Origin2(-3.0, 1.0, -5.0);
    // Corresponding aim vector (local z-axis in reference coordinates)
    Vector3d aim2(0.0, 0.0, 1.0);
    vector_add(1.0, Origin2, 1.0, aim2);

    // **** Setup Elements **** //
    auto el = make_configured_element();
    el->set_reference_frame_geometry(Origin1, aim1, zrot);

    auto st = SolTrace::Data::make_stage(0);
    st->set_reference_frame_geometry(Origin2, aim2, zrot);
    st->add_element(el);

    // **** Tests **** //
    const double TOL = 1e-12;
    Vector3d scratch;
    Vector3d result_vec;
    Matrix3d result_mat;
    // Origin tests
    EXPECT_TRUE(is_identical(el->get_origin_stage(), el->get_origin_ref()));
    EXPECT_TRUE(is_identical(el->get_origin_ref(), Origin1));
    vector_add(1.0, Origin1, 1.0, Origin2, result_vec);
    EXPECT_TRUE(is_identical(el->get_origin_global(), result_vec));

    // Euler angles tests
    result_vec = el->get_euler_angles();
    EXPECT_NEAR(result_vec[0], a, TOL);
    EXPECT_NEAR(result_vec[1], b, TOL);
    EXPECT_NEAR(result_vec[2], g, TOL);
    result_vec = st->get_euler_angles();
    EXPECT_NEAR(result_vec[0], a, TOL);
    EXPECT_NEAR(result_vec[1], b, TOL);
    EXPECT_NEAR(result_vec[2], g, TOL);

    // Aim vector tests
    EXPECT_TRUE(is_identical(el->get_aim_vector_ref(), aim1, TOL));
    EXPECT_TRUE(is_identical(el->get_aim_vector_stage(), el->get_aim_vector_ref(), TOL));
    vector_add(1.0, Origin2, 1.0, aim1, result_vec);
    EXPECT_TRUE(is_identical(el->get_aim_vector_global(), result_vec, TOL));

    Vector3d v_local(-1.0, 2.0, 4.0);
    Vector3d v_stage;
    Vector3d v_global;
    vector_add(1.0, Origin1, 1.0, v_local, v_stage);
    vector_add(1.0, Origin2, 1.0, v_stage, v_global);

    Vector3d result;

    el->convert_local_to_stage(result, v_local);
    EXPECT_TRUE(is_identical(result, v_stage, TOL));
    result.zero();

    el->convert_stage_to_local(result, v_stage);
    EXPECT_TRUE(is_identical(result, v_local, TOL));
    result.zero();

    el->convert_local_to_global(result, v_local);
    EXPECT_TRUE(is_identical(result, v_global, TOL));
    result.zero();

    el->convert_global_to_local(result, v_global);
    EXPECT_TRUE(is_identical(result, v_local, TOL));
}

TEST(Element, CoordinateComputations)
{
    using SolTrace::Data::MatrixTranspose;
    using SolTrace::Data::PI;

    // **** Setup Answers **** //
    // Origin
    Vector3d Origin1(1.0, 2.0, 3.0);
    // Coordinate transform matrix -- stage to local
    Matrix3d Q1;
    Q1.set_value(0, 0, 0.5);
    Q1.set_value(1, 0, sqrt(3.0) / 2.0);
    Q1.set_value(2, 0, 0.0);
    Q1.set_value(0, 1, -1.0 / sqrt(2.0));
    Q1.set_value(1, 1, 1.0 / sqrt(6.0));
    Q1.set_value(2, 1, -1.0 / sqrt(3.0));
    Q1.set_value(0, 2, -0.5);
    Q1.set_value(1, 2, sqrt(3.0) / 6.0);
    Q1.set_value(2, 2, sqrt(6.0) / 3.0);
    // Local to stage matrix
    Matrix3d Q1t;
    MatrixTranspose(Q1.data, 3, Q1t.data);
    // Corresponding Euler angles in radians
    const double a1 = 0.0;
    const double b1 = asin(-1.0 / sqrt(3.0));
    const double g1 = acos(1.0 / cos(b1) * 1.0 / sqrt(6.0)); // approximately 0.615
    // Corresponding aim vector (local z-axis in reference coordinates)
    // Vector3d aim1(0.0, -sqrt(3.0) / 2.0, sqrt(2.0 / 3.0));
    Vector3d aim1(0.0, -1.0 / sqrt(3.0), sqrt(2.0 / 3.0));
    vector_add(1.0, Origin1, 1.0, aim1);

    // Z-Rotation is the last of the Euler angles but in degrees
    const double zrot1 = g1 * 180.0 / PI;

    // Origin
    Vector3d Origin2(-3.0, 1.0, -5.0);
    // Global to stage matrix
    Matrix3d Q2;
    Q2.set_value(0, 0, (sqrt(8.0) + sqrt(6.0)) / 8.0);
    Q2.set_value(0, 1, -0.75);
    Q2.set_value(0, 2, (sqrt(6.0) - sqrt(8.0)) / 8.0);
    Q2.set_value(1, 0, (2.0 * sqrt(6.0) - sqrt(2.0)) / 8.0);
    Q2.set_value(1, 1, sqrt(3.0) / 4.0);
    Q2.set_value(1, 2, (-2.0 * sqrt(6.0) - sqrt(2.0)) / 8.0);
    Q2.set_value(2, 0, sqrt(6.0) / 4.0);
    Q2.set_value(2, 1, 0.5);
    Q2.set_value(2, 2, sqrt(6.0) / 4.0);
    // Stage to global matrix
    Matrix3d Q2t;
    MatrixTranspose(Q2.data, 3, Q2t.data);
    // Corresponding Euler angles in radians
    const double a2 = PI / 4.0;
    const double b2 = PI / 6.0;
    const double g2 = PI / 3.0;
    // Corresponding aim vector (local z-axis in reference coordinates)
    Vector3d aim2(sqrt(3.0 / 8.0), 0.5, sqrt(3.0 / 8.0));
    vector_add(1.0, Origin2, 1.0, aim2);

    // Z-Rotation is the last of the Euler angles but in degrees
    const double zrot2 = 60.0;

    // **** Setup Elements **** //
    auto el = make_configured_element();
    el->set_reference_frame_geometry(Origin1, aim1, zrot1);

    auto st = SolTrace::Data::make_stage(0);
    st->set_reference_frame_geometry(Origin2, aim2, zrot2);
    st->add_element(el);

    // **** Tests **** //
    const double TOL = 1e-12;
    Vector3d scratch;
    Vector3d result_vec;
    Matrix3d result_mat;
    // Origin tests
    EXPECT_TRUE(is_identical(el->get_origin_stage(), el->get_origin_ref()));
    EXPECT_TRUE(is_identical(el->get_origin_ref(), Origin1));
    // vector_add(1.0, Origin1, 1.0, Origin2, result_vec);
    matrix_vector_product(Q2t, Origin1, result_vec);
    vector_add(1.0, Origin2, 1.0, result_vec);
    EXPECT_TRUE(is_identical(el->get_origin_global(), result_vec));

    // Euler angles tests
    result_vec = el->get_euler_angles();
    EXPECT_NEAR(result_vec[0], a1, TOL);
    EXPECT_NEAR(result_vec[1], b1, TOL);
    EXPECT_NEAR(result_vec[2], g1, TOL);
    result_vec = st->get_euler_angles();
    EXPECT_NEAR(result_vec[0], a2, TOL);
    EXPECT_NEAR(result_vec[1], b2, TOL);
    EXPECT_NEAR(result_vec[2], g2, TOL);

    // Aim vector tests
    EXPECT_TRUE(is_identical(el->get_aim_vector_ref(), aim1, TOL));
    EXPECT_TRUE(is_identical(el->get_aim_vector_stage(), el->get_aim_vector_ref(), TOL));
    matrix_vector_product(Q2t, aim1, result_vec);
    vector_add(1.0, Origin2, 1.0, result_vec);
    EXPECT_TRUE(is_identical(el->get_aim_vector_global(), result_vec, TOL));

    Vector3d v_local(-1.0, 2.0, 4.0);
    Vector3d v_stage;
    Vector3d v_global;
    matrix_vector_product(Q1t, v_local, v_stage);
    vector_add(1.0, Origin1, 1.0, v_stage);
    matrix_vector_product(Q2t, v_stage, v_global);
    vector_add(1.0, Origin2, 1.0, v_global);

    Vector3d result;

    el->convert_local_to_stage(result, v_local);
    EXPECT_TRUE(is_identical(result, v_stage, TOL));
    result.zero();

    el->convert_stage_to_local(result, v_stage);
    EXPECT_TRUE(is_identical(result, v_local, TOL));
    result.zero();

    el->convert_local_to_global(result, v_local);
    EXPECT_TRUE(is_identical(result, v_global, TOL));
    result.zero();

    el->convert_global_to_local(result, v_global);
    EXPECT_TRUE(is_identical(result, v_local, TOL));
    result.zero();

    // Note: Global to reference and reference to global
    // use other conversion routines so just need to test
    // that the if-else statement is correct.

    // Stage coordinates are the reference coordinates
    el->convert_global_to_reference(result, v_global);
    EXPECT_TRUE(is_identical(result, v_stage, TOL));
    result.zero();

    el->convert_reference_to_global(result, v_stage);
    EXPECT_TRUE(is_identical(result, v_global, TOL));
    result.zero();

    // Global coordinates are reference coordinates
    st->convert_global_to_reference(result, v_global);
    EXPECT_TRUE(is_identical(result, v_global, TOL));
    result.zero();

    st->convert_reference_to_global(result, v_global);
    EXPECT_TRUE(is_identical(result, v_global, TOL));
    result.zero();
}

TEST(Element, VectorCoordinateComputations)
{
    using SolTrace::Data::MatrixTranspose;
    using SolTrace::Data::PI;

    // **** Setup Answers **** //
    // Origin
    Vector3d Origin1(1.0, 2.0, 3.0);
    // Coordinate transform matrix -- stage to local
    Matrix3d Q1;
    Q1.set_value(0, 0, 0.5);
    Q1.set_value(1, 0, sqrt(3.0) / 2.0);
    Q1.set_value(2, 0, 0.0);
    Q1.set_value(0, 1, -1.0 / sqrt(2.0));
    Q1.set_value(1, 1, 1.0 / sqrt(6.0));
    Q1.set_value(2, 1, -1.0 / sqrt(3.0));
    Q1.set_value(0, 2, -0.5);
    Q1.set_value(1, 2, sqrt(3.0) / 6.0);
    Q1.set_value(2, 2, sqrt(6.0) / 3.0);
    // Local to stage matrix
    Matrix3d Q1t;
    MatrixTranspose(Q1.data, 3, Q1t.data);
    // Corresponding Euler angles in radians
    const double a1 = 0.0;
    const double b1 = asin(-1.0 / sqrt(3.0));
    const double g1 = acos(1.0 / cos(b1) * 1.0 / sqrt(6.0)); // approximately 0.615
    // Corresponding aim vector (local z-axis in reference coordinates)
    // Vector3d aim1(0.0, -sqrt(3.0) / 2.0, sqrt(2.0 / 3.0));
    Vector3d aim1(0.0, -1.0 / sqrt(3.0), sqrt(2.0 / 3.0));
    vector_add(1.0, Origin1, 1.0, aim1);

    // Z-Rotation is the last of the Euler angles but in degrees
    const double zrot1 = g1 * 180.0 / PI;

    // Origin
    Vector3d Origin2(-3.0, 1.0, -5.0);
    // Global to stage matrix
    Matrix3d Q2;
    Q2.set_value(0, 0, (sqrt(8.0) + sqrt(6.0)) / 8.0);
    Q2.set_value(0, 1, -0.75);
    Q2.set_value(0, 2, (sqrt(6.0) - sqrt(8.0)) / 8.0);
    Q2.set_value(1, 0, (2.0 * sqrt(6.0) - sqrt(2.0)) / 8.0);
    Q2.set_value(1, 1, sqrt(3.0) / 4.0);
    Q2.set_value(1, 2, (-2.0 * sqrt(6.0) - sqrt(2.0)) / 8.0);
    Q2.set_value(2, 0, sqrt(6.0) / 4.0);
    Q2.set_value(2, 1, 0.5);
    Q2.set_value(2, 2, sqrt(6.0) / 4.0);
    // Stage to global matrix
    Matrix3d Q2t;
    MatrixTranspose(Q2.data, 3, Q2t.data);
    // Corresponding Euler angles in radians
    const double a2 = PI / 4.0;
    const double b2 = PI / 6.0;
    const double g2 = PI / 3.0;
    // Corresponding aim vector (local z-axis in reference coordinates)
    Vector3d aim2(sqrt(3.0 / 8.0), 0.5, sqrt(3.0 / 8.0));
    vector_add(1.0, Origin2, 1.0, aim2);

    // Z-Rotation is the last of the Euler angles but in degrees
    const double zrot2 = 60.0;

    // **** Setup Elements **** //
    auto el = make_configured_element();
    el->set_reference_frame_geometry(Origin1, aim1, zrot1);

    auto st = SolTrace::Data::make_stage(0);
    st->set_reference_frame_geometry(Origin2, aim2, zrot2);
    st->add_element(el);

    // **** Tests **** //
    const double TOL = 1e-12;
    Vector3d scratch;
    Vector3d result_vec;
    Matrix3d result_mat;

    // Euler angles tests
    result_vec = el->get_euler_angles();
    EXPECT_NEAR(result_vec[0], a1, TOL);
    EXPECT_NEAR(result_vec[1], b1, TOL);
    EXPECT_NEAR(result_vec[2], g1, TOL);
    result_vec = st->get_euler_angles();
    EXPECT_NEAR(result_vec[0], a2, TOL);
    EXPECT_NEAR(result_vec[1], b2, TOL);
    EXPECT_NEAR(result_vec[2], g2, TOL);

    Vector3d v_local(-1.0, 2.0, 4.0);
    Vector3d v_stage;
    Vector3d v_global;
    matrix_vector_product(Q1t, v_local, v_stage);
    matrix_vector_product(Q2t, v_stage, v_global);

    Vector3d result;

    el->convert_vector_local_to_stage(result, v_local);
    EXPECT_TRUE(is_identical(result, v_stage, TOL));
    result.zero();

    el->convert_vector_stage_to_local(result, v_stage);
    EXPECT_TRUE(is_identical(result, v_local, TOL));
    result.zero();

    el->convert_vector_local_to_global(result, v_local);
    EXPECT_TRUE(is_identical(result, v_global, TOL));
    result.zero();

    el->convert_vector_global_to_local(result, v_global);
    EXPECT_TRUE(is_identical(result, v_local, TOL));
    result.zero();

    // Note: Global to reference and reference to global
    // use other conversion routines so just need to test
    // that the if-else statement is correct.

    // Stage coordinates are the reference coordinates
    el->convert_vector_global_to_reference(result, v_global);
    EXPECT_TRUE(is_identical(result, v_stage, TOL));
    result.zero();

    el->convert_vector_reference_to_global(result, v_stage);
    EXPECT_TRUE(is_identical(result, v_global, TOL));
    result.zero();

    // Global coordinates are reference coordinates
    st->convert_vector_global_to_reference(result, v_global);
    EXPECT_TRUE(is_identical(result, v_global, TOL));
    result.zero();

    st->convert_vector_reference_to_global(result, v_global);
    EXPECT_TRUE(is_identical(result, v_global, TOL));
    result.zero();
}

TEST(Element, SingleElementEnforceUserFieldsSet)
{
    auto elem = SolTrace::Data::make_element<SingleElement>();

    // SingleElement requires aperture and surface to be set
    // Test that it throws when aperture is missing
    EXPECT_THROW(elem->enforce_user_fields_set(), std::invalid_argument);

    // Set aperture but not surface - should still throw
    elem->set_aperture(SolTrace::Data::make_aperture<SolTrace::Data::Circle>(1.0));
    EXPECT_THROW(elem->enforce_user_fields_set(), std::invalid_argument);

    // Set both aperture and surface - should not throw
    elem->set_surface(SolTrace::Data::make_surface<SolTrace::Data::Flat>());
    EXPECT_NO_THROW(elem->enforce_user_fields_set());

    // Test with optical properties set as well
    OpticalProperties op(InteractionType::REFLECTION,
                         DistributionType::GAUSSIAN,
                         0.75, 0.25, 0.1, 0.001, 1.0, 1.0);
    elem->set_front_optical_properties(op);
    elem->set_back_optical_properties(op);

    // Should still not throw
    EXPECT_NO_THROW(elem->enforce_user_fields_set());
}

TEST(Element, CompositeElementEnforceUserFieldsSet)
{
    auto comp = SolTrace::Data::make_element<CompositeElement>();

    // CompositeElement requires at least one subelement
    // Test that it throws when no subelements are present
    EXPECT_THROW(comp->enforce_user_fields_set(), std::invalid_argument);

    // Add a child element that is not properly configured
    auto elem1 = SolTrace::Data::make_element<SingleElement>();
    EXPECT_THROW(comp->add_element(elem1), std::invalid_argument);

    // Add properly configured child elements
    auto elem2 = SolTrace::Data::make_element<SingleElement>();
    elem2->set_aperture(SolTrace::Data::make_aperture<SolTrace::Data::Rectangle>(2.0, 3.0));
    elem2->set_surface(SolTrace::Data::make_surface<SolTrace::Data::Parabola>(1.0, 2.0));
    EXPECT_NO_THROW(comp->add_element(elem2));

    // Should still not throw for the CompositeElement
    EXPECT_NO_THROW(comp->enforce_user_fields_set());
}
