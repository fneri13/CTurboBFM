#include <iostream>
#include "../../../include/CMesh.hpp"
#include "../../../include/Config.hpp"
#include "gtest/gtest.h"


TEST(CMeshTest, TestComputeVolumes2D) {
    Config config("input.ini");
    CMesh mesh(config);

    auto volumes = mesh.getVolumes();
    auto ni = volumes.sizeI();
    auto nj = volumes.sizeJ();

    std::vector<FloatType> testVolumes = {volumes(0,0,0),
                                          volumes(1,0,0),
                                          volumes(1,1,0),
                                          volumes(ni-1,nj-1,0),
                                          volumes(ni-2,nj-2,0)};
    
    std::vector<FloatType> expectedVolumes = {0.25,
                                              0.5,
                                              1.0,
                                              0.25,
                                              1.0};                                      
    
    for (int i=0; i<testVolumes.size(); i++){
        ASSERT_DOUBLE_EQ(testVolumes[i], expectedVolumes[i]);
    }
}

TEST(CMeshTest, TestReadMesh) {
    Config config("input.ini");
    CMesh mesh(config);
    
    auto vertices = mesh.getVertices();

    std::vector<size_t> iTest {0,1,2,3,4,5};
    std::vector<size_t> jTest {0,1,2,3,4,5};
    size_t kTest = 0;

    for (size_t i=0; i<iTest.size(); i++){
        for (size_t j=0; j<jTest.size(); j++){
            ASSERT_DOUBLE_EQ(vertices(iTest[i],jTest[j],kTest).x(), static_cast<double>(iTest[i]));
            ASSERT_DOUBLE_EQ(vertices(iTest[i],jTest[j],kTest).y(), static_cast<double>(jTest[j]));
            ASSERT_DOUBLE_EQ(vertices(iTest[i],jTest[j],kTest).z(), static_cast<double>(kTest));
        }
    }
}

TEST(CMeshTest, TestDualGrid2DCartesian) {
    Config config("input.ini");
    CMesh mesh(config);
    
    auto dualNodes = mesh.getDualNodes();

    std::vector<Vector3D> testNodes {dualNodes(0,0,0),
                                     dualNodes(1,0,0),
                                     dualNodes(1,1,0)};

    std::vector<Vector3D> expectedNodes {{0.0, 0.0, -0.5},
                                         {0.5, 0.0, -0.5},
                                         {0.5, 0.5, -0.5}};

    for (size_t i=0; i<testNodes.size(); i++){
        ASSERT_DOUBLE_EQ(testNodes[i].x(), expectedNodes[i].x());
        ASSERT_DOUBLE_EQ(testNodes[i].y(), expectedNodes[i].y());
        ASSERT_DOUBLE_EQ(testNodes[i].z(), expectedNodes[i].z());
    }
}

TEST(CMeshTest, TestComputeSurfaceAndCenters) {
    Config config("input.ini");
    CMesh mesh(config);
    
    // define test points, every new element is a new test
    std::vector<Vector3D> points1 = {Vector3D(0.0, 0.0, 0.0),
                                     Vector3D(0.0, 0.0, 0.0),
                                     Vector3D(0.0, 0.0, 0.0),
                                     Vector3D(0.0, 0.0, 0.0)};

    std::vector<Vector3D> points2 = {Vector3D(1.0, 0.0, 0.0),
                                     Vector3D(0.0, 1.0, 0.0),
                                     Vector3D(0.0, 0.0, 1.0),
                                     Vector3D(0.5, 0.0, 0.0)};

    std::vector<Vector3D> points3 = {Vector3D(1.0, 1.0, 0.0),
                                     Vector3D(0.0, 1.0, 1.0),
                                     Vector3D(0.0, 1.0, 1.0),
                                     Vector3D(0.5, 0.0, 0.5)};

    std::vector<Vector3D> points4 = {Vector3D(0.0, 1.0, 0.0),
                                     Vector3D(0.0, 0.0, 1.0),
                                     Vector3D(0.0, 1.0, 0.0),
                                     Vector3D(0.0, 0.0, 0.5)};

    // defined expected results
    std::vector<Vector3D> expectedSurfaces = {Vector3D(0.0, 0.0, 1.0),
                                              Vector3D(1.0, 0.0, 0.0),
                                              Vector3D(-1.0, 0.0, 0.0),
                                              Vector3D(0.0, -0.25, 0.0)};

    std::vector<Vector3D> expectedCenters = {Vector3D(0.5, 0.5, 0.0),
                                             Vector3D(0.0, 0.5, 0.5),
                                             Vector3D(0.0, 0.5, 0.5),
                                             Vector3D(0.25, 0.0, 0.25)};

    // compute the results
    std::vector<Vector3D> computedSurfaces(points1.size());
    std::vector<Vector3D> computedCenters(points1.size()); 
    for (size_t i = 0; i < points1.size(); ++i) {
        mesh.computeSurfaceVectorAndCG(points1[i], points2[i], points3[i], points4[i], computedSurfaces[i], computedCenters[i]);
    }

    // check the results
    for (size_t i = 0; i < points1.size(); ++i) {
        ASSERT_DOUBLE_EQ(computedSurfaces[i].x(), expectedSurfaces[i].x());    
        ASSERT_DOUBLE_EQ(computedSurfaces[i].y(), expectedSurfaces[i].y());
        ASSERT_DOUBLE_EQ(computedSurfaces[i].z(), expectedSurfaces[i].z());

        ASSERT_DOUBLE_EQ(computedCenters[i].x(), expectedCenters[i].x());
        ASSERT_DOUBLE_EQ(computedCenters[i].y(), expectedCenters[i].y());
        ASSERT_DOUBLE_EQ(computedCenters[i].z(), expectedCenters[i].z());
    } 
}

TEST(CMeshTest, TestComputeInterfaces) {
    Config config("input.ini");
    CMesh mesh(config);

    // check the i surfaces
    auto surfacesI = mesh.getSurfacesI();
    int ni = surfacesI.sizeI();
    for (int i=0; i<ni; i++){
        auto versor = surfacesI(i,0,0) / surfacesI(i,0,0).magnitude();
        ASSERT_DOUBLE_EQ(versor.x(), 1.0);
        ASSERT_DOUBLE_EQ(versor.y(), 0.0);
        ASSERT_DOUBLE_EQ(versor.z(), 0.0);
    }

    // check the j surfaces
    auto surfacesJ = mesh.getSurfacesJ();
    int nj = surfacesJ.sizeJ();
    for (int j=0; j<nj; j++){
        auto versor = surfacesJ(0,j,0) / surfacesJ(0,j,0).magnitude();
        ASSERT_DOUBLE_EQ(versor.x(), 0.0);
        ASSERT_DOUBLE_EQ(versor.y(), 1.0);
        ASSERT_DOUBLE_EQ(versor.z(), 0.0);
    }

    // check the k surfaces
    auto surfacesK = mesh.getSurfacesK();
    int nk = surfacesK.sizeK();
    for (int k=0; k<nk; k++){
        auto versor = surfacesK(0,0,k) / surfacesK(0,0,k).magnitude();
        ASSERT_DOUBLE_EQ(versor.x(), 0.0);
        ASSERT_DOUBLE_EQ(versor.y(), 0.0);
        ASSERT_DOUBLE_EQ(versor.z(), 1.0);
    }
}

TEST(CMeshTest, TestComputeBoundaryAreas){
    Config config("input.ini");
    CMesh mesh(config);

    ASSERT_DOUBLE_EQ(mesh.getBoundaryTotalArea(BoundaryIndices::I_START), 5.0);
    ASSERT_DOUBLE_EQ(mesh.getBoundaryTotalArea(BoundaryIndices::I_END), 5.0);
    ASSERT_DOUBLE_EQ(mesh.getBoundaryTotalArea(BoundaryIndices::J_START), 5.0);
    ASSERT_DOUBLE_EQ(mesh.getBoundaryTotalArea(BoundaryIndices::J_END), 5.0);
    ASSERT_DOUBLE_EQ(mesh.getBoundaryTotalArea(BoundaryIndices::K_START), 25.0);
    ASSERT_DOUBLE_EQ(mesh.getBoundaryTotalArea(BoundaryIndices::K_END), 25.0);
}



int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}