#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkXMLStructuredGridWriter.h>

int main() {
    const int dimX = 10;
    const int dimY = 10;
    const int dimZ = 1;

    auto grid = vtkSmartPointer<vtkStructuredGrid>::New();
    grid->SetDimensions(dimX, dimY, dimZ);

    auto points = vtkSmartPointer<vtkPoints>::New();

    // Insert points
    for (int z = 0; z < dimZ; ++z) {
        for (int y = 0; y < dimY; ++y) {
            for (int x = 0; x < dimX; ++x) {
                points->InsertNextPoint(x, y, z);
            }
        }
    }

    grid->SetPoints(points);

    const int numPoints = dimX * dimY * dimZ;

    // Scalar field: Temperature
    auto temperature = vtkSmartPointer<vtkDoubleArray>::New();
    temperature->SetName("Temperature");
    temperature->SetNumberOfComponents(1);
    temperature->SetNumberOfTuples(numPoints);

    // Vector field: Velocity
    auto velocity = vtkSmartPointer<vtkDoubleArray>::New();
    velocity->SetName("Velocity");
    velocity->SetNumberOfComponents(3);  // X, Y, Z
    velocity->SetNumberOfTuples(numPoints);

    for (int i = 0; i < numPoints; ++i) {
        // Dummy data: linear increase
        temperature->SetValue(i, static_cast<double>(i));

        double vx = 0.1 * i;
        double vy = 0.2 * i;
        double vz = 0.0;
        velocity->SetTuple3(i, vx, vy, vz);
    }

    // Add to grid
    grid->GetPointData()->AddArray(temperature);
    grid->GetPointData()->AddArray(velocity);

    // Set active scalars and vectors
    grid->GetPointData()->SetScalars(temperature);
    grid->GetPointData()->SetVectors(velocity);

    // Write to .vts
    auto writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    writer->SetFileName("output.vts");
    writer->SetInputData(grid);
    writer->Write();

    return 0;
}

