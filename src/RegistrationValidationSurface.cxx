/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <RegistrationValidationSurface.h>
#include <PointWriterXyz.h>

#include <vtkCellLocator.h>
#include <vtkDoubleArray.h>


using namespace points;

namespace registration{

	// constructor
	RegistrationValidationSurface::RegistrationValidationSurface(){

		_HausdorffDistances = vtkDoubleArray::New();
		_average = 0.0;
		_standardDeviation = 0.0; 
		_standardError = 0.0;
	}

	// destructor
	RegistrationValidationSurface::~RegistrationValidationSurface(){

	_HausdorffDistances->Delete();
	}



	// private function
	void RegistrationValidationSurface::findPointToSurfaceDistances(vtkPoints* points, vtkPolyData* surface, double & maxDistance){

		// cell locator
		vtkCellLocator* cellLocator = vtkCellLocator::New();
		cellLocator->SetDataSet(surface);
		cellLocator->SetNumberOfCellsPerBucket(1);
		cellLocator->BuildLocator();

		// find closest points
		vtkPoints* closestPoints = vtkPoints::New();

		for (int a=0; a<points->GetNumberOfPoints(); a++){
				
				double point[3]; double surfacePoint[3];
				vtkIdType cellId; int subId; double dist;
				
				points->GetPoint(a,point);
				cellLocator->FindClosestPoint(point, surfacePoint, cellId, subId, dist);

				closestPoints->InsertNextPoint(surfacePoint);
		}

		// write points
		PointWriterXyz* pointWriter = new PointWriterXyz;
		pointWriter->SetInput(points);
		pointWriter->SetFileName("extractedPoints.txt");
		pointWriter->Update();
		// write closest point
		pointWriter->SetInput(closestPoints);
		pointWriter->SetFileName("closestPoints.txt");
		pointWriter->Update();
		
		
		// calculate distances between points and the closest ones on the surface
		vtkDoubleArray* differences = vtkDoubleArray::New();
		
		for (int a=0; a<points->GetNumberOfPoints(); a++){
			
			double s[3];
			points->GetPoint(a, s);
			double c[3];
			closestPoints->GetPoint(a, c);
		
			differences->InsertNextValue((sqrt( (s[0]-c[0])*(s[0]-c[0]) +  (s[1]-c[1])*(s[1]-c[1]) + (s[2]-c[2])*(s[2]-c[2]) )));
		}

		// find the maximum distance
		double range[2];
		differences->GetRange(range);
		maxDistance = range[1];

		// cleaning up
		cellLocator->Delete();
		differences->Delete();
		delete pointWriter;
		// differences->Delete(); // it crashes
	}
}