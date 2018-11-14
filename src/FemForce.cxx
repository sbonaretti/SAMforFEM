/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <FemForce.h>

#include <iostream>
#include <math.h>

#include <MeshExtractOuterSurface.h>

#include <itkExponentialDeformationFieldImageFilter.h>
#include <itkMultiplyByConstantImageFilter.h>

using namespace mesh;

namespace fem{

	// constructor
	FemForce::FemForce(){

		
	}

	// destructor
	FemForce::~FemForce(){
	}

	
	// member functions
	void FemForce::WalkingFromFemulLenght(){

		// extract points
		double P0[3]; double P1[3]; double P2[3]; double F1[3]; double F3[3]; double F5[3];
		if (_points->GetNumberOfPoints() == 6){ // the reference file contains 7 points; the output files contains 10. See header file
			_points->GetPoint(0, P0); 
			_points->GetPoint(1, P1); 
			_points->GetPoint(2, P2); 
			_points->GetPoint(3, F1); 
			_points->GetPoint(4, F3); 
			_points->GetPoint(5, F5); 
		}
		else std::cout << "wrong number of force and bc application points" << std::endl;

		//std::cout << point1[0] << ' ' << point1[1] << ' ' << point1[2] << ' ' << std::endl;
				
		// calculation of the femur length (F1-F3)
		double femurLength = sqrt( (F1[0]-F3[0]) * (F1[0]-F3[0]) + 
								   (F1[1]-F3[1]) * (F1[1]-F3[1]) + 
								   (F1[2]-F3[2]) * (F1[2]-F3[2]) 
								 );
		femurLength /= 1000; //from mm to m
		std::cout << "femur length: " << femurLength << std::endl;

		// calculation of the body height
		double bodyHeight = femurLength * 100 / 26.75;
		std::cout << "body height: " << bodyHeight << std::endl;

		// calculation of the BMI factor for female
		//std::cout << "BMI for FEMALE: ";
		//double max = 26.36 + 1.5*6.01;
		//double min = 26.36 - 1.5*6.01;
		
		// calculation of the BMI factor for female
		std::cout << "BMI for MALE: ";
		double max = 25.49 + 1.5*4.90;
		double min = 25.49 - 1.5*4.90;
				
		double BMI = fmod(rand(), (max-min)) + min; 
		std::cout << BMI << std::endl;

		// calculation of the body weight
		double bodyWeight = BMI * bodyHeight*bodyHeight;
		std::cout << "body weight: " << bodyWeight << std::endl;
		
		bodyWeight *= 9.8; // gravity

		_points->SetNumberOfPoints(10);

		// calculation of hip contact force
		double x = -54.0 / 100 * bodyWeight; // /100 is due to the %
		double y = -32.8 / 100 * bodyWeight;
		double z = -229.2 / 100 * bodyWeight;
		double hipContact[3];
		hipContact[0] = std::sqrt(x*x + y*y + z*z);
		hipContact[1] = 0.0;
		hipContact[2] = 0.0;
		_points->SetPoint(0, hipContact);
		std::cout << "hipContact: " << hipContact[0] << ' ' << hipContact[1] << ' ' << hipContact[2] << std::endl;

		
		// calculation of the abductor force
		double abductor[3];
		abductor[0] = 58.0 / 100 * bodyWeight;
		abductor[1] = 4.3 / 100 * bodyWeight;
		abductor[2] = 86.5 / 100 * bodyWeight;
		_points->SetPoint(1, abductor);
		std::cout << "abductor: " << abductor[0] << ' ' << abductor[1] << ' ' << abductor[2] << std::endl;


		// calculation of the proximal tensor fascia latae force
		double proxTensor[3];
		proxTensor[0] = 7.2 / 100 * bodyWeight;
		proxTensor[1] = 11.6 / 100 * bodyWeight;
		proxTensor[2] = 13.2 / 100 * bodyWeight;
		_points->SetPoint(2, proxTensor);
		std::cout << "proxTensor: " << proxTensor[0] << ' ' << proxTensor[1] << ' ' << proxTensor[2] << std::endl;


		// calculation of the distal tensor fascia latae force
		double distTensor[3];
		distTensor[0] = -0.5 / 100 * bodyWeight;
		distTensor[1] = -0.7 / 100 * bodyWeight;
		distTensor[2] = -19.0 / 100 * bodyWeight;
		_points->SetPoint(3, distTensor);
		std::cout << "distTensor: " << distTensor[0] << ' ' << distTensor[1] << ' ' << distTensor[2] << std::endl;


		// calculation of the vastus lateralis force
		double vastus[3];
		vastus[0] = -0.9 / 100 * bodyWeight;
		vastus[1] = 18.5 / 100 * bodyWeight;
		vastus[2] = -92.9 / 100 * bodyWeight;
		_points->SetPoint(4, vastus);
		std::cout << "vastus: " << vastus[0] << ' ' << vastus[1] << ' ' << vastus[2] << std::endl;


		// new coordinate system unit vectors for abaqus for the BOUNDARY CONDITIONS
		// A = F3 - P0
		double diff[3];
		diff[0] = F3[0] - P0[0];
		diff[1] = F3[1] - P0[1];
		diff[2] = F3[2] - P0[2];
		double norm = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
		double A[3];
		A[0] = diff[0] / norm;
		A[1] = diff[1] / norm;
		A[2] = diff[2] / norm;
		_points->SetPoint(5, A);
		std::cout << "A: " << A[0] << ' ' << A[1] << ' ' << A[2] << std::endl;

		// B = F5 - P0
		diff[0] = F5[0] - P0[0];
		diff[1] = F5[1] - P0[1];
		diff[2] = F5[2] - P0[2];
		norm = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
		double B[3];
		B[0] = diff[0] / norm;
		B[1] = diff[1] / norm;
		B[2] = diff[2] / norm;
		_points->SetPoint(6, B);
		std::cout << "B: " << B[0] << ' ' << B[1] << ' ' << B[2] << std::endl;


		// new coordinate system unit vectors for abaqus for the FORCES
		// C = F3 - F1
		diff[0] = F3[0] - F1[0];
		diff[1] = F3[1] - F1[1];
		diff[2] = F3[2] - F1[2];
		norm = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
		double C[3];
		C[0] = diff[0] / norm;
		C[1] = diff[1] / norm;
		C[2] = diff[2] / norm;
		_points->SetPoint(7, C);
		std::cout << "C: " << C[0] << ' ' << C[1] << ' ' << C[2] << std::endl;

		// D = F5 - F1
		diff[0] = F5[0] - F1[0];
		diff[1] = F5[1] - F1[1];
		diff[2] = F5[2] - F1[2];
		norm = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
		double D[3];
		D[0] = diff[0] / norm;
		D[1] = diff[1] / norm;
		D[2] = diff[2] / norm;
		_points->SetPoint(8, D);
		std::cout << "D: " << D[0] << ' ' << D[1] << ' ' << D[2] << std::endl;

		// bone data
		double temp[3];
		temp[0] = femurLength;
		temp[1] = bodyHeight;
		temp[2] = BMI;
		_points->SetPoint(9, temp);


		// application nodes id vector 
		double P0id = _idVector->GetValue(0);
		double P1id = _idVector->GetValue(1);
		double P2id = _idVector->GetValue(2);
		double F1id = _idVector->GetValue(3);
		double F3id = _idVector->GetValue(4);
		double F5id = _idVector->GetValue(5);
		
		_idVector->SetNumberOfValues(10);
		
		_idVector->SetValue(0,P0id);
		_idVector->SetValue(1,P1id);
		_idVector->SetValue(2,P1id);
		_idVector->SetValue(3,P1id);
		_idVector->SetValue(4,P2id);
		_idVector->SetValue(5,0); // no application point for A
		_idVector->SetValue(6,0); // no application point for B
		_idVector->SetValue(7,0); // no application point for C
		_idVector->SetValue(8,0); // no application point for D
		bodyWeight /= 10;
		_idVector->SetValue(9,bodyWeight); 
 	
	}

	void FemForce::SideFalling(){

		// extract points
		double P0[3]; double P1[3]; double P2[3]; double F1[3]; double F3[3]; double F5[3];
		if (_points->GetNumberOfPoints() == 6){ // the reference file contains 7 points; the output files contains 10. See header file
			_points->GetPoint(0, P0); 
			_points->GetPoint(1, P1); 
			_points->GetPoint(2, P2); 
			_points->GetPoint(3, F1); 
			_points->GetPoint(4, F3); 
			_points->GetPoint(5, F5); 
		}
		else std::cout << "wrong number of force and bc application points" << std::endl;

		_points->SetNumberOfPoints(3);
		
		// calculation of the body weight
		double force[3];
		force[0] = _bodyWeight * 2.5 * 0.1736 * 9.8; // 2.5 = 2.5 times bodyweight; 0.1736 = cos(80); 9.8 = gravity
		force[1] = _bodyWeight * 2.5 * 0.9659 * 9.8; // 0.9659 = cos(15)
		force[2] = 0.0;
		_points->SetPoint(0, force);
		std::cout << "force: " << force[0] << ' ' << force[1] << ' ' << force[2] << std::endl;

		// new coordinate system unit vectors for abaqus for the FORCES
		// C = F3 - F1
		double diff[3];
		diff[0] = F3[0] - F1[0];
		diff[1] = F3[1] - F1[1];
		diff[2] = F3[2] - F1[2];
		double norm = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
		double C[3];
		C[0] = diff[0] / norm;
		C[1] = diff[1] / norm;
		C[2] = diff[2] / norm;
		_points->SetPoint(1, C);
		std::cout << "C: " << C[0] << ' ' << C[1] << ' ' << C[2] << std::endl;
	
		// F = P0 - F1
		diff[0] = P0[0] - F1[0];
		diff[1] = P0[1] - F1[1];
		diff[2] = P0[2] - F1[2];
		norm = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
		double E[3];
		E[0] = diff[0] / norm;
		E[1] = diff[1] / norm;
		E[2] = diff[2] / norm;
		_points->SetPoint(2, E);
		std::cout << "E: " << E[0] << ' ' << E[1] << ' ' << E[2] << std::endl;
	
		// id vector
		_idVector->SetNumberOfValues(3);
		double P0id = _idVector->GetValue(0);
		_idVector->SetValue(0,P0id);
		_idVector->SetValue(1,0);
		_idVector->SetValue(2,0);


	}
	
	void FemForce::Standing() {

		// for comparison paper - fem validation

		// extract points
		double P0[3]; double P1[3]; double P2[3]; double F1[3]; double F3[3]; double F5[3];
		if (_points->GetNumberOfPoints() == 6){ // the reference file contains 7 points; the output files contains 10. See header file
			_points->GetPoint(0, P0); 
			_points->GetPoint(1, P1); 
			_points->GetPoint(2, P2); 
			_points->GetPoint(3, F1); 
			_points->GetPoint(4, F3); 
			_points->GetPoint(5, F5); 
		}
		else std::cout << "wrong number of force and bc application points" << std::endl;

		_points->SetNumberOfPoints(5);
		
		// calculation of the hip contact force (see fig 3.1 on thesis)
		double bodyWeight = 100 * 9.8; // 9.8 is gravity 
		double x = -54.0 / 100 * bodyWeight; 
		double y = -32.8 / 100 * bodyWeight;
		double z = -229.2 / 100 * bodyWeight;
		double hipContact[3];
		hipContact[0] = std::sqrt(x*x + y*y + z*z);
		hipContact[1] = 0.0;
		hipContact[2] = 0.0;
		_points->SetPoint(0, hipContact);
		std::cout << "hipContact: " << hipContact[0] << ' ' << hipContact[1] << ' ' << hipContact[2] << std::endl;

		// new coordinate system unit vectors for abaqus for the BOUNDARY CONDITIONS
		// A = F3 - P0
		double diff[3];
		diff[0] = F3[0] - P0[0];
		diff[1] = F3[1] - P0[1];
		diff[2] = F3[2] - P0[2];
		double norm = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
		double A[3];
		A[0] = diff[0] / norm;
		A[1] = diff[1] / norm;
		A[2] = diff[2] / norm;
		_points->SetPoint(1, A);
		std::cout << "A: " << A[0] << ' ' << A[1] << ' ' << A[2] << std::endl;

		// B = F5 - P0
		diff[0] = F5[0] - P0[0];
		diff[1] = F5[1] - P0[1];
		diff[2] = F5[2] - P0[2];
		norm = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
		double B[3];
		B[0] = diff[0] / norm;
		B[1] = diff[1] / norm;
		B[2] = diff[2] / norm;
		_points->SetPoint(2, B);
		std::cout << "B: " << B[0] << ' ' << B[1] << ' ' << B[2] << std::endl;


		// new coordinate system unit vectors for abaqus for the FORCES
		// C = F3 - F1
		diff[0] = F3[0] - F1[0];
		diff[1] = F3[1] - F1[1];
		diff[2] = F3[2] - F1[2];
		norm = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
		double C[3];
		C[0] = diff[0] / norm;
		C[1] = diff[1] / norm;
		C[2] = diff[2] / norm;
		_points->SetPoint(3, C);
		std::cout << "C: " << C[0] << ' ' << C[1] << ' ' << C[2] << std::endl;

		// D = F5 - F1
		diff[0] = F5[0] - F1[0];
		diff[1] = F5[1] - F1[1];
		diff[2] = F5[2] - F1[2];
		norm = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
		double D[3];
		D[0] = diff[0] / norm;
		D[1] = diff[1] / norm;
		D[2] = diff[2] / norm;
		_points->SetPoint(4, D);
		std::cout << "D: " << D[0] << ' ' << D[1] << ' ' << D[2] << std::endl;

		// id vector
		_idVector->SetNumberOfValues(5);
		double P0id = _idVector->GetValue(0);
		_idVector->SetValue(0,P0id);
		_idVector->SetValue(1,0);
		_idVector->SetValue(2,0);
		_idVector->SetValue(3,0);
		_idVector->SetValue(4,0);
		/*

		// new coordinate system unit vectors for abaqus for the BOUNDARY CONDITIONS
		// A = F3 - P0
		double diff[3];
		diff[0] = F3[0] - P0[0];
		diff[1] = F3[1] - P0[1];
		diff[2] = F3[2] - P0[2];
		double norm = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
		double A[3];
		A[0] = diff[0] / norm;
		A[1] = diff[1] / norm;
		A[2] = diff[2] / norm;
		_points->SetPoint(1, A);
		std::cout << "A: " << A[0] << ' ' << A[1] << ' ' << A[2] << std::endl;

		// B = F5 - P0
		diff[0] = F5[0] - P0[0];
		diff[1] = F5[1] - P0[1];
		diff[2] = F5[2] - P0[2];
		norm = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
		double B[3];
		B[0] = diff[0] / norm;
		B[1] = diff[1] / norm;
		B[2] = diff[2] / norm;
		_points->SetPoint(2, B);
		std::cout << "B: " << B[0] << ' ' << B[1] << ' ' << B[2] << std::endl;
		
		// calculation of the force direction
		double a[3];
		a[0] = F1[0] - F3[0];
		a[1] = F1[1] - F3[1];
		a[2] = F1[2] - F3[2];

		double b[3];
		b[0] = P0[0] - F3[0];
		b[1] = P0[1] - F3[1];
		b[2] = P0[2] - F3[2];

		double norm_a = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
		double norm_b = sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);

		double cos_alfa = (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]) / (norm_a*norm_b);
		double sin_alfa = sqrt(1 - (cos_alfa*cos_alfa));
		double alfa = acos(cos_alfa)*(180/3.14);

		std::cout << "cos_alfa: " << cos_alfa << std::endl;
		std::cout << "sin_alfa: " << sin_alfa << std::endl;
		std::cout << "alfa: " << alfa << std::endl;
		
		// calculation of the magnitude
		double force[3];
		force[0] = -100 * sin_alfa * 9.8; // 100 = weight, 9.8 = gravity
		force[1] = 0.0; 
		force[2] = -100 * cos_alfa * 9.8;
		_points->SetPoint(0, force);
		std::cout << "force: " << force[0] << ' ' << force[1] << ' ' << force[2] << std::endl;
		
		// new coordinate system unit vectors for abaqus for the FORCES
		// C = F3 - F1
		diff[3];
		diff[0] = F3[0] - F1[0];
		diff[1] = F3[1] - F1[1];
		diff[2] = F3[2] - F1[2];
		norm = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
		double C[3];
		C[0] = diff[0] / norm;
		C[1] = diff[1] / norm;
		C[2] = diff[2] / norm;
		_points->SetPoint(4, C);
		std::cout << "C: " << C[0] << ' ' << C[1] << ' ' << C[2] << std::endl;
	
		// F = P0 - F1
		diff[0] = P0[0] - F1[0];
		diff[1] = P0[1] - F1[1];
		diff[2] = P0[2] - F1[2];
		norm = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
		double E[3];
		E[0] = diff[0] / norm;
		E[1] = diff[1] / norm;
		E[2] = diff[2] / norm;
		_points->SetPoint(5, E);
		std::cout << "E: " << E[0] << ' ' << E[1] << ' ' << E[2] << std::endl;
		*/
			
		

	
	}
	
	void FemForce::PropagateApplicationPoints(){
		
		_idVector = vtkDoubleArray::New();

		// extract surface mesh (points applied on the bone surface)
		MeshExtractOuterSurface* surfaceExtractor = new MeshExtractOuterSurface;
		surfaceExtractor->SetVolumeMesh(_mesh);
		surfaceExtractor->Update();

		// invert the SVF
		
		std::cout << "inverting the SVF" << std::endl; 
		typedef itk::MultiplyByConstantImageFilter< FieldType, double, FieldType > MultiplyFilter;
		MultiplyFilter::Pointer multiplyFilter = MultiplyFilter::New();
		multiplyFilter->SetInput(_field);
		multiplyFilter->SetConstant(-1); // -1 for the reconstructed, +1 for the originals
		multiplyFilter->Update();
		
		// VF to DVF 
		std::cout << "DVF creation (exponentiator)" << std::endl;
		typedef itk::ExponentialDeformationFieldImageFilter< FieldType,FieldType > ExponentialFieldFilterType;
		ExponentialFieldFilterType::Pointer exponentiator = ExponentialFieldFilterType::New();
		exponentiator->SetInput(multiplyFilter->GetOutput());
		exponentiator->Update();
		
		// for each point
		for (int i=0; i<_points->GetNumberOfPoints(); i++){
			
			//std::cout << std::endl;

			// get the point
			double coord[3];
			_points->GetPoint(i, coord);
			std::cout << "reference point" << std::endl;
			//std::cout << coord[0] << ' ' << coord[1] << ' ' << coord[2] << std::endl;

			// divide it by the resolution
			const FieldType::SpacingType& spacing = _field->GetSpacing();
			coord[0] /= spacing[0];
			coord[1] /= spacing[1];
			coord[2] /= spacing[2];
			//std::cout << "reference point divided by the resolution" << std::endl;
			//std::cout << coord[0] << ' ' << coord[1] << ' ' << coord[2] << std::endl;

			// pick the point in the deformation vector field
			FieldType::IndexType pixelIndex;
			pixelIndex[0] = floor (coord[0] + 0.5); 
			pixelIndex[1] = floor (coord[1] + 0.5); 
			pixelIndex[2] = floor (coord[2] + 0.5); 
			//std::cout << "pixel index" << std::endl;
			//std::cout << pixelIndex[0] << ' ' << pixelIndex[1] << ' ' << pixelIndex[2] << std::endl;
			
			// for the reconstructed bones with image-based
			FieldType::PixelType pixelValue = exponentiator->GetOutput()->GetPixel( pixelIndex );
			
			//FieldType::PixelType pixelValue = _field->GetPixel( pixelIndex );
			

			// propagated point
			coord[0] += pixelValue[0];
			coord[1] += pixelValue[1];
			coord[2] += pixelValue[2];
			//std::cout << "pixel value" << std::endl;
			//std::cout << pixelValue[0] << ' ' << pixelValue[1] << ' ' << pixelValue[2] << std::endl;
			//std::cout << "propagated point" << std::endl;
			//std::cout << coord[0] << ' ' << coord[1] << ' ' << coord[2] << std::endl;

			// remultiply it by the spacing
			coord[0] *= spacing[0];
			coord[1] *= spacing[1];
			coord[2] *= spacing[2];
			//std::cout << "remultiplied by the spacing" << std::endl;
			//std::cout << coord[0] << ' ' << coord[1] << ' ' << coord[2] << std::endl;
			
			// find the closest node
			int pointIndex;
			double closestX; double closestY; double closestZ;
			findClosestNode (surfaceExtractor->GetSurfaceMesh(), coord, pointIndex, closestX, closestY, closestZ);
			std::cout << "closest point" << std::endl;
			coord[0] = closestX;
			coord[1] = closestY;
			coord[2] = closestZ;
			//std::cout << closestX << ' ' << closestY << ' ' << closestZ << std::endl;
			//std::cout << pointIndex << std::endl;
			
			// put it back in vtkPoints
			_points->InsertPoint(i,coord);
			// put in the id vector
			_idVector->InsertValue(i,pointIndex);

		}
		
		std::cout << "final points" << std::endl;
		for (int i=0; i<_points->GetNumberOfPoints(); i++){
			double coord[3];
			_points->GetPoint(i, coord);
			
			//std::cout << coord[0] << ' ' << coord[1] << ' ' << coord[2] << std::endl;
		}
	}


	// private
	void FemForce::findClosestNode (vtkPolyData* mesh, double point[3], int &minIndex, double &closestX, double &closestY, double &closestZ){
	
		// find the absolute distances between the input point and all the mesh nodes
		double min = 0.0; //int minIndex;
		int nOfCoord = 3;
		double closestPoint[3];
		
		for (int i=0; i<mesh->GetNumberOfPoints(); i++){
			double coord[3];
			mesh->GetPoint(i, coord);
			//std::cout << std::endl;
			//std::cout << "point: " << point[0] << ' ' << point[1] << ' ' << point[2] << std::endl;
			//std::cout << "mesh point: " << coord[0] << ' ' << coord[1] << ' ' << coord[2] << std::endl;
			

			double distance = 0.0;
			for (int a=0; a<nOfCoord; a++)
				distance += ((coord[a]-(point[a])) * (coord[a]-(point[a])));
			distance = sqrt(distance);

			//std::cout << "distance: " << distance << std::endl;

			if (i==0){
				min = distance;
				minIndex=i+1; // in abaqus node id starts from 1
				closestPoint[0] = coord[0];
				closestPoint[1] = coord[1];
				closestPoint[2] = coord[2];
				//std::cout << "temporary closest point: " << closestPoint[0] << ' ' << closestPoint[1] << ' ' << closestPoint[2] << std::endl;
			}
			else{
				if (distance < min){
					min = distance;
					minIndex=i+1; // in abaqus node id starts from 1
					closestPoint[0] = coord[0];
					closestPoint[1] = coord[1];
					closestPoint[2] = coord[2];
					//std::cout << "temporary closest point: " << closestPoint[0] << ' ' << closestPoint[1] << ' ' << closestPoint[2] << std::endl;
				}
			}
			
		}
		std::cout << closestPoint[0] << ' ' << closestPoint[1] << ' ' << closestPoint[2] << std::endl;
		closestX = closestPoint[0];
		closestY = closestPoint[1];
		closestZ = closestPoint[2];


	}

}