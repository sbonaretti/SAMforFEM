/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <FemAssigner.h>

#include <MeshReaderAbaqus.h>
#include <ImageHandler.h>

#include <itkMinimumMaximumImageCalculator.h>

#include <vtkGenericCell.h>
#include <vtkIdList.h>


namespace fem{

	// constructor
	FemAssigner::FemAssigner(){

		_greyLevel = vtkDoubleArray::New();
		_youngModulus = vtkDoubleArray::New();
	}

	// destructor
	FemAssigner::~FemAssigner(){
		
		_greyLevel->Delete();
		_youngModulus->Delete();
	}

	// member functions
	void FemAssigner::GreyLevelAssignmentUpdate(){

		std::cout << std::endl;
		std::cout << "-- GreyLevelAssignmentUpdate --" << std::endl; 
		
		typedef itk::MinimumMaximumImageCalculator< ImageType > MinMaxCalculator;
		MinMaxCalculator::Pointer calculator = MinMaxCalculator::New();
		calculator->SetImage(_image);
		calculator->Compute();
		std::cout << "min intensity in GreyLevelAssignmentUpdate: " << calculator->GetMinimum() << std::endl;
		std::cout << "max intensity in GreyLevelAssignmentUpdate: " << calculator->GetMaximum() << std::endl;
				
		for (int i=0; i<_mesh->GetNumberOfPoints(); i++){

			// get node
			double pt[3];
			_mesh->GetPoint(i, pt);
			
			// round the coordinates to the lower integer and find all the surrounding values
			const ImageType::SpacingType& resolution = _image->GetSpacing();
			const ImageType::PointType& origin = _image->GetOrigin();
			ImageType::IndexType indexA;	
			for( int dim = 0; dim < 3; dim++ )
				indexA[dim] = std::floor((pt[dim] - origin[dim])/resolution[dim]);
			ImageType::IndexType indexB;
			indexB[0] = indexA[0]+1; indexB[1] = indexA[1]; indexB[2] = indexA[2];
			ImageType::IndexType indexC;
			indexC[0] = indexA[0]+1; indexC[1] = indexA[1]+1; indexC[2] = indexA[2];
			ImageType::IndexType indexD;
			indexD[0] = indexA[0]; indexD[1] = indexA[1]+1; indexD[2] = indexA[2];
			ImageType::IndexType indexE;
			indexE[0] = indexA[0]; indexE[1] = indexA[1]; indexE[2] = indexA[2]+1;
			ImageType::IndexType indexF;
			indexF[0] = indexA[0]+1; indexF[1] = indexA[1]; indexF[2] = indexA[2]+1;
			ImageType::IndexType indexG;
			indexG[0] = indexA[0]+1; indexG[1] = indexA[1]+1; indexG[2] = indexA[2]+1;
			ImageType::IndexType indexH;
			indexH[0] = indexA[0]; indexH[1] = indexA[1]+1; indexH[2] = indexA[2]+1;
	
			// get the intensities
			ImageType::PixelType valueA = _image->GetPixel(indexA);
			ImageType::PixelType valueB = _image->GetPixel(indexB);
			ImageType::PixelType valueC = _image->GetPixel(indexC);
			ImageType::PixelType valueD = _image->GetPixel(indexD);
			ImageType::PixelType valueE = _image->GetPixel(indexE);
			ImageType::PixelType valueF = _image->GetPixel(indexF);
			ImageType::PixelType valueG = _image->GetPixel(indexG);
			ImageType::PixelType valueH = _image->GetPixel(indexH);

			// pt in the image coordinate system 
			for( int dim = 0; dim < 3; dim++ )
				pt[dim] = (pt[dim] - origin[dim])/resolution[dim];
			
			// calculate the distances
			double dPtG = std::sqrt((pt[0]-indexG[0])*(pt[0]-indexG[0]) + (pt[1]-indexG[1])*(pt[1]-indexG[1]) + (pt[2]-indexG[2])*(pt[2]-indexG[2])); 
			double dPtH = std::sqrt((pt[0]-indexH[0])*(pt[0]-indexH[0]) + (pt[1]-indexH[1])*(pt[1]-indexH[1]) + (pt[2]-indexH[2])*(pt[2]-indexH[2]));
			double dPtE = std::sqrt((pt[0]-indexE[0])*(pt[0]-indexE[0]) + (pt[1]-indexE[1])*(pt[1]-indexE[1]) + (pt[2]-indexE[2])*(pt[2]-indexE[2]));
			double dPtF = std::sqrt((pt[0]-indexF[0])*(pt[0]-indexF[0]) + (pt[1]-indexF[1])*(pt[1]-indexF[1]) + (pt[2]-indexF[2])*(pt[2]-indexF[2]));
			double dPtC = std::sqrt((pt[0]-indexC[0])*(pt[0]-indexC[0]) + (pt[1]-indexC[1])*(pt[1]-indexC[1]) + (pt[2]-indexC[2])*(pt[2]-indexC[2]));
			double dPtD = std::sqrt((pt[0]-indexD[0])*(pt[0]-indexD[0]) + (pt[1]-indexD[1])*(pt[1]-indexD[1]) + (pt[2]-indexD[2])*(pt[2]-indexD[2]));
			double dPtA = std::sqrt((pt[0]-indexA[0])*(pt[0]-indexA[0]) + (pt[1]-indexA[1])*(pt[1]-indexA[1]) + (pt[2]-indexA[2])*(pt[2]-indexA[2]));
			double dPtB = std::sqrt((pt[0]-indexB[0])*(pt[0]-indexB[0]) + (pt[1]-indexB[1])*(pt[1]-indexB[1]) + (pt[2]-indexB[2])*(pt[2]-indexB[2]));
		
			double intensity = valueA*dPtG + valueB*dPtH + valueC*dPtE + valueD*dPtF +
							   valueE*dPtC + valueF*dPtD + valueG*dPtA + valueH*dPtB;
			intensity /= (dPtG +  dPtH + dPtE + dPtF + dPtC + dPtD + dPtA + dPtB );
			intensity /= 1000; // this is due to convert from g/mm3 to mg/mm3
			
			_greyLevel->InsertNextValue(intensity);

			/*
			std::cout <<  pt[0] << ' ' << pt[1] << ' ' << pt[2] << std::endl;
			std::cout <<  indexA[0] << ' ' << indexA[1] << ' ' << indexA[2] << std::endl;
			std::cout <<  indexB[0] << ' ' << indexB[1] << ' ' << indexB[2] << std::endl;
			std::cout <<  indexC[0] << ' ' << indexC[1] << ' ' << indexC[2] << std::endl;
			std::cout <<  indexD[0] << ' ' << indexD[1] << ' ' << indexD[2] << std::endl;
			std::cout <<  indexE[0] << ' ' << indexE[1] << ' ' << indexE[2] << std::endl;
			std::cout <<  indexF[0] << ' ' << indexF[1] << ' ' << indexF[2] << std::endl;
			std::cout <<  indexG[0] << ' ' << indexG[1] << ' ' << indexG[2] << std::endl;
			std::cout <<  indexH[0] << ' ' << indexH[1] << ' ' << indexH[2] << std::endl;
			
			std::cout << valueA << ' ' << valueB << ' ' << valueC << ' ' << valueD << ' ' << valueE << ' ' << valueF << ' ' << valueG << ' ' << valueH <<std::endl;
			std::cout << intensity <<std::endl;
			*/

		}
	}

}