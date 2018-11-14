/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <RegistrationValidationVolumeImage.h>

#include <itkImageRegionIterator.h>
#include <itkLabelOverlapMeasuresImageFilter.h>

#include <vnl/vnl_vector.h>

#include <vnlWriterVector.h>

using namespace image;
using namespace vnl;

namespace registration{

	// constructor
	RegistrationValidationVolumeImage::RegistrationValidationVolumeImage(){
	}

	// destructor
	RegistrationValidationVolumeImage::~RegistrationValidationVolumeImage(){
	}

	// member function
	void RegistrationValidationVolumeImage::Update(){

		// read reference mask
		ImageHandler* referenceReader = new ImageHandler;
		std::cout << "reference: " << _referenceMaskFileName.toAscii().data() << std::endl;
		referenceReader->SetImageFileName(_referenceMaskFileName);
		referenceReader->MetafileReaderUpdate();

		// dice coefficient
		vnl_vector<double> corticalDiceCoefficient;
		corticalDiceCoefficient.set_size(_movingMaskFileNames.size());
		vnl_vector<double> trabecularDiceCoefficient;
		trabecularDiceCoefficient.set_size(_movingMaskFileNames.size());
		vnl_vector<double> marrowDiceCoefficient;
		marrowDiceCoefficient.set_size(_movingMaskFileNames.size());

		for (int i=0; i<_movingMaskFileNames.size(); i++){
		
			// read warped-to-reference mask
			ImageHandler* warpedToReferenceReader = new ImageHandler;
			std::cout << "warped to reference " << i+1 << ": " << _movingMaskFileNames[i].toAscii().data() << std::endl;		
			warpedToReferenceReader->SetImageFileName(_movingMaskFileNames[i]);
			warpedToReferenceReader->MetafileReaderUpdate();

			ImageHandler::ImageType::PixelType p = itk::NumericTraits< ImageHandler::ImageType::PixelType >::Zero;
			typedef itk::ImageRegionIterator< ImageHandler::ImageType> IteratorType;
			
			// dice coefficient
			typedef itk::LabelOverlapMeasuresImageFilter<ImageHandler::ImageType> overlapMeasureFilter;
			overlapMeasureFilter::Pointer overlapMeasure = overlapMeasureFilter::New();
			double dice = 0.0;
			
			
			/**** CORTICAL LAYER ****/ 
			//extract cortical layer from the reference mask
			ImageHandler::ImageType::Pointer corticalReferenceMask = ImageHandler::ImageType::New();
			corticalReferenceMask->SetRegions( referenceReader->GetImage()->GetRequestedRegion() );
			corticalReferenceMask->CopyInformation( referenceReader->GetImage() );
			corticalReferenceMask->Allocate();
			corticalReferenceMask->FillBuffer(p);
			
			IteratorType it1 (referenceReader->GetImage(), referenceReader->GetImage()->GetRequestedRegion());
			IteratorType it2 (corticalReferenceMask, corticalReferenceMask->GetRequestedRegion());
			
			for (it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd(); ++it1, ++it2){
				if (it1.Get() == 300)
					it2.Set(1);
			}

			// extract cortical layer from the warped-to-reference one
			ImageHandler::ImageType::Pointer corticalWarpedToReferenceMask = ImageHandler::ImageType::New();
			corticalWarpedToReferenceMask->SetRegions( warpedToReferenceReader->GetImage()->GetRequestedRegion() );
			corticalWarpedToReferenceMask->CopyInformation( warpedToReferenceReader->GetImage() );
			corticalWarpedToReferenceMask->Allocate();
			corticalWarpedToReferenceMask->FillBuffer(p);
			
			IteratorType it3 (warpedToReferenceReader->GetImage(), warpedToReferenceReader->GetImage()->GetRequestedRegion());
			IteratorType it4 (corticalWarpedToReferenceMask, corticalWarpedToReferenceMask->GetRequestedRegion());
			
			for (it3.GoToBegin(), it4.GoToBegin(); !it3.IsAtEnd(); ++it3, ++it4){
				if (it3.Get() == 300)
					it4.Set(1);
			}
			
			// dice coefficient
			overlapMeasure->SetSourceImage(corticalReferenceMask);
			overlapMeasure->SetTargetImage(corticalWarpedToReferenceMask);
			overlapMeasure->Update();
			dice = overlapMeasure->GetDiceCoefficient();
			corticalDiceCoefficient(i) = dice;
			//std::cout << "cortical dice coefficient: " << dice << std::endl;

						
			/**** TRABECULAR LAYER ****/ 
			//extract trabecular layer from the reference mask
			ImageHandler::ImageType::Pointer trabecularReferenceMask = ImageHandler::ImageType::New();
			trabecularReferenceMask->SetRegions( referenceReader->GetImage()->GetRequestedRegion() );
			trabecularReferenceMask->CopyInformation( referenceReader->GetImage() );
			trabecularReferenceMask->Allocate();
			trabecularReferenceMask->FillBuffer(p);
			
			IteratorType it5 (referenceReader->GetImage(), referenceReader->GetImage()->GetRequestedRegion());
			IteratorType it6 (trabecularReferenceMask, trabecularReferenceMask->GetRequestedRegion());

			for (it5.GoToBegin(), it6.GoToBegin(); !it5.IsAtEnd(); ++it5, ++it6){
				if (it5.Get() == 200)
					it6.Set(1);
			}
			
			// extract trabecular layer from the warped-to-reference one
			ImageHandler::ImageType::Pointer trabecularWarpedToReferenceMask = ImageHandler::ImageType::New();
			trabecularWarpedToReferenceMask->SetRegions( warpedToReferenceReader->GetImage()->GetRequestedRegion() );
			trabecularWarpedToReferenceMask->CopyInformation( warpedToReferenceReader->GetImage() );
			trabecularWarpedToReferenceMask->Allocate();

			IteratorType it7 (warpedToReferenceReader->GetImage(), warpedToReferenceReader->GetImage()->GetRequestedRegion());
			IteratorType it8 (trabecularWarpedToReferenceMask, trabecularWarpedToReferenceMask->GetRequestedRegion());
			
			trabecularWarpedToReferenceMask->FillBuffer(p);
			for (it7.GoToBegin(), it8.GoToBegin(); !it7.IsAtEnd(); ++it7, ++it8){
				if (it7.Get() == 200)
					it8.Set(1);
			}
			
			// dice coefficient
			overlapMeasure->SetSourceImage(trabecularReferenceMask);
			overlapMeasure->SetTargetImage(trabecularWarpedToReferenceMask);
			overlapMeasure->Update();
			dice = overlapMeasure->GetDiceCoefficient();
			trabecularDiceCoefficient(i) = dice;
			//std::cout << "trabecular dice coefficient: " << dice << std::endl;


			/**** MARROW LAYER ****/ 
			//extract marrow layer from the reference mask
			ImageHandler::ImageType::Pointer marrowReferenceMask = ImageHandler::ImageType::New();
			marrowReferenceMask->SetRegions( referenceReader->GetImage()->GetRequestedRegion() );
			marrowReferenceMask->CopyInformation( referenceReader->GetImage() );
			marrowReferenceMask->Allocate();
			marrowReferenceMask->FillBuffer(p);
						
			IteratorType it9 (referenceReader->GetImage(), referenceReader->GetImage()->GetRequestedRegion());
			IteratorType it10 (marrowReferenceMask, marrowReferenceMask->GetRequestedRegion());

			for (it9.GoToBegin(), it10.GoToBegin(); !it10.IsAtEnd(); ++it9, ++it10){
				if (it9.Get() == 100)
					it10.Set(1);
			}

			// extract marrow layer from warped-to-reference image
			ImageHandler::ImageType::Pointer marrowWarpedToReferenceMask = ImageHandler::ImageType::New();
			marrowWarpedToReferenceMask->SetRegions( warpedToReferenceReader->GetImage()->GetRequestedRegion() );
			marrowWarpedToReferenceMask->CopyInformation( warpedToReferenceReader->GetImage() );
			marrowWarpedToReferenceMask->Allocate();
			marrowWarpedToReferenceMask->FillBuffer(p);
			
			IteratorType it11 (warpedToReferenceReader->GetImage(), warpedToReferenceReader->GetImage()->GetRequestedRegion());
			IteratorType it12 (marrowWarpedToReferenceMask, marrowWarpedToReferenceMask->GetRequestedRegion());
			
			for (it11.GoToBegin(), it12.GoToBegin(); !it11.IsAtEnd(); ++it11, ++it12){
				if (it11.Get() == 100)
					it12.Set(1);
			}
			
			// dice coefficient
			overlapMeasure->SetSourceImage(marrowReferenceMask);
			overlapMeasure->SetTargetImage(marrowWarpedToReferenceMask);
			overlapMeasure->Update();
			dice = overlapMeasure->GetDiceCoefficient();
			marrowDiceCoefficient(i) = dice;
			std::cout << "marrow dice coefficient: " << dice << std::endl;


			delete warpedToReferenceReader;
		}
		
		// std dev
		double corticalStdDev = 0.0;
		for (int i=0; i<corticalDiceCoefficient.size(); i++)
			corticalStdDev += (corticalDiceCoefficient(i) - corticalDiceCoefficient.mean())*(corticalDiceCoefficient(i) - corticalDiceCoefficient.mean());
		corticalStdDev /= (corticalDiceCoefficient.size()-1);
		corticalStdDev = sqrt (corticalStdDev);

		double trabecularStdDev = 0.0;
		for (int i=0; i<trabecularDiceCoefficient.size(); i++)
			trabecularStdDev += (trabecularDiceCoefficient(i) - trabecularDiceCoefficient.mean())*(trabecularDiceCoefficient(i) - trabecularDiceCoefficient.mean());
		trabecularStdDev /= (trabecularDiceCoefficient.size()-1);
		trabecularStdDev = sqrt (trabecularStdDev);

		double marrowStdDev = 0.0;
		for (int i=0; i<marrowDiceCoefficient.size(); i++)
			marrowStdDev += (marrowDiceCoefficient(i) - marrowDiceCoefficient.mean())*(marrowDiceCoefficient(i) - marrowDiceCoefficient.mean());
		marrowStdDev /= (marrowDiceCoefficient.size()-1);
		marrowStdDev = sqrt (marrowStdDev);

		// standard error
		double corticalStdError = corticalStdDev / sqrt(double(corticalDiceCoefficient.size()));
		double trabecularStdError = trabecularStdDev / sqrt(double(trabecularDiceCoefficient.size()));
		double marrowStdError = marrowStdDev / sqrt(double(marrowDiceCoefficient.size()));
	
		std::cout << "cortical dice: average: " << corticalDiceCoefficient.mean() << " std dev: " << corticalStdDev << " std error: " << corticalStdError << std::endl;
		std::cout << "trabecular dice: average: " << trabecularDiceCoefficient.mean() << " std dev: " << trabecularStdDev << " std error: " << trabecularStdError << std::endl;
		std::cout << "marrow dice: average: " << marrowDiceCoefficient.mean() << " std dev: " << marrowStdDev << " std error: " << marrowStdError << std::endl;

		// total - all layers
		double averageDice = 0.0;
		averageDice = (corticalDiceCoefficient.mean() + trabecularDiceCoefficient.mean() + marrowDiceCoefficient.mean()) / 3.0;
		double stdDevDice = 0.0;
		stdDevDice = (corticalDiceCoefficient.mean() - averageDice) * (corticalDiceCoefficient.mean() - averageDice) +
					 (trabecularDiceCoefficient.mean() - averageDice) * (trabecularDiceCoefficient.mean() - averageDice) +
					 (marrowDiceCoefficient.mean() - averageDice) * (marrowDiceCoefficient.mean() - averageDice);
		stdDevDice /= 2.0;
		stdDevDice = sqrt (stdDevDice);
		double stdErrorDice = stdDevDice / sqrt(3.0);

		std::cout << "total dice: average: " << averageDice << " std dev: " << stdDevDice << " std error: " << stdErrorDice << std::endl;
		
		delete referenceReader;
		

	}

}
