/*
 * University of Bern
 * Institute for Surgical Technology and Biomechanics
 * Christof Seiler
 */

#ifndef __itkFieldPCAShapeModelEstimator_txx
#define __itkFieldPCAShapeModelEstimator_txx

#include "itkFieldPCAShapeModelEstimator.h"

namespace itk
{

template<class TInputImage, class TOutputImage>
FieldPCAShapeModelEstimator<TInputImage,TOutputImage>
::FieldPCAShapeModelEstimator(void):m_NumberOfPixels(0),m_NumberOfTrainingImages(0)
{
  //m_EigenVectors.set_size(0,0);
  m_EigenVectors.resize(0);
  m_EigenValues.set_size(0);

  m_NumberOfPrincipalComponentsRequired = 0;
  this->SetNumberOfPrincipalComponentsRequired( 1 );

}

template<class TInputImage, class TOutputImage>
FieldPCAShapeModelEstimator<TInputImage, TOutputImage>
::~FieldPCAShapeModelEstimator(void)
{

}

/**
 * PrintSelf
 */
template <class TInputImage, class TOutputImage>
void
FieldPCAShapeModelEstimator<TInputImage, TOutputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{  

  os << indent << "                   " << std::endl;
  os << indent << "Shape Models " << std::endl;
  os << indent << "Results printed in the superclass " << std::endl;
  os << indent << "                   " << std::endl;

  Superclass::PrintSelf(os,indent);

  itkDebugMacro(<<"                                    ");
  itkDebugMacro(<<"Results of the shape model algorithms");
  itkDebugMacro(<<"====================================");

  itkDebugMacro(<< "The eigen values new method are: ");
    
  itkDebugMacro(<< m_EigenValues);
  itkDebugMacro(<< m_EigenValuesNormalized);

  itkDebugMacro(<< " ");
  itkDebugMacro(<< "==================   ");

  itkDebugMacro(<< "The eigen vectors new method are: ");

  itkDebugMacro(<< "TODO");
  /*for(unsigned int i = 0; i < m_NumberOfPrincipalComponentsRequired; i++) 
    {
    itkDebugMacro(<< m_EigenVectors[i]);
    }  */

  itkDebugMacro(<< " ");
  itkDebugMacro(<< "+++++++++++++++++++++++++");

  // Print out ivars
  os << indent << "NumberOfPrincipalComponentsRequired: ";
  os << m_NumberOfPrincipalComponentsRequired << std::endl;
  os << indent << "NumberOfTrainingImages: ";
  os << m_NumberOfTrainingImages << std::endl;


}// end PrintSelf

/**
 * Enlarge the output requested region to the largest possible region.
 */
template<class TInputImage, class TOutputImage>
void
FieldPCAShapeModelEstimator<TInputImage,TOutputImage>
::EnlargeOutputRequestedRegion( DataObject * itkNotUsed(output) )
{

  // this filter requires the all of the output images to be in
  // the buffer
  for ( unsigned int idx = 0; idx < this->GetNumberOfOutputs(); ++idx)
    {
    if ( this->GetOutput(idx) )
      {
      this->GetOutput(idx)->SetRequestedRegionToLargestPossibleRegion();
      }
    }

}

/**
 * Requires all of the inputs to be in the buffer up to the
 * LargestPossibleRegion of the first input.
 */
template<class TInputImage, class TOutputImage>
void
FieldPCAShapeModelEstimator<TInputImage,TOutputImage>
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();

  if ( this->GetInput(0) )
    {

    // Set the requested region of the first input to largest possible region
    InputImagePointer input = const_cast< TInputImage *>( this->GetInput(0) );
    input->SetRequestedRegionToLargestPossibleRegion();

    // Set the requested region of the remaining input to the largest possible
    // region of the first input
    unsigned int idx;
    for (idx = 1; idx < this->GetNumberOfInputs(); ++idx)
      {
      if ( this->GetInput(idx) )
        {
        typename TInputImage::RegionType requestedRegion = 
          this->GetInput(0)->GetLargestPossibleRegion();

        typename TInputImage::RegionType largestRegion =
          this->GetInput(idx)->GetLargestPossibleRegion();

        if( !largestRegion.IsInside( requestedRegion ) )
          {
          itkExceptionMacro("LargestPossibleRegion of input " << idx <<
            " is not a superset of the LargestPossibleRegion of input 0" );
          }

        InputImagePointer ptr = const_cast<TInputImage *>( this->GetInput(idx) );
        ptr->SetRequestedRegion( requestedRegion );

        } // if ( this->GetIntput(idx))
      } // for idx
    } // if( this->GetInput(0) )

}

/**
 * Generate data (start the model building process)
 */
template<class TInputImage, class TOutputImage>
void 
FieldPCAShapeModelEstimator<TInputImage, TOutputImage>
::GenerateData( )
{
	std::cout << "GenerateData" << std::endl;

  this->EstimateShapeModels();

  // Allocate memory for each output.
  unsigned int numberOfOutputs = 
    static_cast<unsigned int>( this->GetNumberOfOutputs() );

  InputImagePointer input = const_cast<TInputImage *>( this->GetInput(0) );
  unsigned int j;
  for( j = 0; j < numberOfOutputs; j++ )
    {

    OutputImagePointer output = this->GetOutput( j );
    output->SetBufferedRegion( output->GetRequestedRegion() );
    output->Allocate();

    }

  // Fill the output images.
  VectorOfDoubleType oneEigenVector;
  oneEigenVector.set_size(m_EigenVectors.size());
  typedef ImageRegionIterator<OutputImageType>     OutputIterator;

  //Fill the mean image first

  typename OutputImageType::RegionType region = this->GetOutput( 0 )->GetRequestedRegion();
  OutputIterator outIter( this->GetOutput( 0 ), region ); 

  unsigned int i = 0;
  outIter.GoToBegin();
  while( !outIter.IsAtEnd() )
    {
		typename OutputImageType::PixelType value;
		for(unsigned int d = 0; d < InputImageDimension; d++, i++)
			value[d] = m_Means[i];
		outIter.Set(value);
		++outIter;
		//++i;
    } 

  //Now fill the principal component outputs
  //unsigned int kthLargestPrincipalComp = m_NumberOfTrainingImages;
  unsigned int kthLargestPrincipalComp = m_NumberOfPrincipalComponentsRequired;
  unsigned int numberOfValidOutputs = 
    vnl_math_min( numberOfOutputs, m_NumberOfTrainingImages + 1 );

  for( j = 1; j < numberOfValidOutputs; j++ )
    {

    //Extract one column vector at a time
    for(unsigned int i = 0; i < m_EigenVectors.size(); i++)
		oneEigenVector[i] = m_EigenVectors[i][kthLargestPrincipalComp-1];
	//oneEigenVector = m_EigenVectors.get_column(kthLargestPrincipalComp-1);

    region = this->GetOutput( j )->GetRequestedRegion();
    OutputIterator outIterJ( this->GetOutput( j ), region );

    unsigned int idx = 0;
    outIterJ.GoToBegin();
    while( !outIterJ.IsAtEnd() )
      {
		typename OutputImageType::PixelType value;
		for(unsigned int d = 0; d < InputImageDimension; d++, idx++)
			value[d] = oneEigenVector[idx];
		outIterJ.Set(value);
		++outIterJ;
		//++idx;
      }

    //Decrement to get the next principal component
    --kthLargestPrincipalComp;
    }

  // Fill extraneous outputs with zero
  for(; j < numberOfOutputs; j++ )
    {
    region = this->GetOutput( j )->GetRequestedRegion();
    OutputIterator outIterJ( this->GetOutput( j ), region );

    outIterJ.GoToBegin();
    while( !outIterJ.IsAtEnd() )
      {
		typename OutputImageType::PixelType value;
		value.Fill(0.0);
		outIterJ.Set(value);
		++outIterJ;
      }

    }

}// end Generate data


/**
 * Set the number of required principal components
 */
template<class TInputImage, class TOutputImage>
void
FieldPCAShapeModelEstimator<TInputImage,TOutputImage>
::SetNumberOfPrincipalComponentsRequired( unsigned int n )
{
  if ( m_NumberOfPrincipalComponentsRequired != n )
    {
    m_NumberOfPrincipalComponentsRequired = n;

    this->Modified();

    // Modify the required number of outputs ( 1 extra for the mean image ) 
    this->SetNumberOfRequiredOutputs( m_NumberOfPrincipalComponentsRequired + 1 );

    unsigned int numberOfOutputs = static_cast<unsigned int>( this->GetNumberOfOutputs() );
    unsigned int idx;

    if( numberOfOutputs < m_NumberOfPrincipalComponentsRequired + 1 )
      {
      // Make and add extra outputs
      for( idx = numberOfOutputs; idx <= m_NumberOfPrincipalComponentsRequired; idx++ )
        {
        typename DataObject::Pointer output = this->MakeOutput( idx );
        this->SetNthOutput( idx, output.GetPointer() );
        }
      }
    else if ( numberOfOutputs > m_NumberOfPrincipalComponentsRequired + 1 )
      {
      // Remove the extra outputs
      for ( idx = numberOfOutputs - 1; idx >= m_NumberOfPrincipalComponentsRequired + 1; idx-- )
        {
        typename DataObject::Pointer output = this->GetOutputs()[idx];
        this->RemoveOutput( output );
        }
      }

    }
}

/**
 * Set the number of training images.
 */
template<class TInputImage, class TOutputImage>
void
FieldPCAShapeModelEstimator<TInputImage,TOutputImage>
::SetNumberOfTrainingImages( unsigned int n )
{
  if ( m_NumberOfTrainingImages != n )
    {
    m_NumberOfTrainingImages = n;

    this->Modified();

    // Modify the required number of inputs
    this->SetNumberOfRequiredInputs( m_NumberOfTrainingImages );
    }
}

/**-----------------------------------------------------------------
 * Takes a set of training images and returns the means 
 * and variance of the various classes defined in the
 * training set.
 */
template<class TInputImage, class TOutputImage>
void 
FieldPCAShapeModelEstimator<TInputImage, TOutputImage>
::EstimateShapeModels()
{

  this->CalculateInnerProduct();

  this->EstimatePCAShapeModelParameters(); 

}// end EstimateShapeModels 

/**
 * Calculate the inner product between the input training vector
 * where each image is treated as a vector of n-elements
 */
template<class TInputImage, class TOutputImage>
void 
FieldPCAShapeModelEstimator<TInputImage, TOutputImage>
::CalculateInnerProduct()
{
	std::cout << "CalculateInnerProduct" << std::endl;

  // Get the pointers to the input images and initialize the iterators
  // We use dynamic_cast since inputs are stored as DataObjects.  The
  // ImageToImageFilter::GetInput(int) always returns a pointer to a
  // TInputImage1 so it cannot be used for the second input.

  InputImagePointerArray  inputImagePointerArray( m_NumberOfTrainingImages );
  m_InputImageIteratorArray.resize( m_NumberOfTrainingImages );


  for( unsigned int i=0; i< m_NumberOfTrainingImages; i++ )
    {

    InputImageConstPointer inputImagePtr
      = dynamic_cast<const TInputImage*>(ProcessObject::GetInput(i));
    
    inputImagePointerArray[i] = inputImagePtr;

    InputImageConstIterator inputImageIt( inputImagePtr, inputImagePtr->GetBufferedRegion() );

    m_InputImageIteratorArray[i] = inputImageIt;

    m_InputImageIteratorArray[i].GoToBegin();
    
    }

  //-------------------------------------------------------------------
  // Set up the matrix to hold the inner product and the means from the
  // training data
  //-------------------------------------------------------------------
  m_InputImageSize = (inputImagePointerArray[0])->GetBufferedRegion().GetSize();

  m_NumberOfPixels = 1;
  for( unsigned int i = 0; i < InputImageDimension; i++ )
    {
    m_NumberOfPixels *= m_InputImageSize[i];
    }
  m_NumberOfPixels *= InputImageDimension;

  //-------------------------------------------------------------------------
  //Calculate the Means
  //-------------------------------------------------------------------------
  m_Means.set_size(m_NumberOfPixels);
  m_Means.fill(0);

  InputImageConstIterator tempImageItA;

  for(unsigned int img_number = 0; img_number < m_NumberOfTrainingImages; img_number++)
    {
    tempImageItA = m_InputImageIteratorArray[img_number];

    for(unsigned int band_x = 0; band_x < m_NumberOfPixels; /*band_x++*/)
      {
		  for(unsigned int d = 0; d < InputImageDimension; d++, band_x++)
			m_Means[band_x] += tempImageItA.Get()[d];
		  
		  ++tempImageItA;
      }
    } // end: looping through the image
  //-------------------------------------------------------------------------

  m_Means /= m_NumberOfTrainingImages;

  //-------------------------------------------------------------------------
  // Calculate the inner product
  //-------------------------------------------------------------------------
  m_InnerProduct.set_size( m_NumberOfTrainingImages, m_NumberOfTrainingImages );
  m_InnerProduct.fill( 0 );

  InputImageConstIterator tempImageItB; 

  //------------------------------------------------------------------------- 
  for(unsigned int band_x = 0; band_x < m_NumberOfTrainingImages; band_x++)
    {
    for(unsigned int band_y = 0; band_y <= band_x; band_y++ )
      {
      //Pointer to the vector (in original matrix)
      tempImageItA = m_InputImageIteratorArray[band_x];

      //Pointer to the vector in the transposed matrix
      tempImageItB = m_InputImageIteratorArray[band_y];

      for(unsigned int pix_number = 0; pix_number < m_NumberOfPixels; /*pix_number++*/)
        {
			for(unsigned int d = 0; d < InputImageDimension; d++, pix_number++)
				m_InnerProduct[band_x][band_y] += (tempImageItA.Get()[d] - m_Means[pix_number]) * (tempImageItB.Get()[d] - m_Means[pix_number]);

			++tempImageItA; 
			++tempImageItB;
        }// end: looping through the image
      } // end: band_y loop
    } // end: band_x loop
 
  //--------------------------------------------------------------------- 
  // Fill the rest of the inner product matrix and make it symmetric
  //---------------------------------------------------------------------

  for(unsigned int band_x = 0; band_x < (m_NumberOfTrainingImages - 1); band_x++)
    {
    for(unsigned int band_y = band_x+1; band_y < m_NumberOfTrainingImages; band_y++)
      {  
      m_InnerProduct[band_x][band_y] = m_InnerProduct[band_y][band_x];
      }// end band_y loop
    }// end band_x loop

  if( ( m_NumberOfTrainingImages - 1 ) != 0 )
    {
    m_InnerProduct /= ( m_NumberOfTrainingImages - 1 );
    }
  else
    {
    m_InnerProduct.fill(0);
    }  

}// end CalculateInnerProduct


/*-----------------------------------------------------------------
 *Estimage shape models using PCA.
 *-----------------------------------------------------------------
 */
template<class TInputImage, class TOutputImage>
void 
FieldPCAShapeModelEstimator<TInputImage, TOutputImage>
::EstimatePCAShapeModelParameters()
{
  std::cout << "EstimatePCAShapeModelParameters" << std::endl;

  MatrixOfDoubleType identityMatrix( m_NumberOfTrainingImages, m_NumberOfTrainingImages );
  identityMatrix.set_identity();

  vnl_generalized_eigensystem eigenVectors_eigenValues(m_InnerProduct, identityMatrix);

  /*std::cout << "inner product matrix:" << std::endl;
  std::cout << m_InnerProduct << std::endl;*/

  MatrixOfDoubleType eigenVectorsOfInnerProductMatrix = eigenVectors_eigenValues.V;
  
  /*std::cout << "inner product eigenvectors:" << std::endl;
  std::cout << eigenVectors_eigenValues.V << std::endl;*/

  //--------------------------------------------------------------------
  //Calculate the pricipal shape variations
  //
  //m_EigenVectors capture the principal shape variantions
  //m_EigenValues capture the relative weight of each variation
  //Multiply original image vetors with the eigenVectorsOfInnerProductMatrix
  //to derive the pricipal shapes.
  //--------------------------------------------------------------------

  /*//m_EigenVectors.set_size(m_NumberOfPixels, m_NumberOfTrainingImages);
  m_EigenVectors.set_size(m_NumberOfPixels, m_NumberOfPrincipalComponentsRequired);
  m_EigenVectors.fill(0);*/
  m_EigenVectors.resize(m_NumberOfPixels);
  for(unsigned int row = 0; row < m_EigenVectors.size(); row++) {
  	m_EigenVectors[row].resize(m_NumberOfPrincipalComponentsRequired);
	for(unsigned int col = 0; col < m_NumberOfPrincipalComponentsRequired; col++)
		m_EigenVectors[row][col] = 0;
  }

  double pix_value;
  InputImageConstIterator tempImageItA;

  for(unsigned int img_number =0; img_number < m_NumberOfTrainingImages; img_number++ )
    {
    tempImageItA = m_InputImageIteratorArray[img_number];
    for(unsigned int pix_number =0; pix_number < m_NumberOfPixels; pix_number++ )
      {

	  unsigned int index = pix_number%InputImageDimension;
      pix_value = tempImageItA.Get()[index];
      
      unsigned int shift = m_NumberOfTrainingImages-m_NumberOfPrincipalComponentsRequired;
	  //for( unsigned int vec_number = 0; vec_number < m_NumberOfTrainingImages; vec_number++ )
	  for( unsigned int vec_number = shift; vec_number < m_NumberOfTrainingImages; vec_number++ )
        {
        m_EigenVectors[pix_number][vec_number-shift] += (pix_value * eigenVectorsOfInnerProductMatrix[img_number][vec_number]);
		//m_EigenVectors[pix_number][vec_number] += (pix_value * eigenVectorsOfInnerProductMatrix[img_number][vec_number]);
        }

	  if(index == InputImageDimension-1)
	      ++tempImageItA;

      }
    }

  VectorOfDoubleType oneEigenVector(m_NumberOfPixels);
  oneEigenVector.fill(0);
  for(unsigned int col = 0; col < m_NumberOfPrincipalComponentsRequired; col++) {
  	for(unsigned int row = 0; row < oneEigenVector.size(); row++)
  		oneEigenVector[row] = m_EigenVectors[row][col];
  	oneEigenVector.normalize();
	for(unsigned int row = 0; row < oneEigenVector.size(); row++)
  		m_EigenVectors[row][col] = oneEigenVector[row];
  }
  //m_EigenVectors.normalize_columns();
  
  m_EigenValues.set_size(m_NumberOfTrainingImages);

  //Extract the diagonal elements into the Eigen value vector
  m_EigenValues = (eigenVectors_eigenValues.D).diagonal();
  
  //Flip the eigen values since the eigen vectors output
  //is ordered in decending order of theri corresponding eigen values.
  m_EigenValues.flip();

  //--------------------------------------------------------------------
  //Normalize the eigen values 
  //--------------------------------------------------------------------

  m_EigenValuesNormalized = m_EigenValues;
  m_EigenValuesNormalized.normalize();
  m_EigenValuesNormalized = element_product(m_EigenValuesNormalized, m_EigenValuesNormalized);
  

}// end EstimatePCAShapeModelParameters

//-----------------------------------------------------------------


} // namespace itk

#endif
