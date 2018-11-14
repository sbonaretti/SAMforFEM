/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <RegistrationValidarionVolumeMesh.h>

#include <ImageHandler.h>
#include <MeshReaderMorpherVolume.h>
#include <PointWriterXyz.h>

#include <itkImageRegionIterator.h>

#include <vtkCellLocator.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkSTLWriter.h>


using namespace image;
using namespace mesh;
using namespace points;

namespace mesh{

	// constructor
	RegistrationValidarionVolumeMesh::RegistrationValidarionVolumeMesh(){
	}

	// destructor
	RegistrationValidarionVolumeMesh::~RegistrationValidarionVolumeMesh(){
	}

	// member function
	void RegistrationValidarionVolumeMesh::Update(){

	// REFERENCE //
		// load the reference mesh
		MeshReaderMorpherVolume* referenceMeshReader = new MeshReaderMorpherVolume;
		referenceMeshReader->SetFileName(_referenceMeshFileName);
		std::cout << "reference mesh: " << _referenceMeshFileName.toAscii().data() << std::endl;
		referenceMeshReader->Update();

		// load the reference mask
		ImageHandler *imageHandler = new ImageHandler;
		imageHandler->SetImageFileName(_referenceMaskFileName);
		std::cout << "reference mesh: " << _referenceMaskFileName.toAscii().data() << std::endl; 
		imageHandler->MetafileReaderUpdate();

		// label node IDs in the three parts (cortical, trabecular, marrow)
		vtkIntArray* cortical = vtkIntArray::New();
		vtkIntArray* trabecular = vtkIntArray::New();
		vtkIntArray* marrow = vtkIntArray::New();
		
		// nodes position in the reference image
		ReferenceLabels(referenceMeshReader->GetOutput(), imageHandler->GetImage(), cortical, trabecular, marrow);

		// cleaning up
		delete imageHandler;
		delete referenceMeshReader;
		
	// MOVINGS //
		// evaluate if in the other meshes, nodes are in the same areas
		vtkDoubleArray* corticalDistance = vtkDoubleArray::New();
		vtkDoubleArray* trabecularDistance = vtkDoubleArray::New();
		vtkDoubleArray* marrowDistance = vtkDoubleArray::New();

		vtkIntArray* nOfRightCortical = vtkIntArray::New();
		vtkIntArray* nOfRightTrabecular = vtkIntArray::New();
		vtkIntArray* nOfRightMarrow = vtkIntArray::New();

		for (int i=0; i<_movingMeshFileNames.size(); i++){
			
			// load mesh
			std::cout << std::endl;
			std::cout << "mesh: " << _movingMeshFileNames[i].ascii() << std::endl;
			MeshReaderMorpherVolume* meshReader = new MeshReaderMorpherVolume;
			meshReader->SetFileName(_movingMeshFileNames[i]);
			meshReader->Update();

			// load correspondent image
			QString temp = _movingMeshFileNames[i];
			temp.remove(temp.lastIndexOf("/")+6,temp.size()-1);
			temp.remove(0, temp.lastIndexOf("/")+1);
			ImageHandler* imageHandler = new ImageHandler;
			for (int y=0; y<_movingMeshFileNames.size(); y++){
					if (_movingMaskFileNames[y].contains(temp)){
						std::cout << "mask: " << _movingMaskFileNames[y].ascii() << std::endl;
						imageHandler->SetImageFileName(_movingMaskFileNames[y]);
						imageHandler->MetafileReaderUpdate();
				}
			}
			double distance = 0.0; int nOfRight = 0.0;

			checkMovings(cortical, meshReader->GetOutput(), imageHandler->GetImage(), 1, distance, nOfRight);
			corticalDistance->InsertNextValue(distance); nOfRightCortical->InsertNextValue(nOfRight);
			distance = 0.0; nOfRight = 0.0;
			
			checkMovings(trabecular, meshReader->GetOutput(), imageHandler->GetImage(), 2, distance, nOfRight);
			trabecularDistance->InsertNextValue(distance); nOfRightTrabecular->InsertNextValue(nOfRight);
			distance = 0.0; nOfRight = 0.0;
			
			checkMovings(marrow, meshReader->GetOutput(), imageHandler->GetImage(), 3, distance, nOfRight);
			marrowDistance->InsertNextValue(distance); nOfRightMarrow->InsertNextValue(nOfRight);
			distance = 0.0; nOfRight = 0.0;

			// cleaning up
			delete meshReader;
			delete imageHandler;
		}

	// STATISTICS //
		// measurements - cortical
		// right proportion
		// average
		std::cout << std::endl;
		double average = 0.0;
		for (int a=0; a<nOfRightCortical->GetNumberOfTuples(); a++)
			average += (nOfRightCortical->GetValue(a) * 100 / cortical->GetNumberOfTuples());
		average /= nOfRightCortical->GetNumberOfTuples();
		std::cout << "right cortical average: " << average << "%" << std::endl;
		// standard deviation
		double standardDeviation = 0.0;
		for (int a=0; a<nOfRightCortical->GetNumberOfTuples(); a++)
			standardDeviation += (nOfRightCortical->GetValue(a) * 100 / cortical->GetNumberOfTuples() - average)*(nOfRightCortical->GetValue(a) * 100 / cortical->GetNumberOfTuples() - average);
		standardDeviation /= (nOfRightCortical->GetNumberOfTuples()-1);
		standardDeviation = sqrt (standardDeviation);
		std::cout << "standard deviation: " << standardDeviation << "%" << std::endl;
		// wrong proportion
		// average
		average = 0.0;
		for (int a=0; a<corticalDistance->GetNumberOfTuples(); a++)
			average += corticalDistance->GetValue(a);
		average /= corticalDistance->GetNumberOfTuples();
		std::cout << "average distance of the wrong nodes from the cortical mask: " << average << std::endl;
		// standard deviation
		standardDeviation = 0.0;
		for (int a=0; a<corticalDistance->GetNumberOfTuples(); a++)
			standardDeviation += (corticalDistance->GetValue(a) - average)*(corticalDistance->GetValue(a) - average);
		standardDeviation /= (corticalDistance->GetNumberOfTuples()-1);
		standardDeviation = sqrt (standardDeviation);
		std::cout << "standard deviation: " << standardDeviation << std::endl;
		// standard error
		double standardError = 0.0;
		standardError = standardDeviation / sqrt(double(corticalDistance->GetNumberOfTuples()));
		std::cout << "standard error: " << standardError << std::endl;
		std::cout << std::endl;

		// measurements - trabecular
		// right proportion
		// average
		average = 0.0;
		for (int a=0; a<nOfRightTrabecular->GetNumberOfTuples(); a++)
			average += (nOfRightTrabecular->GetValue(a) * 100 / trabecular->GetNumberOfTuples());
		average /= nOfRightTrabecular->GetNumberOfTuples();
		std::cout << "right trabecular average: " << average << "%" << std::endl;
		// standard deviation
		 standardDeviation = 0.0;
		for (int a=0; a<nOfRightTrabecular->GetNumberOfTuples(); a++)
			standardDeviation += (nOfRightTrabecular->GetValue(a) * 100 / trabecular->GetNumberOfTuples() - average)*(nOfRightTrabecular->GetValue(a) * 100 / trabecular->GetNumberOfTuples() - average);
		standardDeviation /= (nOfRightTrabecular->GetNumberOfTuples()-1);
		standardDeviation = sqrt (standardDeviation);
		std::cout << "standard deviation: " << standardDeviation << "%" << std::endl;
		// wrong proportion
		// average
		average = 0.0;
		for (int a=0; a<trabecularDistance->GetNumberOfTuples(); a++)
			average += trabecularDistance->GetValue(a);
		average /= trabecularDistance->GetNumberOfTuples();
		std::cout << "average distance of the wrong nodes from the trabecular mask: " << average << std::endl;
		// standard deviation
		standardDeviation = 0.0;
		for (int a=0; a<trabecularDistance->GetNumberOfTuples(); a++)
			standardDeviation += (trabecularDistance->GetValue(a) - average)*(trabecularDistance->GetValue(a) - average);
		standardDeviation /= (trabecularDistance->GetNumberOfTuples()-1);
		standardDeviation = sqrt (standardDeviation);
		std::cout << "standard deviation: " << standardDeviation << std::endl;
		// standard error
		standardError = 0.0;
		standardError = standardDeviation / sqrt(double(trabecularDistance->GetNumberOfTuples()));
		std::cout << "standard error: " << standardError << std::endl;
		std::cout << std::endl;

		// measurements - marrow
		// right proportion
		// average
		average = 0.0;
		for (int a=0; a<nOfRightMarrow->GetNumberOfTuples(); a++)
			average += (nOfRightMarrow->GetValue(a) * 100 / marrow->GetNumberOfTuples());
		average /= nOfRightMarrow->GetNumberOfTuples();
		std::cout << "right marrow average: " << average << "%" << std::endl;
		// standard deviation
		standardDeviation = 0.0;
		for (int a=0; a<nOfRightMarrow->GetNumberOfTuples(); a++)
			standardDeviation += (nOfRightMarrow->GetValue(a) * 100 / marrow->GetNumberOfTuples() - average)*(nOfRightMarrow->GetValue(a) * 100 / marrow->GetNumberOfTuples() - average);
		standardDeviation /= (nOfRightMarrow->GetNumberOfTuples()-1);
		standardDeviation = sqrt (standardDeviation);
		std::cout << "standard deviation: " << standardDeviation << "%" << std::endl;
		// wrong proportion
		// average
		average = 0.0;
		for (int a=0; a<marrowDistance->GetNumberOfTuples(); a++)
			average += marrowDistance->GetValue(a);
		average /= marrowDistance->GetNumberOfTuples();
		std::cout << "average distance of the wrong nodes from the marrow mask: " << average << std::endl;
		// standard deviation
		standardDeviation = 0.0;
		for (int a=0; a<marrowDistance->GetNumberOfTuples(); a++)
			standardDeviation += (marrowDistance->GetValue(a) - average)*(marrowDistance->GetValue(a) - average);
		standardDeviation /= (marrowDistance->GetNumberOfTuples()-1);
		standardDeviation = sqrt (standardDeviation);
		std::cout << "standard deviation: " << standardDeviation << std::endl;
		// standard error
		standardError = 0.0;
		standardError = standardDeviation / sqrt(double(marrowDistance->GetNumberOfTuples()));
		std::cout << "standard error: " << standardError << std::endl;
		std::cout << std::endl;	

		// cleaning up
		cortical->Delete();
		trabecular->Delete();
		marrow->Delete();
		corticalDistance->Delete();
		trabecularDistance->Delete();
		marrowDistance->Delete();
		nOfRightCortical->Delete();
		nOfRightTrabecular->Delete();
		nOfRightMarrow->Delete();
				
				
	}

	


	void RegistrationValidarionVolumeMesh::ReferenceLabels(vtkPolyData* mesh, ImageHandler::ImageType::Pointer image, 
		vtkIntArray* &cortical, vtkIntArray* &trabecular, vtkIntArray* &marrow){
		
		vtkPoints* cort = vtkPoints::New();
		vtkPoints* trab = vtkPoints::New();
		vtkPoints* marr = vtkPoints::New();
			
		// iteration on the reference nodes
		for (int i=0; i<mesh->GetNumberOfPoints(); i++){
			
			// get node
			double node[3];
			mesh->GetPoint(i, node);
			
			// resolution (spacing)
			const ImageHandler::ImageType::SpacingType& resolution = image->GetSpacing();
			
			// image origin
			const ImageHandler::ImageType::PointType& origin = image->GetOrigin();
			//std::cout << origin[0] << ' ' << origin[1] << ' ' << origin[2] << std::endl;

			// grid position
			ImageHandler::ImageType::IndexType gridPosition;
			for (int dim = 0; dim < 3; dim++)
				gridPosition[dim] = std::floor((node[dim] - origin[dim]) /resolution[dim] + 0.5) ;
				// the node coordinate is in the physical space. the origin must be subtracted here, since
				// it is still in the physical space. the division by the resolution puts the image
				// in the image space, i.e. the grid

			// get voxel grey value
			ImageHandler::ImageType::ValueType value;
			value = image->GetPixel(gridPosition);

			// check what this voxel belongs to
			if (value == 300){
				cortical->InsertNextValue(i); cort->InsertNextPoint(node);
			}
			else if (value == 200){
				trabecular->InsertNextValue(i); trab->InsertNextPoint(node);
			}
			else if (value == 100){
				marrow->InsertNextValue(i); marr->InsertNextPoint(node);
			}
			else if (value == 0){
				cortical->InsertNextValue(i); cort->InsertNextPoint(node);
			}
		}
		
		std::cout << "number of cortical voxels: " << cortical->GetNumberOfTuples() << std::endl;
		std::cout << "number of trabecular voxels: " << trabecular->GetNumberOfTuples() << std::endl;
		std::cout << "number of marrow voxels: " << marrow->GetNumberOfTuples() << std::endl;
		std::cout << "total number of voxels by the arrays: " << cortical->GetNumberOfTuples() 
			+ trabecular->GetNumberOfTuples() + marrow->GetNumberOfTuples() << std::endl;

		//std::cout << "total number nodes in the mesh: " << mesh->GetNumberOfPoints() << std::endl;

		PointWriterXyz* writer = new PointWriterXyz;
		writer->SetFileName("cortical.txt");
		writer->SetInput(cort);
		//writer->Update();
		writer->SetFileName("trabecular.txt");
		writer->SetInput(trab);
		//writer->Update();
		writer->SetFileName("marrow.txt");
		writer->SetInput(marr);
		//writer->Update();
		
		// cleaning up
		cort->Delete();
		trab->Delete();
		marr->Delete();
		delete writer;

	}


	
	void RegistrationValidarionVolumeMesh::checkMovings(vtkIntArray* index, vtkPolyData* mesh, ImageHandler::ImageType::Pointer image, int part, double &distance, int &nOfRight){
		
		
		// variables that depend on the part of the bone that we are considering
		int greyLevel = 0;
		QString rightNodesFileName; QString wrongNodesFileName; QString closestNodesFileName;
		QString partMaskFileName; QString partMeshFileName; 
		if (part == 1){ // cortical
			std::cout << "evaluating corticals" << std::endl;
			greyLevel = 300;
			rightNodesFileName = "rightCortical.txt"; wrongNodesFileName = "wrongCortical.txt"; closestNodesFileName = "closestCortical.txt";
			partMaskFileName = "cortical.mhd"; partMeshFileName = "cortical.stl";
		}
		if (part == 2){ // trabecular
			std::cout << "evaluating trabecular" << std::endl;
			greyLevel = 200;
			rightNodesFileName = "rightTrabecular.txt";
			wrongNodesFileName = "wrongTrabecular.txt";
			partMaskFileName = "trabecular.mhd"; partMeshFileName = "trabecular.stl"; closestNodesFileName = "closestTrabecular.txt";
		}
		if (part == 3){ // marrow
			std::cout << "evaluating marrow" << std::endl;
			greyLevel = 100;
			rightNodesFileName = "rightMarrow.txt";
			wrongNodesFileName = "wrongMarrow.txt";
			partMaskFileName = "marrow.mhd"; partMeshFileName = "marrow.stl"; closestNodesFileName = "closestMarrow.txt";
		}
		
		// variables that do not dipend on the part of the bone that we are considering
		nOfRight = 0;
		distance = 0.0;
		const ImageHandler::ImageType::SpacingType& spacing = image->GetSpacing();
		const ImageHandler::ImageType::PointType& origin = image->GetOrigin();
		ImageHandler::ImageType::IndexType gridPosition;
		ImageHandler::ImageType::ValueType value;
		vtkPoints* rightNodes = vtkPoints::New();
		vtkPoints* wrongNodes = vtkPoints::New();
		vtkPoints* closestNodes = vtkPoints::New();
		PointWriterXyz* pointWriter = new PointWriterXyz;

		
		// check where the nodes are
		for (int y=0; y<index->GetNumberOfTuples(); y++){
			// node cortical coordinate
			double node[3];
			mesh->GetPoint(index->GetValue(y), node);
			// get voxel grey value
			for (int dim = 0; dim < 3; dim++)
				gridPosition[dim] = std::floor((node[dim] - origin[dim]) /spacing[dim] + 0.5) ; //unique situation in which the resolution is involved
			value = image->GetPixel(gridPosition);
			if (value == greyLevel){  // right nodes
				rightNodes->InsertNextPoint(node);
				nOfRight ++;
			}
			else // wrong nodes
				wrongNodes->InsertNextPoint(node); 
		}
		pointWriter->SetFileName(rightNodesFileName);
		pointWriter->SetInput(rightNodes);
		//pointWriter->Update();
		pointWriter->SetFileName(wrongNodesFileName);
		pointWriter->SetInput(wrongNodes);
		//pointWriter->Update();

		
		// creation of the mask for each part
		ImageHandler::ImageType::Pointer partMask = ImageHandler::ImageType::New();
		partMask->SetOrigin(origin);
		partMask->SetSpacing(spacing);
		ImageHandler::ImageType::RegionType region;
		region.SetSize( image->GetBufferedRegion().GetSize() );
		partMask->SetRegions(region);
		partMask->Allocate();
		partMask->FillBuffer(0);
		typedef itk::ImageRegionIterator< ImageHandler::ImageType> IteratorType;	
		IteratorType it (image, image->GetRequestedRegion());
		IteratorType it2 (partMask, partMask->GetRequestedRegion());
		for (it.GoToBegin(), it2.GoToBegin(); !it.IsAtEnd(); ++it, ++it2){
			if (it.Get() == greyLevel)
				it2.Set(100); // all the part masks put to 100 for the following threshold
		}
		ImageHandler* imageHandler = new ImageHandler;
		imageHandler->SetImage(partMask);
		imageHandler->SetImageFileName(partMaskFileName);
		//imageHandler->MetafileWriterUpdate();
		
		
		// creation of the mesh for each part
		imageHandler->SetThreshold(90);
		imageHandler->MarchingCubesUpdate();
		vtkSTLWriter* writer = vtkSTLWriter::New();
		writer->SetFileName(partMeshFileName);
		writer->SetInput(imageHandler->GetMesh());
		//writer->Update();

		// calculate distance for wrong nodes
		vtkCellLocator* cellLocator = vtkCellLocator::New();
		cellLocator->SetDataSet(imageHandler->GetMesh());
		cellLocator->SetNumberOfCellsPerBucket(1);
		cellLocator->BuildLocator();

		
		// find closest points
		for (int i=0; i<wrongNodes->GetNumberOfPoints(); i++){ // trabecular
			double point[3]; double closestPoint[3]; vtkIdType cellId; int subId; double dist;
			wrongNodes->GetPoint(i,point);
			cellLocator->FindClosestPoint(point, closestPoint, cellId, subId, dist);
			double temp = std::sqrt((point[0]-closestPoint[0]) * (point[0]-closestPoint[0]) +
								    (point[1]-closestPoint[1]) * (point[1]-closestPoint[1]) +
								    (point[2]-closestPoint[2]) * (point[2]-closestPoint[2]) );
			distance += temp;
			closestNodes->InsertNextPoint(closestPoint);

		}
		distance /= wrongNodes->GetNumberOfPoints(); // average distance
		std::cout << "average distance of the wrong nodes: " << distance << " mm" << std::endl; 

		pointWriter->SetFileName(closestNodesFileName);
		pointWriter->SetInput(closestNodes);
		//pointWriter->Update();

		// cleaning up
		rightNodes->Delete();
		wrongNodes->Delete();
		closestNodes->Delete();
		delete pointWriter;
		partMask = NULL;
		delete imageHandler;
		writer->Delete();
		cellLocator->Delete();
	}
	
	/*void RegistrationValidarionVolumeMesh::checkMovings(vtkIntArray* index, vtkPolyData* mesh, ImageHandler::ImageType::Pointer image, int part, double &distance, int &nOfRight){
	
		int corticalIndex = 0; 
		int trabecularIndex = 0;
		int marrowIndex = 0;
		
		const ImageHandler::ImageType::SpacingType& spacing = image->GetSpacing();
		const ImageHandler::ImageType::PointType& origin = image->GetOrigin();
		ImageHandler::ImageType::IndexType gridPosition;
		ImageHandler::ImageType::ValueType value;

		vtkPoints* cort = vtkPoints::New();
		vtkPoints* trab = vtkPoints::New();
		vtkPoints* marr = vtkPoints::New();
		PointWriterXyz* pointWriter = new PointWriterXyz;

		if (part == 3){//
			vtkPoints* ciao = vtkPoints::New();//
			for (int y=0; y<index->GetNumberOfTuples(); y++){//

				double pt[3]; pt[1]=0; pt[2]=0;//
				pt[0]=index->GetValue(y);//
				ciao->InsertNextPoint(pt);//
			}//
			
		pointWriter->SetFileName("index.txt");//
		pointWriter->SetInput(ciao);//
		pointWriter->Update();//

		pointWriter->SetFileName("meshpoints.txt");//
		pointWriter->SetInput(mesh->GetPoints());//
		pointWriter->Update();//
		}//


		// check if the voxels are there
		for (int y=0; y<index->GetNumberOfTuples(); y++){
			
			// node cortical coordinate
			double node[3];
			mesh->GetPoint(index->GetValue(y), node);
			// get voxel grey value
			for (int dim = 0; dim < 3; dim++)
				gridPosition[dim] = std::floor((node[dim] - origin[dim]) /spacing[dim] + 0.5) ; //unique situation in which the resolution is involved
			value = image->GetPixel(gridPosition);
			
			// label nodes
			if (value == 300) {
				corticalIndex ++; 
				cort->InsertNextPoint(node);
			}
			else if (value == 200) {
				trabecularIndex ++;
				trab->InsertNextPoint(node);
			}
			else if (value == 100) {
				marrowIndex ++;
				marr->InsertNextPoint(node);
			}
			else if (value == 0) {
				corticalIndex ++;
				cort->InsertNextPoint(node);
			}
		}			
		//std::cout << std::endl;
		//std::cout << "cortical: " << corticalIndex << std::endl;
		//std::cout << "trabecular: " << trabecularIndex << std::endl;
		//std::cout << "marrow: " << marrowIndex << std::endl;
		//std::cout << "total: " << (corticalIndex + trabecularIndex + marrowIndex) << std::endl; 
		
		if (part == 1){
			pointWriter->SetInput(cort);
			pointWriter->SetFileName("cort.txt");
			nOfRight = corticalIndex;
		}
		else if	(part == 2){
			pointWriter->SetInput(trab);
			pointWriter->SetFileName("trab.txt");
			nOfRight = trabecularIndex;
		}
		else if	(part == 3){
			pointWriter->SetInput(marr);
			pointWriter->SetFileName("marr.txt");
			nOfRight = marrowIndex;
		}
		//pointWriter->Update();
		
		// extract the mask of interest
		ImageHandler::ImageType::Pointer partMask = ImageHandler::ImageType::New();
		partMask->SetOrigin(origin);
		partMask->SetSpacing(spacing);
		ImageHandler::ImageType::RegionType region;
		region.SetSize( image->GetBufferedRegion().GetSize() );
		partMask->SetRegions(region);
		partMask->Allocate();
		partMask->FillBuffer(0);

		typedef itk::ImageRegionIterator< ImageHandler::ImageType> IteratorType;	
		IteratorType it (image, image->GetRequestedRegion());
		IteratorType it2 (partMask, partMask->GetRequestedRegion());
		
		for (it.GoToBegin(), it2.GoToBegin(); !it.IsAtEnd(); ++it, ++it2){
			if (part == 1){ // cortical
				if (it.Get() == 300)
					it2.Set(100);
			}
			if (part == 2){ // trabecular
				if (it.Get() == 200)
					it2.Set(100);
			}
			if (part == 3){ // marrow
				if (it.Get() == 100)
					it2.Set(100);
			}
		}
		

		
		// conversion to surface
		ImageHandler* imageHandler = new ImageHandler;
		imageHandler->SetImage(partMask);
		if (part == 1)
			imageHandler->SetImageFileName("cortical.mhd");
		if (part == 2)
			imageHandler->SetImageFileName("trabecular.mhd");
		if (part == 3)
			imageHandler->SetImageFileName("marrow.mhd");
		//imageHandler->MetafileWriterUpdate();
		
		imageHandler->SetThreshold(90);
		imageHandler->MarchingCubesUpdate();
		
		vtkSTLWriter* writer = vtkSTLWriter::New();
		if (part == 1)
			writer->SetFileName("cortical.stl");
		if (part == 2)
			writer->SetFileName("trabecular.stl");
		if (part == 3)
			writer->SetFileName("marrow.stl");
		writer->SetInput(imageHandler->GetMesh());
		//writer->Update();

		vtkCellLocator* cellLocator = vtkCellLocator::New();
		cellLocator->SetDataSet(imageHandler->GetMesh());
		cellLocator->SetNumberOfCellsPerBucket(1);
		cellLocator->BuildLocator();

		
		// find closest points
		if (part == 1){ // check trabecular and marrow points
			
			vtkPoints* wrongCortical = vtkPoints::New();
			
			distance = 0.0;
							
			for (int a=0; a<trab->GetNumberOfPoints(); a++){ // trabecular
					
				double point[3]; double closestPoint[3];
				vtkIdType cellId; int subId; double dist;
				trab->GetPoint(a,point);
				cellLocator->FindClosestPoint(point, closestPoint, cellId, subId, dist);

				double temp = std::sqrt((point[0]-closestPoint[0]) * (point[0]-closestPoint[0]) +
									    (point[1]-closestPoint[1]) * (point[1]-closestPoint[1]) +
									    (point[2]-closestPoint[2]) * (point[2]-closestPoint[2]) );
				distance += temp;
				wrongCortical->InsertNextPoint(point);
			}

			for (int a=0; a<marr->GetNumberOfPoints(); a++){ // marrow
					
				double point[3]; double closestPoint[3];
				vtkIdType cellId; int subId; double dist;
				marr->GetPoint(a,point);
				cellLocator->FindClosestPoint(point, closestPoint, cellId, subId, dist);

				double temp = std::sqrt((point[0]-closestPoint[0]) * (point[0]-closestPoint[0]) +
									    (point[1]-closestPoint[1]) * (point[1]-closestPoint[1]) +
									    (point[2]-closestPoint[2]) * (point[2]-closestPoint[2]) );
				distance += temp;
				wrongCortical->InsertNextPoint(point);
			}

			distance /= (trab->GetNumberOfPoints()+ marr->GetNumberOfPoints());
			std::cout << "distance wrong cortical: " << distance << std::endl;
			pointWriter->SetInput(wrongCortical);
			pointWriter->SetFileName("wrongCortical.txt");
			pointWriter->Update();

			wrongCortical->Delete();
			
		}

		if (part == 2){ // check cortical and marrow points

			distance = 0.0;
			vtkPoints* wrongTrabecular = vtkPoints::New();
			
			for (int a=0; a<cort->GetNumberOfPoints(); a++){ // cortical
					
				double point[3]; double closestPoint[3];
				vtkIdType cellId; int subId; double dist;
				cort->GetPoint(a,point);
				cellLocator->FindClosestPoint(point, closestPoint, cellId, subId, dist);

				double temp = std::sqrt((point[0]-closestPoint[0]) * (point[0]-closestPoint[0]) +
									    (point[1]-closestPoint[1]) * (point[1]-closestPoint[1]) +
									    (point[2]-closestPoint[2]) * (point[2]-closestPoint[2]) );
				distance += temp;
				wrongTrabecular->InsertNextPoint(point);
			}

			for (int a=0; a<marr->GetNumberOfPoints(); a++){ // marrow
					
				double point[3]; double closestPoint[3];
				vtkIdType cellId; int subId; double dist;
				marr->GetPoint(a,point);
				cellLocator->FindClosestPoint(point, closestPoint, cellId, subId, dist);

				double temp = std::sqrt((point[0]-closestPoint[0]) * (point[0]-closestPoint[0]) +
									    (point[1]-closestPoint[1]) * (point[1]-closestPoint[1]) +
									    (point[2]-closestPoint[2]) * (point[2]-closestPoint[2]) );
				distance += temp;
				wrongTrabecular->InsertNextPoint(point);
			}
			
			distance /= (cort->GetNumberOfPoints()+ marr->GetNumberOfPoints());
			std::cout << "distance wrong trabecular: " << distance << std::endl;
			pointWriter->SetInput(wrongTrabecular);
			pointWriter->SetFileName("wrongTrabecular.txt");
			pointWriter->Update();
			
			wrongTrabecular->Delete();

		}

		if (part == 3){ // check cortical and trabecular points
		
			distance = 0.0;
			vtkPoints* wrongMarrow = vtkPoints::New();
			
			for (int a=0; a<cort->GetNumberOfPoints(); a++){ // cortical
					
				double point[3]; double closestPoint[3];
				vtkIdType cellId; int subId; double dist;
				cort->GetPoint(a,point);
				cellLocator->FindClosestPoint(point, closestPoint, cellId, subId, dist);

				double temp = std::sqrt((point[0]-closestPoint[0]) * (point[0]-closestPoint[0]) +
									    (point[1]-closestPoint[1]) * (point[1]-closestPoint[1]) +
									    (point[2]-closestPoint[2]) * (point[2]-closestPoint[2]) );
				distance += temp;
				wrongMarrow->InsertNextPoint(point);
			}

			for (int a=0; a<trab->GetNumberOfPoints(); a++){ // trabecular
					
				double point[3]; double closestPoint[3];
				vtkIdType cellId; int subId; double dist;
				trab->GetPoint(a,point);
				cellLocator->FindClosestPoint(point, closestPoint, cellId, subId, dist);

				double temp = std::sqrt((point[0]-closestPoint[0]) * (point[0]-closestPoint[0]) +
									    (point[1]-closestPoint[1]) * (point[1]-closestPoint[1]) +
									    (point[2]-closestPoint[2]) * (point[2]-closestPoint[2]) );
				distance += temp;
				wrongMarrow->InsertNextPoint(point);
			}
			distance /= (cort->GetNumberOfPoints()+ trab->GetNumberOfPoints());
			std::cout << "distance wrong marrow: " << distance << std::endl;
			pointWriter->SetInput(wrongMarrow);
			pointWriter->SetFileName("wrongMarrow.txt");
			pointWriter->Update();
		
		// cleaning up
		wrongMarrow->Delete();
		}

		// cleaning up
		cort->Delete();
		trab->Delete();
		marr->Delete();
		delete pointWriter;
		writer->Delete();
		cellLocator->Delete();
		delete imageHandler;
		partMask = NULL;
					

	}
	*/
}