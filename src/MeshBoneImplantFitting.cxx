/*
 * ISTB - University of Bern, Mauricio Reyes, Serena Bonaretti
 */

#include <MeshBoneImplantFitting.h>
#include <RenderingOBB.h>

#include <vtkCellData.h>
#include <vtkCellLocator.h>
#include <vtkDataArray.h>
#include <vtkMath.h>
#include <vtkOBBTree.h>
#include <vtkPoints.h>
#include <vtkPolyDataNormals.h>
#include <vtkProcrustesAlignmentFilter.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

using namespace rendering;

namespace mesh{
	
	// constructor
	MeshBoneImplantFitting::MeshBoneImplantFitting(){

		_boneMesh = vtkPolyData::New();
		_implantMesh = vtkPolyData::New();
		_yResolution = 3.0;
		_zResolution = 3.0;
		_yStep = 1.0;
		_zStep = 1.0;

	}

	// destructor
	MeshBoneImplantFitting::~MeshBoneImplantFitting(){
	}

	

	// member functions
	void MeshBoneImplantFitting::Update(){

		_boneMesh->GetNumberOfPoints();
		_implantMesh->GetNumberOfPoints();
	}

	void MeshBoneImplantFitting::ImplantPositionUpdate(){

		
		// calculate implant OBB
		RenderingOBB* implantOBB = new RenderingOBB;
		implantOBB->SetMesh(_implantMesh);
		implantOBB->CalculateOBB();
		double oneI[3], twoI[3], threeI[3], fourI[3], fiveI[3], sixI[3], sevenI[3], eightI[3];
		implantOBB->GetOBBverteces(oneI, twoI, threeI, fourI, fiveI, sixI, sevenI, eightI);
		double implantCenterOfMass[3];
		implantOBB->GetMeshCenterOfMass(implantCenterOfMass);

		// calculate bone OBB
		RenderingOBB* boneOBB = new RenderingOBB;
		boneOBB->SetMesh(_boneMesh);
		boneOBB->CalculateOBB();
		double oneB[3], twoB[3], threeB[3], fourB[3], fiveB[3], sixB[3], sevenB[3], eightB[3];
		boneOBB->GetOBBverteces(oneB, twoB, threeB, fourB, fiveB, sixB, sevenB, eightB);
		double boneCenterOfMass[3];
		boneOBB->GetMeshCenterOfMass(boneCenterOfMass);

		// rigirare le OBB!!!


		// allignment
								
		//  
		TranslateMesh(_implantMesh, -implantCenterOfMass[0], -implantCenterOfMass[1], -implantCenterOfMass[2]);
		
		// rotation of the implant to the bone coordinate system
		RotateMesh(_implantMesh, oneI, twoI, threeI, fourI, oneB, twoB, fourB, threeB);

		// translation of the implant to the bone center of mass
		TranslateMesh(_implantMesh, boneCenterOfMass[0]+40, boneCenterOfMass[1], boneCenterOfMass[2]);



		
		/*
		vtkOBBTree* implantBoundingBox = vtkOBBTree::New();
		double Icorner[3]; double Imax[3]; double Imid[3]; double Imin[3]; double Isize[3];
		implantBoundingBox->ComputeOBB(_implantMesh->GetPoints(), Icorner, Imax, Imid, Imin, Isize); 
	
		vtkOBBTree* boneBoundingBox = vtkOBBTree::New();
		double Bcorner[3]; double Bmax[3]; double Bmid[3]; double Bmin[3]; double Bsize[3];
		boneBoundingBox->ComputeOBB(_implantMesh->GetPoints(), Bcorner, Bmax, Bmid, Bmin, Bsize); 
		
		//Compute center of mass of bone and implant
		double comBone[3]; double comImplant[3];
		CenterOfMass(_boneMesh,comBone);
		CenterOfMass(_implantMesh,comImplant);

		std::cout << comBone[0] << ' ' << comBone[1] << ' ' << comBone[2] << std::endl;
		std::cout << comImplant[0] << ' ' << comImplant[1] << ' ' << comImplant[2] << std::endl;
		
		double dx; double dy; double dz; 

		if (_positionFlag == 1){
			dx = comBone[0] + comImplant[0];
			dy = comBone[1] + comImplant[1];
			dz = comBone[2] + comImplant[2];
		}
		
		if (_positionFlag == 2){
			dx = comBone[0] - comImplant[0];
			dy = comBone[1] - comImplant[1];
			dz = comBone[2] - comImplant[2];
		}

		// allignment
		TranslateVolume(dx,dy,dz,_implantMesh);
		
		CenterOfMass(_implantMesh,comImplant);
		std::cout << comImplant[0] << ' ' << comImplant[1] << ' ' << comImplant[2] << std::endl;
	/*
		// collision
		vtkPoints* distances = vtkPoints::New();
		CheckCollision (_implantMesh, _boneMesh, distances);
		

		
		for (int i=0; distances->GetNumberOfPoints(); i++){
			double temp[3];
			distances->GetPoint(i, temp);
			std::cout << temp[0] << ' ' << temp[1] << ' ' << temp[2] << std::endl;
		}
		*/
		
		/*
		double xMax = 0.0; double yMax = 0.0; double zMax = 0.0;
		for (int i=0; i<distances->GetNumberOfPoints(); i++){
			double temp[3];
			distances->GetPoint(i, temp);
			if (std::abs(temp[0]) > std::abs(xMax))
				xMax = temp[0];
			if (std::abs(temp[1]) > std::abs(yMax))
				yMax = temp[1];
			if (std::abs(temp[2]) > std::abs(zMax))
				zMax = temp[2];
		}

		std::cout << xMax << ' ' << yMax << ' ' << zMax << std::endl;
		TranslateVolume(xMax,yMax,zMax,_implantMesh);
		*/
		


	}


	// private
	void MeshBoneImplantFitting::TranslateMesh(vtkPolyData* mesh, double dx,double dy,double dz){

		vtkTransform *trInitial =vtkTransform::New();
		trInitial->Translate(dx,dy,dz);

		vtkTransformPolyDataFilter *pdf=vtkTransformPolyDataFilter::New();
		pdf->SetTransform(trInitial);
		pdf->SetInput(mesh);
		pdf->Update();

		mesh->DeepCopy(pdf->GetOutput());

		trInitial->Delete();
		pdf->Delete();
	}

	void MeshBoneImplantFitting::RotateMesh(vtkPolyData* mesh, double oneI[], double twoI[], double threeI[], double fourI[],
															   double oneB[], double twoB[], double threeB[], double fourB[]){
		
		// rotation to the main coordinate system (0,0,0)

		// rotation axis
		double maxVersorT[3], midVersorT[3], minVersorT[3];
		
		for(unsigned int i = 0; i < 3; ++i) {
			maxVersorT[i]=(twoI[i]-oneI[i]);
			midVersorT[i]=(threeI[i]-oneI[i]);
			minVersorT[i]=(fourI[i]-oneI[i]);
		}
		
		// rotation axis norms
		double normMax = vtkMath::Norm(maxVersorT);
		double normMid = vtkMath::Norm(midVersorT);
		double normMin = vtkMath::Norm(minVersorT);
		
		// rotation versors
		for(unsigned int i = 0; i < 3; ++i) {
			maxVersorT[i]=maxVersorT[i]/normMax;
			midVersorT[i]=midVersorT[i]/normMid;
			minVersorT[i]=minVersorT[i]/normMin;
		}	
		
		// rotation matrix
  		vtkMatrix4x4* to = vtkMatrix4x4::New();
		to->Element[0][0]=midVersorT[0]; to->Element[0][1]=midVersorT[1]; to->Element[0][2]=midVersorT[2]; to->Element[0][3]=0.0; 
		to->Element[1][0]=minVersorT[0]; to->Element[1][1]=minVersorT[1]; to->Element[1][2]=minVersorT[2]; to->Element[1][3]=0.0; 
		to->Element[2][0]=maxVersorT[0]; to->Element[2][1]=maxVersorT[1]; to->Element[2][2]=maxVersorT[2]; to->Element[2][3]=0.0; 
		to->Element[3][0]=0.0;			 to->Element[3][1]=0.0;           to->Element[3][2]=0.0;           to->Element[3][3]=1.0;

		// tranformations 
		vtkTransform* rotateImage = vtkTransform::New();
		rotateImage->SetMatrix(to);
		
		// execute rotation
		vtkTransformPolyDataFilter* rotationTo = vtkTransformPolyDataFilter::New();
		rotationTo->SetInput(mesh);
		rotationTo->SetTransform(rotateImage);
		rotationTo->Update();
		mesh->DeepCopy(rotationTo->GetOutput());

		
		
		// rotation to the bone coordinate system

		// rotation to the new position
		for(unsigned int i = 0; i < 3; ++i) {
			maxVersorT[i]=(twoB[i]-oneB[i]);
			midVersorT[i]=(threeB[i]-oneB[i]);
			minVersorT[i]=(fourB[i]-oneB[i]);
		}
		
		// rotation axis norms
		normMax = vtkMath::Norm(maxVersorT);
		normMid = vtkMath::Norm(midVersorT);
		normMin = vtkMath::Norm(minVersorT);
		
		// rotation versors
		for(unsigned int i = 0; i < 3; ++i) {
			maxVersorT[i]=maxVersorT[i]/normMax;
			midVersorT[i]=midVersorT[i]/normMid;
			minVersorT[i]=minVersorT[i]/normMin;
		}	

		vtkMatrix4x4* back = vtkMatrix4x4::New();
		back->Element[0][0]=midVersorT[0]; back->Element[0][1]=minVersorT[0]; back->Element[0][2]=maxVersorT[0]; back->Element[0][3]=0.0; 
		back->Element[1][0]=midVersorT[1]; back->Element[1][1]=minVersorT[1]; back->Element[1][2]=maxVersorT[1]; back->Element[1][3]=0.0; 
		back->Element[2][0]=midVersorT[2]; back->Element[2][1]=minVersorT[2]; back->Element[2][2]=maxVersorT[2]; back->Element[2][3]=0.0; 
		back->Element[3][0]=0.0;		   back->Element[3][1]=0.0;           back->Element[3][2]=0.0;           back->Element[3][3]=1.0;
		
		// tranformations 
		vtkTransform* rotateImageBack = vtkTransform::New();
		rotateImageBack->SetMatrix(back);
		
		// execute rotation
		vtkTransformPolyDataFilter* rotationBack = vtkTransformPolyDataFilter::New();
		rotationBack->SetInput(mesh);
		rotationBack->SetTransform(rotateImageBack);
		rotationBack->Update();
		mesh->DeepCopy(rotationBack->GetOutput());
		
	
	
	}
	/*
	void MeshBoneImplantFitting::CenterOfMass(vtkPolyData *data,double com[3]){
		
		int nbpts=data->GetNumberOfPoints();
		double pt[3];
		
		double xcom=0.0;
		double ycom=0.0;
		double zcom=0.0;

		for(int i = 0 ; i < nbpts ; i++)
		{
			data->GetPoint(i,pt);	
			xcom+=pt[0];
			ycom+=pt[1];
			zcom+=pt[2];
		}

		com[0]=xcom/nbpts;
		com[1]=ycom/nbpts;
		com[2]=zcom/nbpts;
	}
	void MeshBoneImplantFitting::CheckCollision(vtkPolyData* implantMesh, vtkPolyData* boneMesh, vtkPoints* distances) {
		
		//int collision=0;
		
		
		//make temp copy of implant
		//vtkPolyData *myimplant=vtkPolyData::New();
		//myimplant->DeepCopy(implantPoly);

		//for each point in the implant check if it goes inside the bone mesh
		
		
		for (int i=0; i<implantMesh->GetNumberOfPoints(); i=i+4){
			
			double pt[3];
			implantMesh->GetPoints()->GetPoint(i,pt);
			//std::cout << "implant point: " << pt[0] << ' ' << pt[1] << ' ' << pt[2] << std::endl;
			
			double x; double y; double z;
			InVolume(boneMesh, pt, x, y, z);
			//std::cout << "distances: " << x << ' ' << y << ' ' << z << std::endl;

			distances->InsertNextPoint(x, y, z);
			//if ( InVolume (boneMesh,pt, x, y,z) ){
			//	collision=1;
		//		break;
			//}
		}
		//myimplant->Delete();
		//return collision;
	}

	void MeshBoneImplantFitting::InVolume(vtkPolyData* mesh, double point[3], double& x, double& y, double& z){
		
		// closest point
		vtkCellLocator* locator = vtkCellLocator::New();
		locator->SetDataSet(mesh);
		locator->SetNumberOfCellsPerBucket(1);
		locator->BuildLocator();
	  	vtkIdType cell_id; int sub_id; double dist2;
		//double p[3]; p[0]=point[0]; p[1]=point[1]; p[2]=point[2];
		double pt[3]; double normal[3];
		locator->FindClosestPoint(point,pt,cell_id,sub_id,dist2);
	  
		//std::cout << "closest point: " << pt[0] << ' ' << pt[1] << ' ' << pt[2] << std::endl;
	    
		// normals 
		vtkPolyDataNormals* pdnormals=vtkPolyDataNormals::New();
		pdnormals->SetInput(mesh);
		pdnormals->SetFeatureAngle(90.0);
		pdnormals->ConsistencyOn();
		pdnormals->ComputeCellNormalsOn();
		pdnormals->FlipNormalsOn();
		pdnormals->Update();
		
		vtkDataArray* normales= pdnormals->GetOutput()->GetCellData()->GetNormals();
		normal[0] = normales->GetComponent(cell_id,0);
	    normal[1] = normales->GetComponent(cell_id,1);
	    normal[2] = normales->GetComponent(cell_id,2);
	  
		//vector from closest point to input point
		double cp[3];
		cp[0]=point[0]-pt[0];
		cp[1]=point[1]-pt[1];
		cp[2]=point[2]-pt[2];
				
		vtkMath* math=vtkMath::New();
		//int ret;
		if(math->Dot(cp,normal) >= 0){
			//ret=1;
			x = cp[0]; y = cp[1]; z = cp[2];
		}
		else {
			//ret=0;
			x = 0.0; y = 0.0; z = 0.0;
		}

	    math->Delete();
	    locator->Delete();
	    pdnormals->Delete();
		
		//return ret;
	
	}
	*/

} 
