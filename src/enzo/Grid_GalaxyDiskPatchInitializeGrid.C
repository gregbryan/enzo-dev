/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A GALAXY DISK PATCH SIMULATION)
/
/  written by: Miao Li
/  date:       2018
/  modified1:  
/  modified1:  
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "phys_constants.h"

int FindField(int field, int farray[], int numfields);
int GetUnits(float *DensityUnits, float *LengthUnits,
              float *TemperatureUnits, float *TimeUnits,
              float *VelocityUnits, double *MassUnits, FLOAT Time);


int grid::GalaxyDiskPatchInitializeGrid(int GasDistributionType,
                                        float GasDiskScaleHeight_pc,
                                        float MidplaneDensity,
                                        float MidplaneTotalEnergy,
                                        float MidplaneInternalEnergy,
                                        float MidplaneVelocity[],
                                        float UniformBField[])
{
  /* declarations */

  int dim, i, j, k, size,field,MetalNum;

  /* create fields */

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  int vel = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1 || HydroMethod > 2)
    FieldType[NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2 || HydroMethod > 2)
    FieldType[NumberOfBaryonFields++] = Velocity3;

  int  CRNum = NumberOfBaryonFields;
  if( CRModel )
    FieldType[NumberOfBaryonFields++] = CRDensity;

  if( UseMHD ){
    FieldType[NumberOfBaryonFields++] = Bfield1;
    FieldType[NumberOfBaryonFields++] = Bfield2;
    FieldType[NumberOfBaryonFields++] = Bfield3;
  }

  if( HydroMethod == MHD_RK ){
    FieldType[NumberOfBaryonFields++] = PhiField;
  }


  if (TestProblemData.UseMetallicityField) 
    FieldType[MetalNum = NumberOfBaryonFields++] = Metallicity;

  int ColourNum;
  if (UseSNColour){
    FieldType[ColourNum = NumberOfBaryonFields++] = SNColour;
   }

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* compute size of fields */

  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* allocate fields */

  this->AllocateGrids();

  int DenNum,TotalEnergyNum,InternalEnergyNum,ColourNum1;
  TotalEnergyNum = FindField(TotalEnergy, FieldType, NumberOfBaryonFields);
  InternalEnergyNum=FindField(InternalEnergy, FieldType, NumberOfBaryonFields);
  DenNum = FindField(Density,FieldType,NumberOfBaryonFields);
  ColourNum1 = FindField(SNColour,FieldType,NumberOfBaryonFields);

  /* Get Units */

  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1, TimeUnits = 1,
               VelocityUnits = 1, AccelUnits = 1;
  double MassUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
                    ENZO_FAIL("Error in GetUnits.");
  }
  AccelUnits = LengthUnits/TimeUnits/TimeUnits;

  /* Add a uniform or stratified medium. see Grid_GalaxySimulationInitializeGrid.C to add position-dependent BaryonFields at 33%*/

  int index,jndex,kndex; 
  FLOAT xpos,ypos,zpos;  
  FLOAT DomainWidth[3];
  for (dim = 0; dim < GridRank; dim++){
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
  }

  for (i=0;i<size;i++){
    index  = i % GridDimension[0];
    jndex  = (i-index) % (GridDimension[0]*GridDimension[1]);
    kndex  = (i-index - jndex)/(GridDimension[0]*GridDimension[1]);
    jndex /= GridDimension[0];

    xpos  = *(CellLeftEdge[0] + index) + 0.5*(*(CellWidth[0] + index));
    ypos  = *(CellLeftEdge[1] + jndex) + 0.5*(*(CellWidth[1] + jndex));
    zpos  = *(CellLeftEdge[2] + kndex) + 0.5*(*(CellWidth[2] + kndex));

    // GasDistributionType 0: uniform gas distribution 

    if (GasDistributionType == 0) { 

      BaryonField[0][i] = MidplaneDensity;
      BaryonField[1][i] = MidplaneTotalEnergy;

      /* set velocities */

      for (dim = 0; dim < GridRank; dim++)
	BaryonField[vel+dim][i] = MidplaneVelocity[dim];

      /* Set internal energy if necessary. */

      if (DualEnergyFormalism)
	BaryonField[2][i] = MidplaneInternalEnergy;

      if (TestProblemData.UseMetallicityField)
	BaryonField[MetalNum][i] = TestProblemData.MetallicityField_Fraction* BaryonField[0][i];

      if (CRModel)
	BaryonField[CRNum][i] = CREnergyDensityFloor;

      if (UseSNColour)
	BaryonField[ColourNum][i] = 1e-15;

    }//end GasDistributionType ==0

    // GasDistributionType 10:isothermal gas, density linearly decreases with z-position

    else if (GasDistributionType == 10 ){
      BaryonField[0][i] = MidplaneDensity - 0.99*MidplaneDensity * zpos/(DomainWidth[2]/2.);
      BaryonField[1][i] = MidplaneTotalEnergy*BaryonField[0][i]/MidplaneDensity;

      for (dim = 0; dim < GridRank; dim++)
	BaryonField[vel+dim][i] = MidplaneVelocity[dim];

      if (DualEnergyFormalism)
	BaryonField[2][i] = BaryonField[1][i];

      if (TestProblemData.UseMetallicityField)
	BaryonField[MetalNum][i] = TestProblemData.MetallicityField_Fraction* BaryonField[0][i];

      if ( CRModel )
	BaryonField[CRNum][i] = CREnergyDensityFloor;

      if (UseSNColour)
	BaryonField[ColourNum][i] = 1e-15;

    } // end GasDistributionType== 10

//  GasDistributionType 20: isothermal gas in hydrostatic equilibrium with external gravitational field 2

    else if (GasDistributionType == 20 ){

      float z_star,vel_dispersion_star_square,alpha,Sigma_star;
      z_star =  ExternalGravityScaleHeight_pc*pc_cm; //cm
      Sigma_star =  MassDensityForExternalGravity * SolarMass /POW(pc_cm,2); // g/cm^2
      vel_dispersion_star_square = z_star *3.141593 *GravConst * Sigma_star /POW(VelocityUnits,2);//code unit
      alpha = vel_dispersion_star_square /(StellarMassFraction*(Gamma -1.)*MidplaneTotalEnergy );
      BaryonField[0][i] =max(3e-28/DensityUnits, MidplaneDensity * POW(cosh(zpos/(z_star/LengthUnits) )  , -2*alpha) );
      BaryonField[1][i] = MidplaneTotalEnergy;
  
      for (dim = 0; dim < GridRank; dim++)
	BaryonField[vel+dim][i] = MidplaneVelocity[dim];

      if (DualEnergyFormalism)
	BaryonField[2][i] = BaryonField[1][i];

      if (TestProblemData.UseMetallicityField)
	BaryonField[MetalNum][i] = TestProblemData.MetallicityField_Fraction* BaryonField[0][i];

      if ( CRModel )
	BaryonField[CRNum][i] = CREnergyDensityFloor;

      if (UseSNColour)
	BaryonField[ColourNum][i] = 1e-15;

    }    // end GasDistributionType== 20  

//   GasDistributionType 30: SLICC Walch et al 15. MN454,238

    else if (GasDistributionType == 30){

      float hz = GasDiskScaleHeight_pc*pc_cm/LengthUnits;
           
      BaryonField[0][i] = max(MidplaneDensity * exp(-POW(zpos/hz,2)), 1e-27/DensityUnits);
      BaryonField[1][i] = MidplaneTotalEnergy;

      for (dim = 0; dim < GridRank; dim++)
	BaryonField[vel+dim][i] = MidplaneVelocity[dim];

      if (DualEnergyFormalism)
	BaryonField[2][i] = BaryonField[1][i];

      if (TestProblemData.UseMetallicityField)
	BaryonField[MetalNum][i] = TestProblemData.MetallicityField_Fraction* BaryonField[0][i];

      if ( CRModel )
	BaryonField[CRNum][i] = CREnergyDensityFloor;

      if (UseSNColour)
	BaryonField[ColourNum][i] = 1e-15;

    } // end GasDistributionType == 30

  } //end i

  if (UseMHD){

    for ( int Bfield=0; Bfield < 3; Bfield++ ){
      for ( int k=0; k<MagneticSize[Bfield]; k++ ){
        MagneticField[Bfield][k] = UniformBField[Bfield];
      }
    }  // end for field

    this->CenterMagneticField();

    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, B1Num, B2Num, B3Num;
    this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num,
				     TENum, B1Num, B2Num, B3Num);

    for (i=0;i<size;i++){
      BaryonField[ TENum][i] += 0.5*(BaryonField[B1Num][i]*BaryonField[B1Num][i]+
				     BaryonField[B2Num][i]*BaryonField[B2Num][i]+
				     BaryonField[B3Num][i]*BaryonField[B3Num][i])/BaryonField[ DensNum][i];
    } // end i 

  }//end if use MHD

  float zd = 100.*pc_cm/LengthUnits;
  float g0 = 2.63e-9/AccelUnits;
  printf ("scale factor_Miao: %e\n", cosh(1./2./zd ));
  printf ("scale factor_Miao: %e\n",-2.*zd*g0/MidplaneTotalEnergy );
  printf ("scale factor_Miao: %e\n",-2.*zd*g0/(Gamma-1.0)/MidplaneTotalEnergy );
  printf ("scale factor_Miao: %e\n",POW( cosh(1./2./zd), -2.*zd*g0/(Gamma-1.0)/MidplaneTotalEnergy));

  /*end of setting BaryonField*/

  TotalEnergyNum = FindField(TotalEnergy, FieldType, NumberOfBaryonFields);
  InternalEnergyNum=FindField(InternalEnergy, FieldType, NumberOfBaryonFields);

  return SUCCESS;

}

