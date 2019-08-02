/***********************************************************************
/
/  GRID CLASS (APPLY A SN EXPLOSION TIME-ACTION TO THIS GRID)
/
/  written by: Miao Li
/  date:       2018
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "CosmologyParameters.h"
#include "phys_constants.h"

#define MAX_TEMPERATURE 1e8

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, double *MassUnits, FLOAT Time);

int FindField(int field, int farray[], int numfields);



int grid::ApplySNExplosion_Miao(FLOAT pos[3], float SN_radius)
{
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  int dim;
  float DomainWidth[3], GridLeft[3], GridRight[3];
  for (dim = 0; dim < GridRank; dim++){
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
    GridLeft[dim] = CellLeftEdge[dim][0];  //including ghost zones
    GridRight[dim] = CellLeftEdge[dim][GridDimension[dim]-1] + CellWidth[dim][0]; //including ghost zones
  }

  /*mirroring 8 SN*/
  /*check if 8 sphere overlaps with this grid*/

  for (dim = 0; dim < GridRank-1; dim++) {
    if ( (pos[dim] - SN_radius > GridRight[dim] ||
        pos[dim] + SN_radius < GridLeft[dim] ) &&
        (pos[dim] + DomainWidth[dim] - SN_radius > GridRight[dim] ||
        pos[dim] + DomainWidth[dim] + SN_radius < GridLeft[dim]) &&
       ( pos[dim] - DomainWidth[dim] - SN_radius > GridRight[dim] ||
        pos[dim] - DomainWidth[dim] + SN_radius < GridLeft[dim]) ) 
      return SUCCESS;
  }

  /*check sphere in z-direction, if overlapping with this grid*/

  dim = 2;
  if (pos[dim] - SN_radius > GridRight[dim] ||
      pos[dim] + SN_radius < GridLeft[dim] )
    return SUCCESS;

  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num,CRNum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum,CRNum) == FAIL) {
        ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  /* Find Metallicity or SNColour field and set flag. */

  int SNColourNum, MetalNum, Metal2Num, MBHColourNum, Galaxy1ColourNum,
    Galaxy2ColourNum, MetalIaNum, MetalIINum;
  int MetallicityField = FALSE;

  if (this->IdentifyColourFields(SNColourNum, Metal2Num, MetalIaNum,
                                 MetalIINum, MBHColourNum, Galaxy1ColourNum,
                                 Galaxy2ColourNum) == FAIL)
    ENZO_FAIL("Error in grid->IdentifyColourFields.\n");


/***********************************************************************
 *                                 SUPERNOVAE
 ***********************************************************************/

  float DensityUnits, LengthUnits,TemperatureUnits, TimeUnits,VelocityUnits;
  double MassUnits;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, 
	       InitialTimeInCodeUnits) == FAIL) {
    ENZO_FAIL("Error in GetUnits in CheckForTimeAction.C.");
  }

  int i,j,k,index;

  float maxGE, delz,sz,dely,sy,delx,sx,radius2;
  maxGE = MAX_TEMPERATURE / (TemperatureUnits * (Gamma-1.0) * 0.6);

  double SN_mass, SN_energy, den_add, energy_den_add, old_den, SN_colour_add;
  SN_mass = 10; //solar masses
  SN_energy = 1e51; //erg
  den_add = SN_mass*SolarMass/MassUnits/(4*pi/3.*POW(SN_radius,3));
  energy_den_add = SN_energy/(MassUnits* POW(VelocityUnits,2))/(4*pi/3.*POW(SN_radius,3));
  SN_colour_add = 1./(4*pi/3.*POW(SN_radius,3));

  printf("den_add=%e\n",den_add);
  printf("energy_den_add=%e\n",energy_den_add);
  printf("SN_colour_add=%e\n", SN_colour_add);
  
  for (k = 0; k < GridDimension[2]; k++) {

    delz = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - pos[2];
    sz = sign(delz);
    delz = fabs(delz);

//if no periodic boundary condition, comment the following out
//      delz = min(delz, DomainWidth[2]-delz);

    for (j = 0; j < GridDimension[1]; j++) {

      dely = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - pos[1];
      sy = sign(dely);
      dely = fabs(dely);
      dely = min(dely, DomainWidth[1]-dely);

      index = (k*GridDimension[1] + j)*GridDimension[0];
      for (i = 0; i < GridDimension[0]; i++, index++) {

	delx = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - pos[0];
	sx = sign(delx);
	delx = fabs(delx);
	delx = min(delx, DomainWidth[0]-delx);

	radius2 = delx*delx + dely*dely + delz*delz;
	if (radius2 <= pow(SN_radius,2)) {

	  // add SN: update density, energy, etc. 

	  printf("Inside SN sphere,i,j,k=%i, %i, %i.\n",i,j,k);
	  old_den = BaryonField[DensNum][index];
	  BaryonField[DensNum][index]= BaryonField[DensNum][index]+den_add  ;
	  float den_ratio;
	  den_ratio = old_den/BaryonField[DensNum][index];

	  if (CRModel){
	    BaryonField[CRNum][index] = SNEnergyCRFraction * energy_den_add;
	  }

	  if (SNThermalFeedback) {
	    if (GENum>0){
	      BaryonField[GENum][index]= ( old_den* BaryonField[GENum][index]+(1.0-SNEnergyCRFraction*CRModel)*energy_den_add) / BaryonField[DensNum][index];

	      BaryonField[Vel1Num][index] =  BaryonField[Vel1Num][index]*pow(den_ratio,0.5);
	      BaryonField[Vel2Num][index] =  BaryonField[Vel2Num][index]*pow(den_ratio,0.5);
	      BaryonField[Vel3Num][index] =  BaryonField[Vel3Num][index]*pow(den_ratio,0.5);

	      BaryonField[TENum][index] = BaryonField[GENum][index]+0.5*( pow(BaryonField[Vel1Num][index],2)+pow(BaryonField[Vel2Num][index],2)+pow(BaryonField[Vel3Num][index],2));

	      float e_ratio = 0.5*( pow(BaryonField[Vel1Num][index],2)+pow(BaryonField[Vel2Num][index],2)+pow(BaryonField[Vel3Num][index],2))/ BaryonField[GENum][index];
	      if (e_ratio>1.)
		printf("e_ratio,i,j,k,old_den=%e,%i, %i, %i,%e.\n",e_ratio,i,j,k,old_den);
	    }
	    else {
	      BaryonField[TENum][index]= ( old_den* BaryonField[TENum][index]+ (1.0- SNEnergyCRFraction*CRModel)*energy_den_add) / BaryonField[DensNum][index];
	    }
	  } // end if SNThermalFeedback

	  if (UseSNColour){
	    BaryonField[SNColourNum][index]= BaryonField[SNColourNum][index] + SN_colour_add;
	  } // end if UseSNColour

	} // end if (radius2 <= pow(SN_radius,2))
	
      }//end i
    }// end j
  } //end k

   return SUCCESS;
}

