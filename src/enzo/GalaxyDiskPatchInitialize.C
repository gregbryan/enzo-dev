/***********************************************************************
/
/  INITIALIZE a patch of galaxy disk
/
/  written by: Miao Li
/  adapted by: Nicole Melso
/
/  date:       Nov, 2017
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.

#include <string.h>
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
#include "Hierarchy.h"
#include "TopGridData.h"
#include "phys_constants.h"

void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
                     LevelHierarchyEntry *LevelArray[], int level);


int GetUnits(float *DensityUnits, float *LengthUnits,
       float *TemperatureUnits, float *TimeUnits,
       float *VelocityUnits, FLOAT Time);


int GalaxyDiskPatchInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			 TopGridData &MetaData, ExternalBoundary &Exterior)
{
const  char *DensName = "Density";
const  char *TEName   = "TotalEnergy";
const  char *GEName   = "GasEnergy";
const  char *Vel1Name = "x-velocity";
const  char *Vel2Name = "y-velocity";
const  char *Vel3Name = "z-velocity";
const  char *CRName   = "CREnergyDensity";

char *BxName = "Bx";
char *ByName = "By";
char *BzName = "Bz";
char *PhiName = "Phi";
const char *ColourName = "SN_Colour";

  /* parameter declarations */

  
  /* local declarations */

  char line[MAX_LINE_LENGTH];
  int  dim, ret, NumberOfSubgridZones[MAX_DIMENSION],
                          SubgridDims[MAX_DIMENSION];
 

  /* make sure this is 2D or 3D */

  if (MetaData.TopGridRank < 3 || MetaData.TopGridRank > 3) {
    ENZO_VFAIL("Cannot model GalaxyDiskPatch in %"ISYM" dimension(s)\n", MetaData.TopGridRank)
  }    

  /* parameters */

  int  GalaxyDiskPatchRefineAtStart         = 1;
  int  GalaxyDiskPatchGasDistributionType   = 0;      
  float GalaxyDiskPatchGasScaleHeight_pc    = 60.;           
  float GalaxyDiskPatchMidplaneDensity            = 1.0;
  float GalaxyDiskPatchTotalEnergy                = 1e-5;
  float GalaxyDiskPatchInternalEnergy             = 1e-5;
  float GalaxyDiskPatchVelocity[3]                = {0.0,0.0,0.0};
  float GalaxyDiskPatchBField[3]                  = {1e-16,0.0,0.0};
  
  /* nicole: boundary parameters (top & bottom) */ 
  int UseInflow = TRUE; 
  float BoundaryDensity[2] = {1.0, 1.0}; 
  float BoundaryTemperature[2] = {10000.0, 10000.0}; 
  float BoundaryVelocity[2] = {0.0, 0.0}; 
  /* end nicole */

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line,"GalaxyDiskPatchRefineAtStart  = %"ISYM, &GalaxyDiskPatchRefineAtStart);
    ret += sscanf(line,"GalaxyDiskPatchGasDistributionType  = %"ISYM, &GalaxyDiskPatchGasDistributionType);
    ret += sscanf(line,"GalaxyDiskPatchGasScaleHeight_pc  = %"FSYM, &GalaxyDiskPatchGasScaleHeight_pc);
    ret += sscanf(line,"GalaxyDiskPatchMidplaneDensity  = %"FSYM,&GalaxyDiskPatchMidplaneDensity);
    ret += sscanf(line,"GalaxyDiskPatchTotalEnergy  = %"FSYM,&GalaxyDiskPatchTotalEnergy);    
    ret += sscanf(line,"GalaxyDiskPatchInternalEnergy  = %"FSYM,&GalaxyDiskPatchInternalEnergy);    

    /* nicole: read boundary parameters */ 
    ret += sscanf(line,"GalaxyDiskPatchUseInflow  = %"ISYM, &UseInflow);
    
    if (UseInflow == TRUE) { 
    	ret += sscanf(line,"GalaxyDiskPatchTopDensity  = %"ESYM, &BoundaryDensity[0]);
    	ret += sscanf(line,"GalaxyDiskPatchBottomDensity  = %"ESYM, &BoundaryDensity[1]);
    	ret += sscanf(line,"GalaxyDiskPatchTopTemperature  = %"ESYM, &BoundaryTemperature[0]);
    	ret += sscanf(line,"GalaxyDiskPatchBottomTemperature  = %"ESYM, &BoundaryTemperature[1]);
    	ret += sscanf(line,"GalaxyDiskPatchTopVelocity  = %"ESYM, &BoundaryVelocity[0]);
    	ret += sscanf(line,"GalaxyDiskPatchBottomVelocity  = %"ESYM, &BoundaryVelocity[1]);
    }
    /* end nicole */
    
    /*read in metal field parameters by Miao  */
    ret += sscanf(line, "TestProblemUseMetallicityField  = %"ISYM, &TestProblemData.UseMetallicityField);
    ret += sscanf(line, "TestProblemInitialMetallicityFraction  = %"FSYM, &TestProblemData.MetallicityField_Fraction);
    /*end read in metal field*/

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "GalaxyDiskPatch") &&
        line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr,
         "warning: the following parameter line was not interpreted:\n%s\n",
              line);

  } // end input from parameter file

  printf("GDP_Density, GDP_TotalEnergy=%e,%e\n",GalaxyDiskPatchMidplaneDensity,GalaxyDiskPatchTotalEnergy);

  if (TopGrid.GridData->GalaxyDiskPatchInitializeGrid(GalaxyDiskPatchGasDistributionType,
                                               GalaxyDiskPatchGasScaleHeight_pc,
                                               GalaxyDiskPatchMidplaneDensity,
                                               GalaxyDiskPatchTotalEnergy,
                                              GalaxyDiskPatchInternalEnergy,
                                              GalaxyDiskPatchVelocity,
                                              GalaxyDiskPatchBField) == FAIL) {

    ENZO_FAIL("Error in GalaxyDiskPatchInitialize[Sub]Grid.");
  }  // end TopGrid if


  /* If requested, refine the grid to the desired level. */

  if (GalaxyDiskPatchRefineAtStart) {

    /* Declare, initialize and fill out the LevelArray. */
    int level;
    LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      LevelArray[level] = NULL;
    AddLevel(LevelArray, &TopGrid, 0);

    /* Add levels to the maximum depth or until no new levels are created,
         and re-initialize the level after it is created. */

    for (level = 0; level < MaximumRefinementLevel; level++) {
      if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
        fprintf(stderr, "Error in RebuildHierarchy.\n");
        return FAIL;
      }
      if (LevelArray[level+1] == NULL)
        break;
      LevelHierarchyEntry *Temp = LevelArray[level+1];
      while (Temp != NULL) {

         if (Temp->GridData->GalaxyDiskPatchInitializeGrid(GalaxyDiskPatchGasDistributionType,
                                               GalaxyDiskPatchGasScaleHeight_pc,
                                               GalaxyDiskPatchMidplaneDensity,
                                               GalaxyDiskPatchTotalEnergy,
                                              GalaxyDiskPatchInternalEnergy,
                                              GalaxyDiskPatchVelocity,
                                              GalaxyDiskPatchBField) == FAIL) {

            ENZO_FAIL("Error in GalaxyDiskPatchInitialize[Sub]Grid.");
        }  // end subgrid if

        Temp = Temp->NextGridThisLevel;
      }
    } // end: loop over levels

    if (UseMHD){
       MHD_ProjectB=TRUE;
       MHD_ProjectE=FALSE;
    }

    /* Loop back from the bottom, restoring the consistency among levels. */

    for (level = MaximumRefinementLevel; level > 0; level--) {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
        if (Temp->GridData->ProjectSolutionToParentGrid(
                                   *LevelArray[level-1]->GridData) == FAIL) {
          fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
          return FAIL;
        }
        Temp = Temp->NextGridThisLevel;
      }
    }

  } // end: if (RefineAtStart)

  /* set up field names and units */

  int i = 0;
  int DensNum = i;      // nicole 
  DataLabel[i++] = (char*) DensName;
  int TENum = i;       // nicole 
  DataLabel[i++] = (char*) TEName;
  if (DualEnergyFormalism)
    DataLabel[i++] = (char*) GEName;
  DataLabel[i++] = (char*) Vel1Name;
  DataLabel[i++] = (char*) Vel2Name;
  int Vel3Num = i;     // nicole 
  DataLabel[i++] = (char*) Vel3Name;
  if(CRModel)
    DataLabel[i++] =(char*) CRName;
  if( UseMHD ){
    DataLabel[i++] = BxName;
    DataLabel[i++] = ByName;
    DataLabel[i++] = BzName;
    DataLabel[i++] = PhiName;
  }

  if ( UseMHDCT ){
    MHDLabel[0] = "BxF";
    MHDLabel[1] = "ByF";
    MHDLabel[2] = "BzF";

    MHDeLabel[0] = "Ex";
    MHDeLabel[1] = "Ey";
    MHDeLabel[2] = "Ez";

    MHDUnits[0] = "None";
    MHDUnits[1] = "None";
    MHDUnits[2] = "None";

    MHDeUnits[0] = "None";
    MHDeUnits[1] = "None";
    MHDeUnits[2] = "None";
  }
  
  if (TestProblemData.UseMetallicityField){
    DataLabel[i++] = "Metal_Density";
  }

  if (UseSNColour){
    DataLabel[i++] = (char *) ColourName;
  }

  for (int j=0; j< i; j++) 
    DataUnits[j] = NULL;

  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {

    fprintf(Outfptr, "GalaxyDiskPatchGasDistributionType         = %"ISYM"\n"  , GalaxyDiskPatchGasDistributionType);
    fprintf(Outfptr, "GalaxyDiskPatchMidplaneDensity         = %"FSYM"\n"  , GalaxyDiskPatchMidplaneDensity);
    fprintf(Outfptr, "GalaxyDiskPatchGasScaleHeight_pc         = %"FSYM"\n"  , GalaxyDiskPatchGasScaleHeight_pc);
    fprintf(Outfptr, "GalaxyDiskPatchTotalEnergy        = %"FSYM"\n"  , GalaxyDiskPatchTotalEnergy);
    fprintf(Outfptr, "GalaxyDiskPatchInternalEnergy     = %"FSYM"\n"  , GalaxyDiskPatchInternalEnergy);
    fprintf(Outfptr, "TestProblemUseMetallicityField  = %"ISYM"\n", TestProblemData.UseMetallicityField);
    fprintf(Outfptr, "TestProblemInitialMetallicityFraction  = %"FSYM"\n", TestProblemData.MetallicityField_Fraction);
    
    /* nicole: add output */
    fprintf(Outfptr, "GalaxyDiskPatchUseInflow = %"ESYM"\n", UseInflow);
    
    if (UseInflow == TRUE) {
    	fprintf(Outfptr, "GalaxyDiskPatchTopDensity = %"ESYM"\n", BoundaryDensity[0]); 
    	fprintf(Outfptr, "GalaxyDiskPatchBottomDensity = %"ESYM"\n", BoundaryDensity[1]);
    	fprintf(Outfptr, "GalaxyDiskPatchTopTemperature = %"ESYM"\n", BoundaryTemperature[0]);
    	fprintf(Outfptr, "GalaxyDiskPatchBottomTemperature = %"ESYM"\n", BoundaryTemperature[1]);
    	fprintf(Outfptr, "GalaxyDiskPatchTopVelocity = %"ESYM"\n", BoundaryVelocity[0]);
    	fprintf(Outfptr, "GalaxyDiskPatchTopVelocity = %"ESYM"\n", BoundaryVelocity[1]);
    }
   /* end nicole */
   }

  /* nicole addition convert T */
  /* Get units so we can convert temperature later. */

  float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1,
    VelocityUnits=1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits,
               MetaData.Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }
  /* end nicole */ 

  /* nicole initialize the exterior */
  if (UseInflow == TRUE) { 
    Exterior.Prepare(TopGrid.GridData); 
   
    float InflowValueTop[i], InflowValueBottom[i], Dummy[i]; 
    for (int j = 0; j < i; j++){
      InflowValueTop[j] = 0.0; 
      InflowValueBottom[j] = 0.0;
      Dummy[j] = 0.0; 
    }
    
    InflowValueTop[DensNum] = BoundaryDensity[0]; 
    InflowValueBottom[DensNum] = BoundaryDensity[1]; 
    InflowValueTop[Vel3Num] = BoundaryVelocity[0]; 
    InflowValueBottom[Vel3Num] = BoundaryVelocity[1]; 
    InflowValueTop[TENum] = BoundaryTemperature[0]/TemperatureUnits/(Gamma - 1); 
    if (HydroMethod != Zeus_Hydro) 
      InflowValueTop[TENum] += 0.5*BoundaryVelocity[0]*BoundaryVelocity[0];
    InflowValueBottom[TENum] = BoundaryTemperature[1]/TemperatureUnits/(Gamma - 1);
   
    if (HydroMethod != Zeus_Hydro)
      InflowValueBottom[TENum] += 0.5*BoundaryVelocity[1]*BoundaryVelocity[1]; 

    Exterior.InitializeExternalBoundaryFace(0, periodic, periodic, Dummy, Dummy); 
    Exterior.InitializeExternalBoundaryFace(1, periodic, periodic, Dummy, Dummy);

    if (Exterior.InitializeExternalBoundaryFace(2, inflow, inflow, InflowValueTop, InflowValueBottom) == FAIL){
      ENZO_FAIL("Error Inflow InitizlizeExternalBoundaryFace.\n"); 
    }
  }  
  /* end nicole */
 
 return SUCCESS;

}
