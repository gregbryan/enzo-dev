/***********************************************************************
/
/  CHECK FOR TIME ACTION
/
/  written by: Greg Bryan
/  date:       September, 2000
/  modified1:
/
/  PURPOSE:
/    This routine checks to see if the time has arrived for a given
/      "time-action" to occur (i.e. some action that is applied to
/       the data at a given time/redshift).
/
************************************************************************/
#include <stdlib.h> 
#include <math.h>
#include <stdio.h>
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
#include "LevelHierarchy.h"
#include "CosmologyParameters.h"
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, double *MassUnits, FLOAT Time);

int CommunicationBroadcastValues(Eflt64 *Values, int Number, int BroadcastProcessor);
 
int CheckForTimeAction(LevelHierarchyEntry *LevelArray[],
		       TopGridData &MetaData)
{
 
  /* Declarations. */
 
  int i, level;
 
  /* Check for time against list of actions. */
 
  for (i = 0; i < MAX_TIME_ACTIONS; i++)
    if (MetaData.Time >= TimeActionTime[i] && TimeActionTime[i] > 0) {
 
      if (debug)
	printf("Applying TimeAction %"ISYM" at t=%"GOUTSYM"\n", i, TimeActionTime[i]);

      float DensityUnits, LengthUnits,TemperatureUnits, TimeUnits,VelocityUnits;
      double MassUnits;
      if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                &TimeUnits, &VelocityUnits, &MassUnits, InitialTimeInCodeUnits) == FAIL) {
        ENZO_FAIL("Error in GetUnits in CheckForTimeAction.C.");
      }

      /* If Type 2 (SN explosion), then handle seperately. */
	 
      if (TimeActionType[i] == 2) {

	TimeActionTime[i] += TimeActionParameter[i]*3.15e7/TimeUnits;  // re-occur

	/* First come up with location of SNR */

	FLOAT SNPosition[3]={0.5,0.5,0.5};
        FLOAT DomainWidth[3];
        int dim;
        float SNScaleHeight_CodeUnits;
        SNScaleHeight_CodeUnits = SNScaleHeight_pc*pc_cm/LengthUnits;

        for (dim=0;dim<3;dim++){
          DomainWidth[dim] = DomainRightEdge[dim]-DomainLeftEdge[dim];
	}
	if (MyProcessorNumber == ROOT_PROCESSOR) {
 
          if (SNDistribution_Miao == 1){

              SNPosition[0] = DomainWidth[0]* rand()*1.0/RAND_MAX;    
              SNPosition[1] = DomainWidth[1]* rand()*1.0/RAND_MAX;

              float SNIaScaleHeight_CodeUnits, SNIILowScaleHeight_CodeUnits,SNIIHighScaleHeight_CodeUnits, SN_zpos,ran1,ran2;
              int sign;
              SNIaScaleHeight_CodeUnits = SNIaScaleHeight_pc*pc_cm/LengthUnits;
              SNIILowScaleHeight_CodeUnits = SNIILowScaleHeight_pc*pc_cm/LengthUnits;
              SNIIHighScaleHeight_CodeUnits = SNIIHighScaleHeight_pc*pc_cm/LengthUnits;

              /*decide if this is a Type Ia or Type II SN*/

              if (rand()*1.0/RAND_MAX < SNIaFraction) {
                 printf("This is Ia\n");
                /*Exponential distribution of Type Ia SN*/
                 sign = ( rand()*1.0/RAND_MAX>0.5 ) ? 1: -1 ; 
                 SN_zpos = sign*  SNIaScaleHeight_CodeUnits* log( rand()*1.0/RAND_MAX ) ;
	      }  //end of SNIa 

              else{
                 printf("This is II");
		 /*Gaussian distribution of Type II SN*/
		 ran1 = rand()*1.0/RAND_MAX;
		 ran2 = rand()*1.0/RAND_MAX;
                    
		 if (rand()*1.0/RAND_MAX < SNIILowFraction ) { //if this type II SN is in the low scale height group
		   SN_zpos = SNIILowScaleHeight_CodeUnits * pow(-2 *log(ran1),0.5) * cos(2*3.1415926535*ran2); //Box-Muller Method to generate random numbers with normal distribution 
		 } // end of rand()*1.0/RAND_MAX < SNIILowFraction

		 else {
		   SN_zpos = SNIIHighScaleHeight_CodeUnits * pow(-2 *log(ran1),0.5) * cos(2*3.1415926535*ran2);
		 } // end of ELSE of (rand()*1.0/RAND_MAX < SNIILowFraction)

	      }// end of Type II

              SNPosition[2] = SN_zpos;

	  } // end of SNDistribution_Miao == 1

          else { // a uniformly random position within a cubic box
          SNPosition[0] = DomainWidth[0]* rand()*1.0/RAND_MAX;	  
          SNPosition[1] = DomainWidth[1]* rand()*1.0/RAND_MAX;	  

// following 4 lines: if SN only confined to a fraction of x-y plane
//          float x = 300.*pc_cm /LengthUnits;
//          SNPosition[0] = 0.5*DomainWidth[0] -0.5 *x+ x * rand()*1.0/RAND_MAX; 
//          SNPosition[1] = 0.5*DomainWidth[1] -0.5 *x+ x * rand()*1.0/RAND_MAX; 

          SNPosition[2] = SNScaleHeight_CodeUnits * 2*(rand()*1.0/RAND_MAX-0.5);

// For fixed position:
//          SNPosition[0] = 0.05;
//          SNPosition[1] = 0.05;
//          SNPosition[2] = 0.00;

	  } // end of ELSE of (  SNDistribution_Miao == 1 )

	}
       
	CommunicationBroadcastValues(SNPosition, 3, ROOT_PROCESSOR);
        if(MyProcessorNumber == ROOT_PROCESSOR)
	  printf("SNPosition random numbers =\n9999	%f	%f	%f\n",SNPosition[0],SNPosition[1],SNPosition[2]);

        float SNRadius;
        SNRadius = SingleSNRadius_pc*pc_cm/LengthUnits;

	/* Now apply to all grids. */
 
	for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
	  LevelHierarchyEntry *Temp = LevelArray[level];
	  while (Temp != NULL) {
	    
	    if (Temp->GridData->ApplySNExplosion_Miao(SNPosition, SNRadius)  == FAIL) {
	      ENZO_FAIL("Error in grid->ApplyExplosion\n");

	    }
	    
	    Temp = Temp->NextGridThisLevel;
	  }
	}

	 /* Loop back from the bottom, restoring the consistency among levels. */

	for (level = MaximumRefinementLevel; level > 0; level--) {
	  LevelHierarchyEntry *Temp = LevelArray[level];
	  while (Temp != NULL) {
	    if (Temp->GridData->ProjectSolutionToParentGrid(
				   *Temp->GridHierarchyEntry->ParentGrid->GridData) == FAIL) {
	      fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
	      return FAIL;
	    }
	    Temp = Temp->NextGridThisLevel;
	  }
	}

// just comment out    }


      } else {
 
	/* Done, turn it off (-1 in redshift indicates off). */
 
	TimeActionTime[i] = 0;
	TimeActionRedshift[i] = -1;
 
	/* Now apply to all grids. */
 
	for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
//	for (level = MaximumRefinementLevel; level < MAX_DEPTH_OF_HIERARCHY; level++) { //to add SN only at finest level
	  LevelHierarchyEntry *Temp = LevelArray[level];
	  while (Temp != NULL) {
	    if (Temp->GridData->ApplyTimeAction(TimeActionType[i],
						TimeActionParameter[i]) == FAIL) {
	      ENZO_FAIL("Errot in grid->ApplyTimeAction\n");

	    }
	    Temp = Temp->NextGridThisLevel;
	  }
	}
      } // end: ifelse TimeActionType == 2
 
    }
 
  return SUCCESS;
}		
