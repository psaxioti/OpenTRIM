/*********************************************************************

    Copyright 2019, Christian Borschel

    This file is part of iradina.

    iradina is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    iradina is distributed WITHOUT ANY WARRANTY; without even the implied
    warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with iradina.  If not, see <http://www.gnu.org/licenses/>.

***********************************************************************/




/*****************************************************************************/
/* This module can be specialized to incorporate special geometries          */
/*                                                                           */
/* This is the NANOWIRE version of the geometry module                       */
/*                                                                           */
/*****************************************************************************/


#include "geometry.h"

int InitSpecialGeometry(){
  /* iradina calls this function if the special-geometry-parameter is activated.
     It can be used to load additional parameters from other files, if desired */

  int result;
  char filename[100];
  
  strcpy(filename,"nw_structure.in");

  /* Load additional nanowire structure file */
  result=IniFileReader(NW_DataBlockReader, NW_DataReader, filename);
  if(result!=0){
    message_error(result,"Error reading NW structure file %s.\n",filename);
    return result;
  }
  
  /* Calculate some more parameters */
  NW_radius     = NW_diameter*0.5;
  NW_radius_sqr = NW_radius * NW_radius;
  x_half        = target_size_x/2.0; /* Necessary to place NW in the middle */
  y_half        = target_size_y/2.0;

  return 0;
}

int special_CheckSurfaceAtom(int cellindex, float x, float y, float z, float spacing){
  /* like GetMaterial, but compare to a reduced radius */
  float temp;
  temp = (NW_radius-spacing);
  temp = temp * temp;
  /*  if( ((x-NW_radius)*(x-NW_radius)+(y-NW_radius)*(y-NW_radius)) < temp ) {*/
  if( ((x-x_half)*(x-x_half)+(y-y_half)*(y-y_half)) < temp ) {
    return 0;
  } else {
    return 1;
  }
}

int GetMaterialFromPosition(int cell, float x, float y, float z){
  /* if( ((x-NW_radius)*(x-NW_radius)+(y-NW_radius)*(y-NW_radius)) < NW_radius_sqr ) {*/ /* we're inside the NW */
  if( ((x-x_half)*(x-x_half)+(y-y_half)*(y-y_half)) < NW_radius_sqr ) { /* we're inside the NW */
    return 0;
  } else {
    return 1;
  }
}

float special_CalcSurfaceNormal(int new_cell, int cell,float x, float y, float z, float* nx, float* ny, float* nz, float factor){
  /* if the ion moves from material or outside, we need surface normal */
  /* Factor can be used to change its direction */
  float temp;
  /*  x -= NW_radius;
      y -= NW_radius;*/
  x -= x_half;
  y -= y_half;
  temp = factor/sqrtf(x*x + y*y); /* inverse factorized length */
  *nz=0;
  *nx=x*temp;
  *ny=y*temp;
  return 0;
}

float special_CalcDirectionalFractionSqr(float x, float y, float z, float vx, float vy, float vz){
  /* This function returns the square directional fraction of the energy which is perpendicular to the surface.
     This function should return 0, if the velocity does not point into vacuum but into the solid. This way we
     know if we must apply the E_surf or E_latt */
  float temp;
  
  /*  x -= NW_radius;
      y -= NW_radius;*/
  x -= x_half;
  y -= y_half;
  temp = x*vx+y*vy;
  if(temp<0){
    return 0;
  } else {
    return (vx*vx+vy*vy)/(x*x+y*y)*temp*temp;
  }
}


int NW_DataBlockReader(char* BlockName){
  return 0;
}

int NW_DataReader(char* ParName, char* ParValue){
  /* reads the data from the NW structure input file */
  if(strcmp(ParName,"NW_diameter")==0){ /* The NW diameter */
    sscanf(ParValue,"%f",&NW_diameter);
    message(1,"NW diameter:\t %g nm\n",NW_diameter);
  }
  return 0;
}
