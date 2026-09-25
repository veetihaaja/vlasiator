/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifndef LCCP_RECONNECTION_H
#define LCCP_RECONNECTION_H

#include "../../definitions.h"
#include "../project.h"

namespace projects {

   struct LCCP_ReconnectionSpeciesParameters {
      Real rho;
      Real T;
   };

   class LCCP_Reconnection: public Project {
    public:
      LCCP_Reconnection();
      virtual ~LCCP_Reconnection();
      
      virtual bool initialize(void) override;
      virtual void addParameters(void) override;
      virtual void getParameters(void) override;
      virtual void setProjectBField(
         fsgrids::perbspan perb,
         fsgrids::bgbspan bgb,
         fsgrids::technicalspan technical, FieldSolverGrid &fsgrid
      ) override;
      //ARCH_HOSTDEV inline Realf MaxwellianPhaseSpaceDensity_LCCP_reconnection(
      //   creal& vx, creal& vy, creal& vz, creal& v_ts, creal& rho) const;

      virtual Realf fillPhaseSpace(spatial_cell::SpatialCell *cell,
                                  const uint popID,
                                  const uint nRequested) const override;
      virtual void calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) override;
      
      Real B0;
      Real ionInertialLength;
      Real electronMass; // these are here globally for reasons
      Real ionMass; 
      Real electronT;
      Real ionT;
      Real ionRho;
      Real Theta;
      std::vector<LCCP_ReconnectionSpeciesParameters> speciesParams;
      std::vector<LCCP_ReconnectionSpeciesParameters*> speciesParamsRead;
   } ; // class LCCP_Reconnection
} // namespace projects

#endif
