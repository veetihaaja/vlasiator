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

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <numbers>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../object_wrapper.h"

#include "LCCP_reconnection.h"

using namespace std;

namespace projects {
   LCCP_Reconnection::LCCP_Reconnection(): Project() { }
   LCCP_Reconnection::~LCCP_Reconnection() { }

   bool LCCP_Reconnection::initialize(void) {
      bool success = Project::initialize();

      return success;
   }

   void LCCP_Reconnection::addParameters() {
      typedef Readparameters RP;
      RP::add<Real>("LCCP_Reconnection.B0", "B0", this->B0,1.0e-10);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;
         LCCP_ReconnectionSpeciesParameters* sP=new LCCP_ReconnectionSpeciesParameters();

         this->speciesParamsRead.push_back(sP);
         RP::add<Real>(pop + "_LCCP_Reconnection.rho", "Number density (m^-3)", sP->rho,10.e8);
         RP::add<Real>(pop + "_LCCP_Reconnection.Temperature", "Temperature (K)", sP->T,0.86456498092);
         if (pop == "proton") {
            this->ionT = sP->T;
            this->ionMass = getObjectWrapper().particleSpecies[i].mass;
            this->ionRho = sP->rho;
         }
         if (pop == "electron") {
            this->electronT = sP->T;
            this->electronMass = getObjectWrapper().particleSpecies[i].mass;
         }

      }
      // now that we have values, we can calculate stuff
      this->ionInertialLength = physicalconstants::LIGHT_SPEED / sqrt(this->ionRho * physicalconstants::CHARGE * physicalconstants::CHARGE / (physicalconstants::EPS_0 * this->ionMass));
   }

   void LCCP_Reconnection::getParameters(){
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
        this->speciesParams.push_back(*this->speciesParamsRead.at(i));
      }
   }

   // Realf LCCP_Reconnection::MaxwellianPhaseSpaceDensity_LCCP_reconnection(
   //    creal& vx, creal& vy, creal& vz, creal& v_ts, creal& rho) const {
   //    return (rho / pow((sqrt(2.0 * numbers::pi) * v_ts), 3.0)) * exp(- (vx*vx + vy*vy + vz*vz)/(2.0 * v_ts * v_ts));
   // }

   Realf LCCP_Reconnection::fillPhaseSpace(spatial_cell::SpatialCell *cell,
                                       const uint popID,
                                       const uint nRequested
      ) const {
      const LCCP_ReconnectionSpeciesParameters& sP = this->speciesParams[popID];

      // Fetch spatial cell center coordinates
      const Real x  = cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX];
      const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      // const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];

      creal mass = getObjectWrapper().particleSpecies[popID].mass;
      creal charge = getObjectWrapper().particleSpecies[popID].charge;
      creal mu0 = physicalconstants::MU_0;

      creal rho = sP.rho * (1.0/cosh(y/(this->ionInertialLength*0.5))) * (1.0/cosh(y/(this->ionInertialLength*0.5)));

      Real Theta = sP.T / (this->ionT + this->electronT);
      Real v_ts = sqrt((Theta * this->B0 * this->B0)/ (2.0 * physicalconstants::MU_0 * sP.rho * mass));
      Real v_ds = -1.0 * 2.0 * mass * v_ts * v_ts / (charge * this->B0 * this->ionInertialLength / 2.0);

      #ifdef USE_GPU
      vmesh::VelocityMesh *vmesh = cell->dev_get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* VBC = cell->dev_get_velocity_blocks(popID);
      #else
      vmesh::VelocityMesh *vmesh = cell->get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* VBC = cell->get_velocity_blocks(popID);
      #endif
      // Loop over blocks
      Realf rhosum = 0;
      arch::parallel_reduce<arch::null>(
         {WID, WID, WID, nRequested},
         ARCH_LOOP_LAMBDA (const uint i, const uint j, const uint k, const uint initIndex, Realf *lsum ) {
            vmesh::GlobalID *GIDlist = vmesh->getGrid()->data();
            Realf* bufferData = VBC->getData();
            const vmesh::GlobalID blockGID = GIDlist[initIndex];
            // Calculate parameters for new block
            Real blockCoords[6];
            vmesh->getBlockInfo(blockGID,&blockCoords[0]);
            creal vxBlock = blockCoords[0];
            creal vyBlock = blockCoords[1];
            creal vzBlock = blockCoords[2];
            creal dvxCell = blockCoords[3];
            creal dvyCell = blockCoords[4];
            creal dvzCell = blockCoords[5];
            ARCH_INNER_BODY(i, j, k, initIndex, lsum) {
               creal vx = vxBlock + (i+0.5)*dvxCell;
               creal vy = vyBlock + (j+0.5)*dvyCell;
               creal vz = vzBlock + (k+0.5)*dvzCell - v_ds;
               //const Realf value = MaxwellianPhaseSpaceDensity(vx, vy, vz, sP.T, rho, mass);
               const Realf value = MaxwellianPhaseSpaceDensity_LCCP_reconnection(vx, vy, vz, v_ts, rho);
               bufferData[initIndex*WID3 + k*WID2 + j*WID + i] = value;
               //lsum[0] += value;
            };
         }, rhosum);
      return rhosum;
   }

   void LCCP_Reconnection::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
      //Real* cellParams = cell->get_cell_parameters();
      //creal x = cellParams[CellParams::XCRD];
      //creal dx = cellParams[CellParams::DX];
      //creal y = cellParams[CellParams::YCRD];
      //creal dy = cellParams[CellParams::DY];
      //
      //Real ksi = ((x + 0.5 * dx)  * cos(this->ALPHA) + (y + 0.5 * dy) * sin(this->ALPHA)) / this->WAVELENGTH;
      //Real dBxavg = sin(2.0 * M_PI * ksi);
      //Real dByavg = sin(2.0 * M_PI * ksi);
      //Real dBzavg = cos(2.0 * M_PI * ksi);
   }

   void LCCP_Reconnection::setProjectBField(fsgrids::perbspan perb,
                                 fsgrids::bgbspan bgb,
                                 fsgrids::technicalspan technical, FieldSolverGrid &fsgrid) {
      setBackgroundFieldToZero(fsgrid, technical, bgb);

      if (!P::isRestart) {
         // local copies for lambda capture
         const Real ionInertialLength_lcl = this->ionInertialLength;
         const Real B0_lcl = this->B0;

         fsgrid.parallel_for([](int timerId) -> phiprof::Timer { return phiprof::Timer{timerId}; },
                             phiprof::initializeTimer("setProjectBField"), technical,
                             [=](const fsgrid::Coordinates &coordinates, const fsgrid::FsStencil& stencil, cuint sysBoundaryFlag, cuint sysBoundaryLayer) {
            const std::array<Real, 3> xyz = coordinates.getPhysicalCoords(stencil.i, stencil.j, stencil.k);
            const std::array<Real, 3> gridSpacing = coordinates.physicalGridSpacing;
            auto& cell = perb[stencil.ooo()];

            const Real dx = gridSpacing[0];
            const Real dy = gridSpacing[1];

            cell[fsgrids::bfield::PERBX] = B0_lcl * tanh(xyz[1] / (ionInertialLength_lcl * 0.5)); // Bx for the reconnection test
            cell[fsgrids::bfield::PERBY] = 0.0;
            cell[fsgrids::bfield::PERBZ] = 0.0;
         });
      }
   }

} // namespace projects
