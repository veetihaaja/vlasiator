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

#include <algorithm>
#include <cmath>
#include <utility>

#include <Eigen/Geometry>
#include <Eigen/Core>

#include "cpu_acc_semilag.hpp"
#include "cpu_acc_transform.hpp"
#include "cpu_acc_intersections.hpp"
#include "cpu_acc_map.hpp"

using namespace std;
using namespace spatial_cell;
using namespace Eigen;

/*!
  Prepare to accelerate species in cell. Sets the maximum allowed dt to the
  correct value.

 * @param spatial_cell Spatial cell containing the accelerated population.
 * @param popID ID of the accelerated particle species.
*/

void prepareAccelerateCell(
   SpatialCell* spatial_cell,
   const uint popID){   
   updateAccelerationMaxdt(spatial_cell, popID);
}

/*!
  Compute the number of subcycles needed for the acceleration of the particle
  species in the spatial cell. Note that one should first prepare to
  accelerate the cell with prepareAccelerateCell.

 * @param spatial_cell Spatial cell containing the accelerated population.
 * @param popID ID of the accelerated particle species.
*/

uint getAccelerationSubcycles(SpatialCell* spatial_cell, Real dt, const uint popID)
{
   //return max( convert<uint>(ceil(dt*spatial_cell->CellParams[CELLPARAMS::TIMECLASSDT] / spatial_cell->get_max_v_dt(popID))), 1u);
   return max( convert<uint>(ceil(dt / spatial_cell->get_max_v_dt(popID))), 1u);
}

/*!
  Propagates the distribution function in velocity space of given real
  space cell.

  Based on SLICE-3D algorithm: Zerroukat, M., and T. Allen. "A
  three‐dimensional monotone and conservative semi‐Lagrangian scheme
  (SLICE‐3D) for transport problems." Quarterly Journal of the Royal
  Meteorological Society 138.667 (2012): 1640-1651.

 * @param spatial_cell Spatial cell containing the accelerated population.
 * @param popID ID of the accelerated particle species.
 * @param vmesh Velocity mesh.
 * @param blockContainer Velocity block data container.
 * @param map_order Order in which vx,vy,vz mappings are performed. 
 * @param dt Time step of one subcycle.
*/

void cpu_accelerate_cell(SpatialCell* spatial_cell,
                         const uint popID,     
                         const uint map_order,
                         const Real& dt,
                         int tc_delta) {
   //double t1 = MPI_Wtime();
   vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>* vmeshPtr = NULL; 
   vmesh::VelocityBlockContainer<vmesh::LocalID>* blockContainerPtr = NULL;
   // Main branch: check if time-ghost data required, handle that and recurse as needed
   // .... of course this is now substepping unnecessarily...
   // std::cout << spatial_cell->parameters[CellParams::CELLID] << "c Accelerate, initial refs" << " vmesh " << &vmesh << " blockContainer " << &blockContainer << "\n";

   if(tc_delta == 0)
   {
      for(auto i : spatial_cell->requested_timeclass_ghosts) {
         if(i - spatial_cell->get_tc() > 0) {
            spatial_cell->set_velocity_mesh_ghost(popID, i);
            spatial_cell->set_velocity_blocks_ghost(popID, i);
            int tc_d = i-spatial_cell->get_tc();
            cpu_accelerate_cell(spatial_cell, popID, map_order, dt/pow(2,tc_d), tc_d);
         }
         else if (i - spatial_cell->get_tc() < 0){
            spatial_cell->set_velocity_mesh_ghost(popID, i);
            spatial_cell->set_velocity_blocks_ghost(popID, i);
            int tc_d = i-spatial_cell->get_tc();
            // cpu_accelerate_cell(spatial_cell, popID, map_order, dt/pow(2,tc_d), tc_d);
         }
      }
      vmeshPtr = &spatial_cell->get_velocity_mesh(popID);
      blockContainerPtr = &spatial_cell->get_velocity_blocks(popID);
   }
   // Ghost branch: select the ghost vmesh instead for acc
   // dt we already have adjusted as needed
   else{
      std::cout << "Request at tc " << spatial_cell->get_tc() + tc_delta << " at dt = " << dt << " being run at cell " << 
      spatial_cell->parameters[CellParams::CELLID] << "\n";
      vmeshPtr = &spatial_cell->get_velocity_mesh_ghost(popID,spatial_cell->get_tc() + tc_delta);
      blockContainerPtr = &spatial_cell->get_velocity_blocks_ghost(popID,spatial_cell->get_tc() + tc_delta);
    }
    
   vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh = *vmeshPtr;
   vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = *blockContainerPtr;

   std::cout << spatial_cell->parameters[CellParams::CELLID] << "c Accelerate tc " << spatial_cell->get_tc() << " with tcdelta " 
   << tc_delta << " at dt = " << dt <<"; vmesh " << &vmesh << " blockContainer " << &blockContainer << "\n";

   // vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = spatial_cell->get_velocity_blocks(popID);

   // compute transform, forward in time and backward in time
   phiprof::Timer transformTimer {"compute-transform"};

   //compute the transform performed in this acceleration
   Transform<Real,3,Affine> fwd_transform= compute_acceleration_transformation(spatial_cell,popID,dt);
   Transform<Real,3,Affine> bwd_transform= fwd_transform.inverse();
   transformTimer.stop();

   // if (spatial_cell->parameters[CellParams::CELLID] == 16)
   // #pragma omp critical(output)
   // {
   //    std::cout<< "Cellid " << spatial_cell->parameters[CellParams::CELLID] << 
   //    " at t=" << spatial_cell->parameters[CellParams::TIME_V] <<": dt " << 
   //    dt <<  "\n";
   //    std::stringstream ss;
   //    ss.precision(64);
   //    ss << std::scientific;
   //    ss << fwd_transform(0,0) << " " << fwd_transform(0,1) << " " << fwd_transform(0,2) << "\n" <<
   //          fwd_transform(1,0) << " " << fwd_transform(1,1) << " " << fwd_transform(1,2) << "\n" <<
   //          fwd_transform(2,0) << " " << fwd_transform(2,1) << " " << fwd_transform(2,2) << "\n";
   //    std::cout << ss.str() << "\n";
   // }


   const uint8_t refLevel = 0;
   Real intersection_z,intersection_z_di,intersection_z_dj,intersection_z_dk;
   Real intersection_x,intersection_x_di,intersection_x_dj,intersection_x_dk;
   Real intersection_y,intersection_y_di,intersection_y_dj,intersection_y_dk;
   int intersections_id {phiprof::initializeTimer("compute-intersections")};
   int mapping_id {phiprof::initializeTimer("compute-mapping")};
   switch(map_order){
      case 0: {
         //Map order XYZ
         phiprof::Timer intersectionsTimer {intersections_id};
         compute_intersections_1st(vmesh,bwd_transform, fwd_transform, 0, refLevel,
                                   intersection_x,intersection_x_di,intersection_x_dj,intersection_x_dk);
         compute_intersections_2nd(vmesh,bwd_transform, fwd_transform, 1, refLevel,
                                   intersection_y,intersection_y_di,intersection_y_dj,intersection_y_dk);
         compute_intersections_3rd(vmesh,bwd_transform, fwd_transform, 2, refLevel,
                                   intersection_z,intersection_z_di,intersection_z_dj,intersection_z_dk);
         intersectionsTimer.stop();
         phiprof::Timer mappingTimer {mapping_id};
         map_1d(spatial_cell, popID, intersection_x,intersection_x_di,intersection_x_dj,intersection_x_dk,0,vmesh,blockContainer); // map along x
         map_1d(spatial_cell, popID, intersection_y,intersection_y_di,intersection_y_dj,intersection_y_dk,1,vmesh,blockContainer); // map along y
         map_1d(spatial_cell, popID, intersection_z,intersection_z_di,intersection_z_dj,intersection_z_dk,2,vmesh,blockContainer); // map along z
         mappingTimer.stop();
         break;
      }
         
      case 1: {
         //Map order YZX
         phiprof::Timer intersectionsTimer{intersections_id};
         compute_intersections_1st(vmesh, bwd_transform, fwd_transform, 1, refLevel,
                                   intersection_y,intersection_y_di,intersection_y_dj,intersection_y_dk);
         compute_intersections_2nd(vmesh, bwd_transform, fwd_transform, 2, refLevel,
                                   intersection_z,intersection_z_di,intersection_z_dj,intersection_z_dk);
         compute_intersections_3rd(vmesh, bwd_transform, fwd_transform, 0, refLevel,
                                   intersection_x,intersection_x_di,intersection_x_dj,intersection_x_dk);
      
         intersectionsTimer.stop();
         phiprof::Timer mappingTimer {mapping_id};
         map_1d(spatial_cell, popID, intersection_y,intersection_y_di,intersection_y_dj,intersection_y_dk,1,vmesh,blockContainer); // map along y
         map_1d(spatial_cell, popID, intersection_z,intersection_z_di,intersection_z_dj,intersection_z_dk,2,vmesh,blockContainer); // map along z
         map_1d(spatial_cell, popID, intersection_x,intersection_x_di,intersection_x_dj,intersection_x_dk,0,vmesh,blockContainer); // map along x
         mappingTimer.stop();
         break;
      }

      case 2: {
         phiprof::Timer intersectionsTimer{intersections_id};
         //Map order Z X Y
         compute_intersections_1st(vmesh, bwd_transform, fwd_transform, 2, refLevel,
                                   intersection_z,intersection_z_di,intersection_z_dj,intersection_z_dk);
         compute_intersections_2nd(vmesh, bwd_transform, fwd_transform, 0, refLevel,
                                   intersection_x,intersection_x_di,intersection_x_dj,intersection_x_dk);
         compute_intersections_3rd(vmesh, bwd_transform, fwd_transform, 1, refLevel,
                                   intersection_y,intersection_y_di,intersection_y_dj,intersection_y_dk);
         intersectionsTimer.stop();
         phiprof::Timer mappingTimer {mapping_id};
         map_1d(spatial_cell, popID, intersection_z,intersection_z_di,intersection_z_dj,intersection_z_dk,2,vmesh,blockContainer); // map along z
         map_1d(spatial_cell, popID, intersection_x,intersection_x_di,intersection_x_dj,intersection_x_dk,0,vmesh,blockContainer); // map along x
         map_1d(spatial_cell, popID, intersection_y,intersection_y_di,intersection_y_dj,intersection_y_dk,1,vmesh,blockContainer); // map along y
         mappingTimer.stop();
         break;
      }
   }

   if (Parameters::prepareForRebalance == true) {
//       spatial_cell->parameters[CellParams::LBWEIGHTCOUNTER] += (MPI_Wtime() - t1);
   }
}
