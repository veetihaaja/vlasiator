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
#include <vector>

#include "../definitions.h"
#include "cuda_acc_sort_blocks.hpp"

// Extra profiling stream synchronizations?
//#define SSYNC HANDLE_ERROR( cudaStreamSynchronize(stream) );
#define SSYNC

// Ensure printing of CUDA runtime errors to console
#define CUB_STDERR
#include <cub/device/device_radix_sort.cuh>

using namespace std;
using namespace spatial_cell;

// Kernels for converting GIDs to dimension-sorted indices
__global__ void blocksID_mapped_dim0_kernel(
   const vmesh::VelocityMesh* vmesh,
   vmesh::GlobalID *blocksID_mapped,
   vmesh::LocalID *blocksLID_unsorted,
   const uint nBlocks
   // const uint refL=0; //vAMR
   // const vmesh::LocalID D0 = vmesh->getGridLength(refL)[0];
   // const vmesh::LocalID D1 = vmesh->getGridLength(refL)[1];
   // const vmesh::LocalID D2 = vmesh->getGridLength(refL)[2];
   ) {
   const int cudaBlocks = gridDim.x * gridDim.y * gridDim.z;
   const uint warpSize = blockDim.x * blockDim.y * blockDim.z;
   const int blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   for (vmesh::LocalID index=blocki*warpSize; index<nBlocks; index += cudaBlocks*warpSize) {
      const vmesh::LocalID LID = (index+ti);
      if (LID < nBlocks) {
         blocksID_mapped[LID] = vmesh->getGlobalID(LID);
         blocksLID_unsorted[LID]=LID;
      }
   }
}

__global__ void blocksID_mapped_dim1_kernel(
   const vmesh::VelocityMesh* vmesh,
   vmesh::GlobalID *blocksID_mapped,
   vmesh::LocalID *blocksLID_unsorted,
   const uint nBlocks
   ) {
   const int cudaBlocks = gridDim.x * gridDim.y * gridDim.z;
   const uint warpSize = blockDim.x * blockDim.y * blockDim.z;
   const int blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   const uint refL=0; //vAMR
   const vmesh::LocalID D0 = vmesh->getGridLength(refL)[0];
   const vmesh::LocalID D1 = vmesh->getGridLength(refL)[1];
   // const vmesh::LocalID D2 = vmesh->getGridLength(refL)[2];
   for (vmesh::LocalID index=blocki*warpSize; index<nBlocks; index += cudaBlocks*warpSize) {
      const vmesh::LocalID LID = (index+ti);
      if (LID < nBlocks) {
         const vmesh::GlobalID GID = vmesh->getGlobalID(LID);
         const vmesh::LocalID x_index = GID % D0;
         const vmesh::LocalID y_index = (GID / D0) % D1;
         blocksID_mapped[LID] = GID - (x_index + y_index*D0) + y_index + x_index * D1;
         blocksLID_unsorted[LID]=LID;
      }
   }
}

__global__ void blocksID_mapped_dim2_kernel(
   const vmesh::VelocityMesh* vmesh,
   vmesh::GlobalID *blocksID_mapped,
   vmesh::LocalID *blocksLID_unsorted,
   const uint nBlocks
   ) {
   const int cudaBlocks = gridDim.x * gridDim.y * gridDim.z;
   const uint warpSize = blockDim.x * blockDim.y * blockDim.z;
   const int blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   const uint refL=0; //vAMR
   const vmesh::LocalID D0 = vmesh->getGridLength(refL)[0];
   const vmesh::LocalID D1 = vmesh->getGridLength(refL)[1];
   const vmesh::LocalID D2 = vmesh->getGridLength(refL)[2];
   for (vmesh::LocalID index=blocki*warpSize; index<nBlocks; index += cudaBlocks*warpSize) {
      const vmesh::LocalID LID = (index+ti);
      if (LID < nBlocks) {
         const vmesh::GlobalID GID = vmesh->getGlobalID(LID);
         const vmesh::LocalID x_index = GID % D0;
         const vmesh::LocalID y_index = (GID / D0) % D1;
         const vmesh::LocalID z_index = (GID / (D0*D1));
         blocksID_mapped[LID] = z_index + y_index*D2 + x_index*D1*D2;
         blocksLID_unsorted[LID]=LID;
      }
   }
}

// LIDs are already in order.
// Now also order GIDS. (can be ridiculously parallel, minus memory access patterns)
__global__ void order_GIDs_kernel(
   const vmesh::VelocityMesh* vmesh,
   vmesh::GlobalID *blocksLID,
   vmesh::GlobalID *blocksGID,
   const uint nBlocks
   ) {
   const int cudaBlocks = gridDim.x * gridDim.y * gridDim.z;
   const uint warpSize = blockDim.x * blockDim.y * blockDim.z;
   const int blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   for (vmesh::LocalID index=blocki*warpSize; index<nBlocks; index += cudaBlocks*warpSize) {
      const vmesh::LocalID i = (index+ti);
      if (i < nBlocks) {
         blocksGID[i]=vmesh->getGlobalID(blocksLID[i]);
      }
   }
}

// Kernel for scanning columnsets for block counts
__global__ void scan_blocks_for_columns_kernel(
   const vmesh::VelocityMesh* vmesh,
   const uint dimension,
   vmesh::GlobalID *blocksID_mapped_sorted,
   vmesh::LocalID *dev_columnNBlocks,
   const uint nBlocks
   ) {
   const int cudaBlocks = gridDim.x * gridDim.y * gridDim.z;
   const uint warpSize = blockDim.x * blockDim.y * blockDim.z;
   const int blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   const uint refL=0; //vAMR
   vmesh::LocalID DX;
   switch (dimension) {
      case 0:
         DX = vmesh->getGridLength(refL)[0];
         break;
      case 1:
         DX = vmesh->getGridLength(refL)[1];
         break;
      case 2:
         DX = vmesh->getGridLength(refL)[2];
         break;
      default:
         printf("Incorrect dimension in __FILE__ __LINE__\n");
   }
   for (vmesh::LocalID LID=blocki*warpSize; LID<nBlocks; LID += cudaBlocks*warpSize) {
      if (LID+ti < nBlocks) {
         vmesh::LocalID column_id = blocksID_mapped_sorted[LID+ti] / DX;
         // Increment number of blocks in column
         const vmesh::LocalID old  = atomicAdd(&dev_columnNBlocks[column_id],1);
         // // Evaluate smallest GID in column
         // old = atomicMin(&columnMinBlock[columnid],GID);
         // // Evaluate largest GID in colum
         // old = atomicMax(&columnMaxBlock[columnid],GID);
      }
   }
}

/*** Kernel for constructing columns
 Checks if all blocks in a columnset belong to a single kernel, and
 can quickly jump through the whole column.
 For columnsets containing several columns, it trials blocks by scanning
 warpSize GIDs at a time.

 Still probably room for memory optimization.
**/

__global__ void construct_columns_kernel(
   const vmesh::VelocityMesh* vmesh,
   const uint dimension,
   vmesh::GlobalID *blocksID_mapped_sorted,
   vmesh::LocalID *dev_columnNBlocks,
   ColumnOffsets* columnData,
   const uint nBlocks
   ) {
   const int cudaBlocks = gridDim.x * gridDim.y * gridDim.z;
   const uint warpSize = blockDim.x * blockDim.y * blockDim.z;
   //const int blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   if (cudaBlocks!=1) {
      printf("Error in construct_columns_kernel; unsafe gridDim\n");
      return;
   }
   const uint refL=0; //vAMR
   vmesh::LocalID DX;
   switch (dimension) {
      case 0:
         DX = vmesh->getGridLength(refL)[0];
         break;
      case 1:
         DX = vmesh->getGridLength(refL)[1];
         break;
      case 2:
         DX = vmesh->getGridLength(refL)[2];
         break;
      default:
         printf("Incorrect dimension in __FILE__ __LINE__\n");
   }
   vmesh::LocalID prev_column_id, prev_dimension_id;

   __shared__ vmesh::LocalID i;
   __shared__ vmesh::LocalID blocks_in_columnset;
   if (ti==0) {
      i = 0;
      blocks_in_columnset = 0;
      // Put in the sorted blocks, and also compute column offsets and lengths:
      columnData->columnBlockOffsets.device_push_back(0); //first offset
      columnData->setColumnOffsets.device_push_back(0); //first offset
   }
   __syncthreads();
   while (i < nBlocks) {
      // identifies a particular column
      const vmesh::LocalID column_id = blocksID_mapped_sorted[i] / DX;
      // identifies a particular block in a column (along the dimension)
      const vmesh::LocalID dimension_id = blocksID_mapped_sorted[i] % DX;
      // How many blocks in this (new) column(set)?
      if ((ti==0) && (blocks_in_columnset==0)) {
         blocks_in_columnset = dev_columnNBlocks[column_id];
      }
      // Trial: new column?
      if ( (ti==0) && (i > 0) &&  ( (column_id != prev_column_id) || (dimension_id != (prev_dimension_id + 1) ))) {
         //encountered new column! For i=0, we already entered the correct offset (0).
         //We also identify it as a new column if there is a break in the column (e.g., gap between two populations)
         //add offset where the next column will begin
         columnData->columnBlockOffsets.device_push_back(i);
         //add length of the current column that now ended
         columnData->columnNumBlocks.device_push_back(columnData->columnBlockOffsets[columnData->columnBlockOffsets.size()-1] - columnData->columnBlockOffsets[columnData->columnBlockOffsets.size()-2]);

         if (column_id != prev_column_id ){
            //encountered new set of columns, add offset to new set starting at present column
            columnData->setColumnOffsets.device_push_back(columnData->columnBlockOffsets.size() - 1);
            //add length of the previous column set that ended
            columnData->setNumColumns.device_push_back(columnData->setColumnOffsets[columnData->setColumnOffsets.size()-1] - columnData->setColumnOffsets[columnData->setColumnOffsets.size()-2]);
         }
      }
      __syncthreads();
      // Trial if only one column in columnset?
      if ( ( (blocksID_mapped_sorted[i+blocks_in_columnset-1] % DX) == (dimension_id + blocks_in_columnset - 1) ) &&
           ( (blocksID_mapped_sorted[i+blocks_in_columnset-1] / DX) == column_id ) ) {
         // skip to end of column
         if (ti==0) {
            i += blocks_in_columnset;
            blocks_in_columnset = 0;
         }
         // push_back to vectors happens at next loop
         __syncthreads();
      } else {
         // More than one column in columnset
         vmesh::LocalID this_col_length = 0;
         // Now trial by warpSize to see where column ends
         for (vmesh::LocalID ci=0; ci<blocks_in_columnset; ci += warpSize) {
            int notInColumn = 1;
            if (ci+ti < blocks_in_columnset) {
               // This evaluates if the block at the target point is no longer within the same column
               if ( (blocksID_mapped_sorted[i+ci+ti] % DX) == (dimension_id + ci+ti) &&
                    ( (blocksID_mapped_sorted[i+ci+ti] / DX) == column_id ) ) {
                  notInColumn = 0;
               }
            }
            // Warp vote to find first index (potentially) outside old column
            unsigned ballot_result = __ballot_sync(FULL_MASK, notInColumn);
            vmesh::LocalID minstep = __ffs(ballot_result); // Find first significant
            if (minstep==0) {
               minstep +=32; // no value found, jump whole warpSize
            } else {
               minstep--; // give actual index
            }
            this_col_length += minstep;
            // Exit this for loop if we reached the end of a column within a set
            if (minstep!=warpSize) {
               if (ti==0) {
                  // skip to end of column
                  i += this_col_length;
                  // Decrease number of free blocks in column(set)
                  blocks_in_columnset -= this_col_length;
               }
               __syncthreads();
               // push_back to vectors happens at next loop
               break;
            }
         }
      }
      prev_column_id = column_id;
      prev_dimension_id = dimension_id;
      __syncthreads();
   }
   // Add offsets for final column
   if (ti==0) {
      columnData->columnNumBlocks.device_push_back(nBlocks - columnData->columnBlockOffsets[columnData->columnBlockOffsets.size()-1]);
      columnData->setNumColumns.device_push_back(columnData->columnNumBlocks.size() - columnData->setColumnOffsets[columnData->setColumnOffsets.size()-1]);
   }
}

/*
   This function returns a sorted list of blocks in a cell.

   The sorted list is sorted according to the location, along the given dimension.
   This version uses triplets internally and also returns the LIDs of the sorted blocks.
*/
void sortBlocklistByDimension( //const spatial_cell::SpatialCell* spatial_cell,
                               vmesh::VelocityMesh* vmesh,
                               const vmesh::LocalID nBlocks,
                               const uint dimension,
                               vmesh::GlobalID *blocksID_mapped,
                               vmesh::GlobalID *blocksID_mapped_sorted,
                               vmesh::GlobalID *blocksGID,
                               vmesh::LocalID *blocksLID_unsorted,
                               vmesh::LocalID *blocksLID,
                               vmesh::LocalID *dev_columnNBlocks,
                               ColumnOffsets* columnData,
   // split::SplitVector<uint> columnBlockOffsets; // indexes where columns start (in blocks, length totalColumns)
   // split::SplitVector<uint> columnNumBlocks; // length of column (in blocks, length totalColumns)
   // split::SplitVector<uint> setColumnOffsets; // index from columnBlockOffsets where new set of columns starts (length nColumnSets)
   // split::SplitVector<uint> setNumColumns; // how many columns in set of columns (length nColumnSets)
                               const uint cuda_async_queue_id,
                               cudaStream_t stream
   ) {

   columnData->columnBlockOffsets.clear();
   columnData->columnNumBlocks.clear();
   columnData->setColumnOffsets.clear();
   columnData->setNumColumns.clear();

   // Ensure at least one launch block
   uint nCudaBlocks  = (nBlocks/CUDATHREADS) > CUDABLOCKS ? CUDABLOCKS : std::ceil((Real)nBlocks/(Real)CUDATHREADS);

   phiprof::start("calc new dimension id");
   // Map blocks to new dimensionality
   switch( dimension ) {
      case 0: {
         blocksID_mapped_dim0_kernel<<<nCudaBlocks, CUDATHREADS, 0, stream>>> (
            vmesh,
            blocksID_mapped,
            blocksLID_unsorted,
            nBlocks
            );
         break;
      }
      case 1: {
         blocksID_mapped_dim1_kernel<<<nCudaBlocks, CUDATHREADS, 0, stream>>> (
            vmesh,
            blocksID_mapped,
            blocksLID_unsorted,
            nBlocks
            );
         break;
      }
      case 2: {
         blocksID_mapped_dim2_kernel<<<nCudaBlocks, CUDATHREADS, 0, stream>>> (
            vmesh,
            blocksID_mapped,
            blocksLID_unsorted,
            nBlocks
            );
         break;
      }
      default:
         printf("Incorrect dimension in cuda_acc_sort_blocks.cpp\n");
   }
   HANDLE_ERROR( cudaPeekAtLastError() );
   SSYNC
   phiprof::stop("calc new dimension id");

   phiprof::start("CUB sort");
   // Sort (with CUB)

   // Determine temporary device storage requirements
   void     *dev_temp_storage = NULL;
   size_t   temp_storage_bytes = 0;
   cub::DeviceRadixSort::SortPairs(dev_temp_storage, temp_storage_bytes,
                                   blocksID_mapped, blocksID_mapped_sorted,
                                   blocksLID_unsorted, blocksLID, nBlocks,
                                   0, sizeof(vmesh::GlobalID)*8, stream);
   HANDLE_ERROR( cudaPeekAtLastError() );
   phiprof::start("cudamallocasync");
   HANDLE_ERROR( cudaMallocAsync((void**)&dev_temp_storage, temp_storage_bytes, stream) );
   SSYNC
   phiprof::stop("cudamallocasync");
   //printf("allocated %d bytes of temporary memory for CUB SortPairs\n",temp_storage_bytes);

   // Now sort
   cub::DeviceRadixSort::SortPairs(dev_temp_storage, temp_storage_bytes,
                                   blocksID_mapped, blocksID_mapped_sorted,
                                   blocksLID_unsorted, blocksLID, nBlocks,
                                   0, sizeof(vmesh::GlobalID)*8, stream);
   HANDLE_ERROR( cudaPeekAtLastError() );
   HANDLE_ERROR( cudaStreamSynchronize(stream) ); // In case SortPairs won't like the free below too soon
   HANDLE_ERROR( cudaFreeAsync(dev_temp_storage, stream) );
   phiprof::stop("CUB sort");

   // Gather GIDs in order
   phiprof::start("reorder GIDs");
   order_GIDs_kernel<<<nCudaBlocks, CUDATHREADS, 0, stream>>> (
      vmesh,
      blocksLID,
      blocksGID,
      nBlocks
      );
   HANDLE_ERROR( cudaPeekAtLastError() );
   phiprof::stop("reorder GIDs");

   phiprof::start("Scan for column block counts");
   scan_blocks_for_columns_kernel<<<nCudaBlocks, CUDATHREADS, 0, stream>>> (
      vmesh,
      dimension,
      blocksID_mapped,
      dev_columnNBlocks,
      nBlocks
      );
   HANDLE_ERROR( cudaPeekAtLastError() );
   phiprof::stop("Scan for column block counts");

   phiprof::start("construct columns");
   // Construct columns. To ensure order,
   // these are done serially, but still form within a kernel.
   construct_columns_kernel<<<1, CUDATHREADS, 0, stream>>> (
      vmesh,
      dimension,
      blocksID_mapped_sorted,
      dev_columnNBlocks,
      columnData,
      nBlocks
      );
   HANDLE_ERROR( cudaPeekAtLastError() );
   SSYNC
   phiprof::stop("construct columns");
   // printf("\n Output for dimension %d ",dimension);
   // printf("\nColumnBlockOffsets %d\n", columnData->columnBlockOffsets.size());
   // //for (auto i : columnData->columnBlockOffsets) printf("%d ",i);
   // printf("\ncolumnNumBlocks %d\n", columnData->columnNumBlocks.size());
   // //for (auto i : columnData->columnNumBlocks) printf("%d ",i);
   // printf("\nsetColumnOffsets %d\n", columnData->setColumnOffsets.size());
   // //for (auto i : columnData->setColumnOffsets) printf("%d ",i);
   // printf("\nsetNumColumns %d\n", columnData->setNumColumns.size());
   // //for (auto i : columnData->setNumColumns) printf("%d ",i);
   // printf("\n \n",dimension);

}
