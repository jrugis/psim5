/*
 * global_defs.hpp
 *
 *	Created on: 26/04/2018
 *			Author: jrugis
 */

#ifndef DEFS_H_
#define DEFS_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <mpi.h>
#include <string>

//************************************************************************
// model sizes
//	NOTE: a single acinus with lumen (for now)
#define CELLS_COUNT 7

//************************************************************************
// MPI defs, macros

#define MPI_ERROR_ABORT -101
#define MPI_NODES_ABORT -102
#define ACINUS_CELL_TAG 400
#define CELL_CELL_TAG 401

// mpi shutdown on error
#define mpi_abort(err)                                                \
  {                                                                   \
    std::cout << "Program shutdown error code: " << err << std::endl; \
    MPI_Abort(MPI_COMM_WORLD, err);                                   \
  }

// mpi error checking macro
#define MPI_CHECK(call)                                      \
  if ((call) != MPI_SUCCESS) {                               \
    std::cout << "MPI error calling: " << call << std::endl; \
    mpi_abort(MPI_ERROR_ABORT);                              \
  }

//************************************************************************
// acinus to cell mpi message
//  delta time, current time, solver error
enum acinus2cell { dTime, cTime, sError, ACCOUNT };

//************************************************************************
// acinus to lumen mpi message
//  delta time, current time, apical calcium, basal calcium, cell volume
enum acinus2lumen { dlTime, clTime, aC, Bc, cV, ALCOUNT };

//************************************************************************
// cell to cell mpi message
// triangle index, value
enum cell2cell { tInd, tVal, C2CCOUNT };

//************************************************************************
// cell to cell connectivity
//  this_triangle, other_cell, other_triangle
enum cell_conn { tTri, oCell, oTri, CCONNCOUNT };

//************************************************************************
enum parameter_types { calciumParms, flowParms, PTCOUNT };

// the calcium model parameters
enum calcium_parameters {
  delT,
  totalT,
  Tstride,
  PLCsrt,
  PLCfin,
  APICALds,
  APICALdl,
  c0,
  ip0,
  ce0,
  Gamma,
  Dc,
  Dp,
  De,
  Fc,
  Fip,
  d_RyR,
  V_RyR,
  K_RyR,
  K_RyR2,
  m_RyR,
  n_RyR,
  k_beta,
  K_p,
  K_c,
  K_h,
  k_IPR,
  V_p,
  k_p,
  K_bar,
  PLCds,
  PLCdl,
  V_3K,
  V_5K,
  K_PLC,
  K3K,
  V_PLC,
  h0,
  K_tau,
  tau_max,
  g0,
  K_hRyR,
  tau,
  PCOUNT
};

//************************************************************************
// the fluid flow model parameters
enum fluid_flow_parameters {
  aNkcc1,
  a1,
  a2,
  a3,
  a4,
  r,
  alpha1,
  aNaK,
  GtNa,
  GtK,
  GCl,
  KCaCC,
  eta1,
  GK,
  KCaKC,
  eta2,
  G1,
  KNa,
  KH,
  G4,
  KCl,
  KB,
  GB,
  kn,
  kp,
  pHl,
  pHi,
  pHe,
  HCO3l,
  CO20,
  Ul,
  Cle,
  Nae,
  Ke,
  HCO3e,
  CO2e,
  Hl,
  CO2l,
  Hy,
  La,
  Lb,
  Lt,
  FPCOUNT
};

//************************************************************************
// some convenience typedefs
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> MatrixN1d;
typedef Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajorBit> MatrixN3d;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixNNd;
typedef Eigen::SparseMatrix<double> SparceMatrixd;

typedef Eigen::Matrix<int, Eigen::Dynamic, 1> MatrixN1i;
typedef Eigen::Matrix<int, Eigen::Dynamic, 2, Eigen::RowMajorBit> MatrixN2i;
typedef Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajorBit> MatrixN3i;
typedef Eigen::Matrix<int, Eigen::Dynamic, 4, Eigen::RowMajorBit> MatrixN4i;
typedef Eigen::Matrix<int, Eigen::Dynamic, CCONNCOUNT, Eigen::RowMajorBit> MatrixNCi;
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> MatrixNNi;

typedef Eigen::Vector3d Vector3d;
typedef Eigen::Vector3i Vector3i;
typedef Eigen::Vector4i Vector4i;

//************************************************************************
// cell mesh values
struct sMeshVals {
  int vertices_count;
  int surface_triangles_count;
  int tetrahedrons_count;
  MatrixN3d vertices;          // 3x coordinate
  MatrixN3i surface_triangles; // 3x vertex
  MatrixN4i tetrahedrons;      // 4x vertex
};

//************************************************************************
// thermodynamic constants
#define R 8.314462100000000
#define T 310
#define F 9.645833650000000e4
#define RTF 26.721207772435513

//************************************************************************
// other constants
const double third = 1.0 / 3.0;

//************************************************************************
//************************************************************************

#endif /* DEFS_H_ */
