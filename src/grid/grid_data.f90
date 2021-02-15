#include "../util/macro.f90.inc"

module grid_data_m

  use error_m
  use grid_m, only: grid

  implicit none

  private
  public allocate_grid_data
  public grid_data_int,  grid_data_log,  grid_data_real,  grid_data_cmplx
  public grid_data1_int, grid_data1_log, grid_data1_real, grid_data1_cmplx
  public grid_data2_int, grid_data2_log, grid_data2_real, grid_data2_cmplx
  public grid_data3_int, grid_data3_log, grid_data3_real, grid_data3_cmplx
  public grid_data4_int, grid_data4_log, grid_data4_real, grid_data4_cmplx
  public grid_data5_int, grid_data5_log, grid_data5_real, grid_data5_cmplx
  public grid_data6_int, grid_data6_log, grid_data6_real, grid_data6_cmplx
  public grid_data7_int, grid_data7_log, grid_data7_real, grid_data7_cmplx
  public grid_data8_int, grid_data8_log, grid_data8_real, grid_data8_cmplx

#define T int
#define TT integer
#include "grid_data_def.f90.inc"

#define T log
#define TT logical
#define TLOG
#include "grid_data_def.f90.inc"

#define T real
#define TT real
#include "grid_data_def.f90.inc"

#define T cmplx
#define TT complex
#include "grid_data_def.f90.inc"

contains

#define T int
#define TT integer
#include "grid_data_imp.f90.inc"

#define T log
#define TT logical
#define TLOG
#include "grid_data_imp.f90.inc"

#define T real
#define TT real
#include "grid_data_imp.f90.inc"

#define T cmplx
#define TT complex
#include "grid_data_imp.f90.inc"

end module
