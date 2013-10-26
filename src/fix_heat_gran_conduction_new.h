/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(heat/gran/conduction_new,FixHeatGranCondNew)
FixStyle(heat/gran,FixHeatGranCondNew)

#else

#ifndef LMP_FIX_HEATGRAN_CONDUCTION_NEW_H
#define LMP_FIX_HEATGRAN_CONDUCTION_NEW_H

#include "fix_heat_gran.h"

namespace LAMMPS_NS {

  class FixHeatGranCondNew : public FixHeatGran {
  public:
    FixHeatGranCondNew(class LAMMPS *, int, char **);
    ~FixHeatGranCondNew();
    void pre_delete(bool);

    int setmask();
    void init();
    void post_force(int);

    void cpl_evaluate(class ComputePairGranLocal *);
    void register_compute_pair_local(ComputePairGranLocal *);
    void unregister_compute_pair_local(ComputePairGranLocal *);

  private:
    template <int> void post_force_eval(int,int);

    class FixPropertyGlobal* fix_conductivity;
    double *conductivity;

    // for heat transfer area correction
    int area_correction_flag;
    double const* const* deltan_ratio;
  };

}

#endif
#endif

