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

/* ----------------------------------------------------------------------
   Contributing authors:
   Andrei Malinouski
------------------------------------------------------------------------- */

#include "fix_heat_gran_conduction_new.h"

#include "atom.h"
#include "compute_pair_gran_local.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"
#include "force.h"
#include "math_extra.h"
#include "mech_param_gran.h"
#include "modify.h"
#include "neigh_list.h"
#include "pair_gran.h"
#include "stdlib.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixHeatGranCondNew::FixHeatGranCondNew(class LAMMPS *lmp, int narg, char **arg) : FixHeatGran(lmp, narg, arg)
{
  int iarg = 5;

  area_correction_flag = 0;
  cond_fluid = 0.; //default value - non-conducting fluid
  bool hasargs = true;
  while(iarg < narg && hasargs)
  {
    hasargs = false;
    if (strcmp(arg[iarg],"area_correction") == 0)
    {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments for keyword 'area_correction'");
      if(strcmp(arg[iarg+1],"yes") == 0)
        area_correction_flag = 1;
      else if(strcmp(arg[iarg+1],"no") == 0)
        area_correction_flag = 0;
      else error->fix_error(FLERR,this,"");
      iarg += 2;
      hasargs = true;
    } 
    else if (strcmp(arg[iarg],"fluid_thermal_conduction") == 0)
    {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      cond_fluid = atof(arg[iarg+1]);
      if(cond_fluid < 0. ) error->fix_error(FLERR,this,"invalid fluid_thermal_conduction");
      iarg += 2;
      hasargs = true;
    }
   else if(strcmp(style,"heat/gran/conduction_new") == 0)
        error->fix_error(FLERR,this,"unknown keyword or wrong keyword order");
  }

  fix_conductivity = NULL;
  conductivity = NULL;

  fix_heatCond = fix_heatConv = fix_heatRadiat = fix_heatExtern = NULL;
}

/* ---------------------------------------------------------------------- */

FixHeatGranCondNew::~FixHeatGranCondNew()
{

  if (conductivity)
    delete []conductivity;
}

/* ---------------------------------------------------------------------- */

void FixHeatGranCondNew::post_create()
{
  // register new fluxes
  if(!fix_heatCond)
  {
    char* fixarg[11];
    fixarg[0]="conductiveHeatFlux";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="conductiveHeatFlux";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.";
    fixarg[9]="0.";
    fixarg[10]="0.";
    fix_heatCond = modify->add_fix_property_atom(11,fixarg,style);
  }
  if(!fix_heatConv)
  {
    char* fixarg[11];
    fixarg[0]="convectiveHeatFlux";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="convectiveHeatFlux";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.";
    fixarg[9]="0.";
    fixarg[10]="0.";
    fix_heatConv = modify->add_fix_property_atom(11,fixarg,style);
  }
  if(!fix_heatRadiat)
  {
    char* fixarg[11];
    fixarg[0]="radiativeHeatFlux";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="radiativeHeatFlux";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.";
    fixarg[9]="0.";
    fixarg[10]="0.";
    fix_heatRadiat = modify->add_fix_property_atom(11,fixarg,style);
  }
  if(!fix_heatExtern)
  {
    char* fixarg[11];
    fixarg[0]="externalHeatFlux";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="externalHeatFlux";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.";
    fixarg[9]="0.";
    fixarg[10]="0.";
    fix_heatExtern = modify->add_fix_property_atom(11,fixarg,style);
  }
}

/* ---------------------------------------------------------------------- */

void FixHeatGranCondNew::altern_updatePtrs()
{

  Temp = fix_temp->vector_atom;
  vector_atom = Temp; 

  heatFlux = fix_heatFlux->vector_atom;
  heatSource = fix_heatSource->vector_atom;
  directionalHeatFlux = fix_directionalHeatFlux->array_atom;

  conductiveHeatFlux = fix_heatCond->vector_atom;
  convectiveHeatFlux = fix_heatConv->vector_atom;
  radiativeHeatFlux = fix_heatRadiat->vector_atom;
  externalHeatFlux = fix_heatExtern->vector_atom;

}

/* ---------------------------------------------------------------------- */

void FixHeatGranCondNew::pre_delete(bool unfixflag)
{

  // tell cpl that this fix is deleted
  if(cpl && unfixflag) cpl->reference_deleted();

}

/* ---------------------------------------------------------------------- */

int FixHeatGranCondNew::setmask()
{
  int mask = FixHeatGran::setmask();
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixHeatGranCondNew::init()
{
  const double *Y, *nu, *Y_orig;
  double expo, Yeff_ij, Yeff_orig_ij, ratio;
  Fix *ymo_fix;

  if (FHG_init_flag == false){
    FixHeatGran::init();
  }

  fix_heatCond = static_cast<FixPropertyAtom*>(modify->find_fix_property("conductiveHeatFlux","property/atom","scalar",0,0,style));
  fix_heatConv = static_cast<FixPropertyAtom*>(modify->find_fix_property("convectiveHeatFlux","property/atom","scalar",0,0,style));
  fix_heatRadiat = static_cast<FixPropertyAtom*>(modify->find_fix_property("radiativeHeatFlux","property/atom","scalar",0,0,style)); 
  fix_heatExtern = static_cast<FixPropertyAtom*>(modify->find_fix_property("externalHeatFlux","property/atom","scalar",0,0,style));
  altern_updatePtrs();

  int max_type = pair_gran->mpg->max_type();
  if (conductivity) delete []conductivity;
  conductivity = new double[max_type];
  fix_conductivity = static_cast<FixPropertyGlobal*>(modify->find_fix_property("thermalConductivity","property/global","peratomtype",max_type,0,style));

  // pre-calculate conductivity for possible contact material combinations
  for(int i=1;i< max_type+1; i++)
      for(int j=1;j<max_type+1;j++)
      {
          conductivity[i-1] = fix_conductivity->compute_vector(i-1);
          if(conductivity[i-1] < 0.) error->all(FLERR,"Fix heat/gran/conduction_new: Thermal conductivity must not be < 0");
      }

  // calculate heat transfer correction

  ymo_fix = NULL;
  if(area_correction_flag)
  {
    ymo_fix = modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",0,0,style);

    if(force->pair_match("gran/hooke",0)) expo = 1.;
    else if(force->pair_match("gran/hertz",0)) expo = 2./3.;
    else error->fix_error(FLERR,this,"area correction could not identify the granular pair style you are using, supported are hooke and hertz types");

    Y = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0,style))->get_values();
    nu = static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0,style))->get_values();
    Y_orig = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style))->get_values();

    // allocate a new array within youngsModulusOriginal
    static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style))->new_array(max_type,max_type);

    // feed deltan_ratio into this array
    for(int i = 1; i < max_type+1; i++)
    {
      for(int j = 1; j < max_type+1; j++)
      {
        Yeff_ij      = 1./((1.-pow(nu[i-1],2.))/Y[i-1]     +(1.-pow(nu[j-1],2.))/Y[j-1]);
        Yeff_orig_ij = 1./((1.-pow(nu[i-1],2.))/Y_orig[i-1]+(1.-pow(nu[j-1],2.))/Y_orig[j-1]);
        ratio = pow(Yeff_ij/Yeff_orig_ij,expo);
        
        static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style))->array_modify(i-1,j-1,ratio);
      }
    }

    // get reference to deltan_ratio
    deltan_ratio = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style))->get_array_modified();
  }

  altern_updatePtrs();

  // error checks on coarsegraining
  if(force->cg_active())
    error->cg(FLERR,this->style);
}

/* ---------------------------------------------------------------------- */

void FixHeatGranCondNew::post_force(int vflag){

  //template function for using touchflag or not
  if(history_flag == 0) post_force_eval<0>(vflag,0);
  if(history_flag == 1) post_force_eval<1>(vflag,0);

}

/* ---------------------------------------------------------------------- */

void FixHeatGranCondNew::cpl_evaluate(ComputePairGranLocal *caller)
{
  if(caller != cpl) error->all(FLERR,"Illegal situation in FixHeatGranCondNew::cpl_evaluate");
  if(history_flag == 0) post_force_eval<0>(0,1);
  if(history_flag == 1) post_force_eval<1>(0,1);
}

/* ---------------------------------------------------------------------- */

template <int HISTFLAG>
void FixHeatGranCondNew::post_force_eval(int vflag,int cpl_flag)
{
  double hc1,hc2,contactArea,delta_n,flux,dirFlux[3];
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv,tcoi,tcoj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double  displ, alpha, beta, cosTeta0, cosTetaC,rij;//

  int newton_pair = force->newton_pair;

  if (strcmp(force->pair_style,"hybrid")==0)
    error->warning(FLERR,"Fix heat/gran/conduction_new implementation may not be valid for pair style hybrid");
  if (strcmp(force->pair_style,"hybrid/overlay")==0)
    error->warning(FLERR,"Fix heat/gran/conduction_new implementation may not be valid for pair style hybrid/overlay");

  inum = pair_gran->list->inum;
  ilist = pair_gran->list->ilist;
  numneigh = pair_gran->list->numneigh;
  firstneigh = pair_gran->list->firstneigh;
  if(HISTFLAG) firsttouch = pair_gran->listgranhistory->firstneigh;

  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  altern_updatePtrs();

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    if(HISTFLAG) touch = firsttouch[i];

    conductiveHeatFlux[i] = 0.;
    radiativeHeatFlux[i] = 0.;
    convectiveHeatFlux[i] = 0.;//resetting from the previous step. heatFlux is set in ste module. This step  may be moved there later
    for (jj = 0; jj < jnum; jj++)
      {
      j = jlist[jj];
      j &= NEIGHMASK;

      if (!(mask[i] & groupbit) && !(mask[j] & groupbit)) continue;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        radj = radius[j];
        radsum = radi + radj;
        r = sqrt(rsq);
	displ = 0.5*(r - radsum);
	rij = radi; //////////////////////////////////////////solid angle?
        tcoi = conductivity[type[i]-1];
    	tcoj = conductivity[type[j]-1];
	alpha = (tcoi*tcoj*cond_fluid)*(4*radi*radj)/(radsum* (tcoi*tcoj+tcoj*cond_fluid+tcoi*cond_fluid) );//with effective radius for polydisperse
	beta = cond_fluid*((2*radi*radj)/(radsum)+displ);
	cosTeta0 = ((2*radi*radj)/(radsum)+displ)/sqrt(((2*radi*radj)/(radsum)+displ)*((2*radi*radj)/(radsum)+displ)+rij*rij);

      if (rsq < radsum*radsum)
	     {  //contact
 		   if(area_correction_flag)
    		    {
    		      delta_n = radsum - r;
    		      delta_n *= deltan_ratio[type[i]-1][type[j]-1];
   		      r = radsum - delta_n;
    		    }
	
    		    contactArea = - M_PI/4 * ( (r-radi-radj)*(r+radi-radj)*(r-radi+radj)*(r+radi+radj) )/(r*r); //contact area of the two spheres

     		   if (tcoi < SMALL || tcoj < SMALL) hc1 = 0.;
     		   else hc1 = 4.*tcoi*tcoj/(tcoi+tcoj)*sqrt(contactArea);//through contact spot

   		   conductiveHeatFlux[i]+= (Temp[j]-Temp[i])*hc1;
		   cosTetaC = ((2*radi*radj)/(radsum)+displ)/sqrt(((2*radi*radj)/(radsum)+displ)*((2*radi*radj)/(radsum)+displ)+contactArea);
		   hc2 = M_PI*beta*log((beta-alpha*cosTeta0)/(beta-alpha*cosTetaC));//conduction through fluid layer
		   convectiveHeatFlux[i] += (Temp[j]-Temp[i])*hc2;
    	     }
	     else//no direct contact; only through gas layer
	     {
		hc2 = M_PI*beta*log((beta-alpha*cosTeta0)/(beta-alpha));
		conductiveHeatFlux[i] += 0.;
		convectiveHeatFlux[i] += (Temp[j]-Temp[i])*hc2;
	      }
		flux = (Temp[j]-Temp[i])*(hc1+hc2);
	             dirFlux[0] = flux*delx;
   		     dirFlux[1] = flux*dely;
   		     dirFlux[2] = flux*delz;
    		    if(!cpl_flag)
   		     {
   		     	  //Add half of the flux (located at the contact) to each particle in contact
    		    	  heatFlux[i] += flux;
    		    	  directionalHeatFlux[i][0] += 0.50 * dirFlux[0];
     		   	  directionalHeatFlux[i][1] += 0.50 * dirFlux[1];
     		   	  directionalHeatFlux[i][2] += 0.50 * dirFlux[2];
    		   	   if (newton_pair || j < nlocal)
    		   	   {
        		    heatFlux[j] -= flux;
   		     	    directionalHeatFlux[j][0] += 0.50 * dirFlux[0];
   		   	    directionalHeatFlux[j][1] += 0.50 * dirFlux[1];
  		            directionalHeatFlux[j][2] += 0.50 * dirFlux[2];
    		           }

      		     }

        	     if(cpl_flag && cpl) cpl->add_heat(i,j,flux);
     }//end of jj cycle
  }//end of ii cycle

  if(newton_pair) fix_heatFlux->do_reverse_comm();
  if(newton_pair) fix_directionalHeatFlux->do_reverse_comm();
}

/* ----------------------------------------------------------------------
   register and unregister callback to compute
------------------------------------------------------------------------- */

void FixHeatGranCondNew::register_compute_pair_local(ComputePairGranLocal *ptr)
{
   
   if(cpl != NULL)
      error->all(FLERR,"Fix heat/gran/conduction_new allows only one compute of type pair/local");
   cpl = ptr;
}

void FixHeatGranCondNew::unregister_compute_pair_local(ComputePairGranLocal *ptr)
{
   
   if(cpl != ptr)
       error->all(FLERR,"Illegal situation in FixHeatGranCondNew::unregister_compute_pair_local");
   cpl = NULL;
}
