

int Check_CellPBC(int i, int iMax)
{
  int j;
  
  if (i > (iMax-1))
    j = i - iMax;
  else if (i < 0)
    j = i + iMax;
  else
    j = i;

  return j;
}


void Check_NeighborsInCell (int n0, int i)
{
  int p, Idx, k;
  double dx[3], DX[3], dr;
  
  for (p = 0; p < (Cell[n0].AtomIdxMax); p++)
    {
      Idx = Cell[n0].AtomIdx[p];

      if (Idx > i)
        {
          /* clatoms.pos, DX is in Bohr */
          V3V3subV3(&clatoms.pos[3*i], &clatoms.pos[3*Idx], DX);

          Cart2Red(DX, invH, dx);   /* invH in 1/Bohr, dx in reduced units */

          dx[0] = checkPBC_Shifted (dx[0]);
          dx[1] = checkPBC_Shifted (dx[1]);
          dx[2] = checkPBC_Shifted (dx[2]);

          Red2Cart(dx, H, DX); /* H in Bohr */

          DX[0] *= BOHR_ANGS;  /* now DX in Angstroms */
          DX[1] *= BOHR_ANGS;  /* now DX in Angstroms */
          DX[2] *= BOHR_ANGS;  /* now DX in Angstroms */

          dr = NormV3(DX);     /* dr in Angstroms */
          
          if (dr < NeighborList_Cutoff_Global)
            {
              k = AtomsNeighborList[i].NeighbourMax;

              AtomsNeighborList[i].NList[k] = Idx;

              AtomsNeighborList[i].NDist[3*k+0] = DX[0];
              AtomsNeighborList[i].NDist[3*k+1] = DX[1];
              AtomsNeighborList[i].NDist[3*k+2] = DX[2];

              AtomsNeighborList[i].NeighbourMax++;
            }
        }
      if (AtomsNeighborList[i].NeighbourMax >= MAX_Neighbour_Count)
        printf("\nSERIOUS PROBLEM : Need to increase the MAX_Neighbour_Count\n\n");
    }
  
  return ;
}

void InitializeCell()
{
  /*
    Cells have dimensions of Nx x Ny x Nz, where Nx = H[0][0]/CutOFF_Global, Ny = H[1][1]/CutOFF_Global, etc.
    Total number of cells Ncells = det(H)/(CutOFF_Global*CutOFF_Global*CutOFF_Global)
    Each cell can have an average of (Natoms/NCells) atoms
  */

  int i, k, nx, ny, nz, n0, mx, my, mz, m0, ix, iy, iz;
  double dx[3];

  //  NeighborList_Cutoff_Global = 1.10*CutOFF_Global;   /* Angstroms */


  /* Define Cell parameters, allocate memory */
  CellNx = (int) floor(H[0][0]/(ANGS_BOHR*CutOFF_Global*1.10)) + 1;  /* Hmatrix is in Bohr, CutOFF_Global is in Angstroms */
  CellNy = (int) floor(H[1][1]/(ANGS_BOHR*CutOFF_Global*1.10)) + 1;  /* Hmatrix is in Bohr, CutOFF_Global is in Angstroms */
  CellNz = (int) floor(H[2][2]/(ANGS_BOHR*CutOFF_Global*1.10)) + 1;  /* Hmatrix is in Bohr, CutOFF_Global is in Angstroms */

  //  printf("%f %f %f %f\n", H[0][0], H[1][1], H[2][2], (ANGS_BOHR*CutOFF_Global*1.10));

  CellX = 1.00/((double) CellNx);
  CellY = 1.00/((double) CellNy);
  CellZ = 1.00/((double) CellNz);
  
  Ncells = CellNx*CellNy*CellNz;
  Natoms_per_cell = (int) (3.0*Natoms)/(1.0*Ncells);
  
  printf("The cell sizes are %d x %d x %d and total cells are %d\n", CellNx, CellNy, CellNz, Ncells);
  printf("Dimensions of each cell (reduced) : (%f %f %f), atoms per cell %d\n", CellX, CellY, CellZ, Natoms_per_cell);

  Cell = (Global_Cell *) malloc (Ncells*sizeof(Global_Cell));
  for (i = 0; i < Ncells; i++)
    {
      Cell[i].AtomIdx = (int *) malloc (Natoms_per_cell*sizeof(int));
      Cell[i].AtomIdxMax = 0;
    }

  for (i = 0; i < Natoms; i++)
    {
      AtomsNeighborList[i].NeighbourMax = 0;

      Cart2Red(&clatoms.pos[3*i], invH, dx);   /* invH in 1/Bohr, dx in reduced units */
      
      dx[0] = checkPBC(dx[0]);
      dx[1] = checkPBC(dx[1]);
      dx[2] = checkPBC(dx[2]);
      
      nx = (int) floor(dx[0]/CellX);
      ny = (int) floor(dx[1]/CellY);
      nz = (int) floor(dx[2]/CellZ);

      n0 = nx + ny*CellNx + nz*CellNx*CellNy;

      //      printf("%d %d (%f %f %f)\n", n0, Cell[n0].AtomIdxMax, clatoms.pos[3*i+0], clatoms.pos[3*i+1], clatoms.pos[3*i+2]);

      Cell[n0].AtomIdx[Cell[n0].AtomIdxMax++] = i;
    }

  for (i = 0; i < Natoms; i++)
    {

      Cart2Red(&clatoms.pos[3*i], invH, dx);   /* invH in 1/Bohr, dx in reduced units */
      
      dx[0] = checkPBC(dx[0]);
      dx[1] = checkPBC(dx[1]);
      dx[2] = checkPBC(dx[2]);

      nx = (int) floor(dx[0]/CellX);
      ny = (int) floor(dx[1]/CellY);
      nz = (int) floor(dx[2]/CellZ);

      n0 = nx + ny*CellNx + nz*CellNx*CellNy;

      k = 0;
      for (ix = -1; ix < 2; ix++)
        for (iy = -1; iy < 2; iy++)
          for (iz = -1; iz < 2; iz++)
            {
              /* shift along x,y,z-axis */
              mx = Check_CellPBC(nx+ix, CellNx);
              my = Check_CellPBC(ny+iy, CellNy);
              mz = Check_CellPBC(nz+iz, CellNz);

              m0 = mx + my*CellNx + mz*CellNx*CellNy;

              Check_NeighborsInCell(m0, i);

              Cell[n0].CellNeighborIdx[k] = m0;
              k++;
            }
    }

  /*
  printf("finished memory allocation and cell initialization, %d\n", AtomsNeighborList[11].NeighbourMax);
  for(i = 0; i < AtomsNeighborList[11].NeighbourMax; i++)
    printf("%d  %f  %f  %f\n", AtomsNeighborList[11].NList[i], AtomsNeighborList[11].NDist[3*i+0], 
       AtomsNeighborList[11].NDist[3*i+1], AtomsNeighborList[11].NDist[3*i+2]);
  */
  
  return;
}


void UpdateCell()
{
  int i, nx, ny, nz, n0, mx, my, mz, m0, ix, iy, iz;
  double dx[3];

  for (i = 0; i < Ncells; i++)
    Cell[i].AtomIdxMax = 0;

  for (i = 0; i < Natoms; i++)
    {
      AtomsNeighborList[i].NeighbourMax = 0;

      Cart2Red(&clatoms.pos[3*i], invH, dx);   /* invH in 1/Bohr, dx in reduced units */
      
      dx[0] = checkPBC(dx[0]);
      dx[1] = checkPBC(dx[1]);
      dx[2] = checkPBC(dx[2]);
      
      nx = (int) floor(dx[0]/CellX);
      ny = (int) floor(dx[1]/CellY);
      nz = (int) floor(dx[2]/CellZ);

      n0 = nx + ny*CellNx + nz*CellNx*CellNy;

      Cell[n0].AtomIdx[Cell[n0].AtomIdxMax++] = i;
    }

  for (i = 0; i < Natoms; i++)
    {

      Cart2Red(&clatoms.pos[3*i], invH, dx);   /* invH in 1/Bohr, dx in reduced units */
      
      dx[0] = checkPBC(dx[0]);
      dx[1] = checkPBC(dx[1]);
      dx[2] = checkPBC(dx[2]);

      nx = (int) floor(dx[0]/CellX);
      ny = (int) floor(dx[1]/CellY);
      nz = (int) floor(dx[2]/CellZ);

      n0 = nx + ny*CellNx + nz*CellNx*CellNy;

      for (ix = -1; ix < 2; ix++)
        for (iy = -1; iy < 2; iy++)
          for (iz = -1; iz < 2; iz++)
            {
              /* shift along x,y,z-axis */
              mx = Check_CellPBC(nx+ix, CellNx);
              my = Check_CellPBC(ny+iy, CellNy);
              mz = Check_CellPBC(nz+iz, CellNz);

              m0 = mx + my*CellNx + mz*CellNx*CellNy;

              Check_NeighborsInCell(m0, i);
            }
    }

  /*
  printf("finished memory allocation and cell initialization, %d\n", AtomsNeighborList[11].NeighbourMax);
  for(i = 0; i < AtomsNeighborList[11].NeighbourMax; i++)
    printf("%d  %f  %f  %f\n", AtomsNeighborList[11].NList[i], AtomsNeighborList[11].NDist[3*i+0], 
    AtomsNeighborList[11].NDist[3*i+1], AtomsNeighborList[11].NDist[3*i+2]);
  */

  return;
}



void Create_NeighborList()
{
  int i, j;
  double dx[3], dr, DX[3];

  for (i = 0; i < Natoms; i++)
    AtomsNeighborList[i].NeighbourMax = 0;
  
  for (i = 0; i < (Natoms-1); i++)
    {
      for (j = (i+1); j < Natoms; j++)
        {
	  /* clatoms.pos, DX is in Bohr */
          V3V3subV3(&clatoms.pos[3*i], &clatoms.pos[3*j], DX);

          Cart2Red(DX, invH, dx);   /* invH in 1/Bohr, dx in reduced units */

          dx[0] = checkPBC_Shifted (dx[0]);
          dx[1] = checkPBC_Shifted (dx[1]);
          dx[2] = checkPBC_Shifted (dx[2]);

          Red2Cart(dx, H, DX); /* H in Bohr */

          DX[0] *= BOHR_ANGS;  /* now DX in Angstroms */
          DX[1] *= BOHR_ANGS;  /* now DX in Angstroms */
          DX[2] *= BOHR_ANGS;  /* now DX in Angstroms */

	  dr = NormV3(DX);     /* dr in Angstroms */
	  
	  if (dr < NeighborList_Cutoff_Global)
	    {
	      AtomsNeighborList[i].NList[AtomsNeighborList[i].NeighbourMax] = j;
	      AtomsNeighborList[i].NDist[3*AtomsNeighborList[i].NeighbourMax + 0] = DX[0];
	      AtomsNeighborList[i].NDist[3*AtomsNeighborList[i].NeighbourMax + 1] = DX[1];
	      AtomsNeighborList[i].NDist[3*AtomsNeighborList[i].NeighbourMax + 2] = DX[2];
	      AtomsNeighborList[i].NeighbourMax++;
	    }
	}
      if (AtomsNeighborList[i].NeighbourMax >= MAX_Neighbour_Count)
	printf("\nSERIOUS PROBLEM : Need to increase the MAX_Neighbour_Count\n\n");
    }

  return;
}


double EnergyForce_EAM()
{
  double *density = NULL, EnPair, EnEmbed, Etotal;
  double dr, Fx[3], dx[3], DX[3], F, argLoc, pvten_tmp[9];
  int i, j, Ni;

  //  update_H(par_rahman.hmat, par_rahman.hmati, &(par_rahman.vol));

  density = (double *) malloc (Natoms*sizeof(double));

  for (i = 0; i < Natoms; i++)
    density[i] = 0;

  EnPair = 0;
  for (i = 0; i < Natoms; i++)
    {
      for (Ni = 0; Ni < AtomsNeighborList[i].NeighbourMax; Ni++)
        {
	  j = AtomsNeighborList[i].NList[Ni];

	  /* clatoms.pos, DX in Bohr */
          V3V3subV3(&clatoms.pos[3*i], &clatoms.pos[3*j], DX);

          Cart2Red(DX, invH, dx); /* invH in 1/Bohr */

          dx[0] = checkPBC_Shifted (dx[0]);
          dx[1] = checkPBC_Shifted (dx[1]);
          dx[2] = checkPBC_Shifted (dx[2]);

          Red2Cart(dx, H, DX); /* DX, H in Bohr */

          DX[0] *= BOHR_ANGS;  /* now DX in Angstroms */
          DX[1] *= BOHR_ANGS;  /* now DX in Angstroms */
          DX[2] *= BOHR_ANGS;  /* now DX in Angstroms */

	  dr = NormV3(DX);     /* dr in Angstroms */

	  /* update the distance in the Neighbor list */
	  AtomsNeighborList[i].NDist[3*Ni+0] = DX[0];
	  AtomsNeighborList[i].NDist[3*Ni+1] = DX[1];
	  AtomsNeighborList[i].NDist[3*Ni+2] = DX[2];

	  if (dr < CutOFF_Global)  /* CutOFF_Global in Angstroms */
	    {
	      argLoc = gsl_spline_eval (SplineGSL.splineDen, dr, SplineGSL.acc1);
	      density[i] += argLoc;
	      density[j] += argLoc;
	      
	      EnPair += gsl_spline_eval (SplineGSL.splinePair, dr, SplineGSL.acc3);
	      /* EnPair is in eV */
	    }
        }
    }
  
  EnEmbed = 0;
  for (i = 0;i < Natoms; i++)
    EnEmbed += gsl_spline_eval (SplineGSL.splineEmbed, density[i], SplineGSL.acc2);
  /* EnEmbed is in eV */
  
  Etotal = (EnEmbed + EnPair)*EV_HARTR;  /* in Hartree */
  //  printf("Energy is %.12f, %.12f, %.12f, %.12f \n", Etotal, Etotal/((double) Natoms), EnEmbed, EnPair);

  InitV9(pvten_tmp);
  InitV9(&(ptens.pvten_Global[0]));

  for (i = 0; i < 3*Natoms; i++)
    clatoms.force[i] = 0;

  for (i = 0; i < Natoms; i++)
    {
      for (Ni = 0; Ni < AtomsNeighborList[i].NeighbourMax; Ni++)
        {
          j = AtomsNeighborList[i].NList[Ni];

	  DX[0] = AtomsNeighborList[i].NDist[3*Ni+0];
	  DX[1] = AtomsNeighborList[i].NDist[3*Ni+1];
	  DX[2] = AtomsNeighborList[i].NDist[3*Ni+2];
	  
	  dr = NormV3(DX);

	  if (dr < CutOFF_Global)
	    {
	      F = gsl_spline_eval_deriv(SplineGSL.splinePair, dr, SplineGSL.acc3);

	      F += gsl_spline_eval_deriv(SplineGSL.splineEmbed, density[i], SplineGSL.acc2)*
		gsl_spline_eval_deriv(SplineGSL.splineDen, dr, SplineGSL.acc1);

	      F += gsl_spline_eval_deriv(SplineGSL.splineEmbed, density[j], SplineGSL.acc2)*
		gsl_spline_eval_deriv(SplineGSL.splineDen, dr, SplineGSL.acc1);

	      F *= EV_HARTR/ANGS_BOHR;       /* F in eV/Ang -> Hartree/Bohr */

	      V3SmulV3(DX, (F/dr), Fx);      /* F, Fx in Hartree/Bohr */

	      V3SmulV3addV3(DX, (-F/dr), &(clatoms.force[3*i]));  /* F, Fx in Hartree/Bohr */
	      V3SmulV3addV3(DX, (F/dr), &(clatoms.force[3*j]));   /* F, Fx in Hartree/Bohr */

	      V3SmulV3(DX, (-1*ANGS_BOHR), DX); 
	      V3V3diadicV9addV9(DX, Fx, pvten_tmp);
	    }
        }
    }
  
  V9symmetrizeV9(pvten_tmp, &(ptens.pvten_Global[0]));
  free(density);

  return Etotal;
}


double ql_pot_forc_calc(COLLC_VAR *Ql)
{
  int i, j, k, Ni;
  
  int l_ordr            = steinhdt_pkg.l_ordr;
  
  double arg, arg1, f_pre, DX[3], dr, r12t, r12_inv, r12_inv3;
  double r12xy, r12xy_inv, r12xy_inv3, dx12t, dy12t, dz12t;
  double cos_th, dcos_th_x, dcos_th_y, dcos_th_z, fc, dfc_x, dfc_y, dfc_z, ql_v;
  double tmpalp[l_ordr+1], tmpalp_d[l_ordr+1], ForceLoc[3], diff, Energy_ExtVar;
  double pvten_tmp[9];
  
  double *pre_sph2      = steinhdt_pkg.pre_sph2;

  double rmin_ql        = steinhdt_pkg.rmin_ql;
  double rmax_ql        = steinhdt_pkg.rmax_ql;

  double *cos_mphi      = steinhdt_pkg.cos_mphi;
  double *sin_mphi      = steinhdt_pkg.sin_mphi;

  double *dcos_mphi_x   = steinhdt_pkg.dcos_mphi_x;
  double *dsin_mphi_x   = steinhdt_pkg.dsin_mphi_x;
  double *dcos_mphi_y   = steinhdt_pkg.dcos_mphi_y;
  double *dsin_mphi_y   = steinhdt_pkg.dsin_mphi_y;

  double *alp           = steinhdt_pkg.alp;
  double *alp_d         = steinhdt_pkg.alp_d;

  double *qlm_re        = steinhdt_pkg.qlm_re;
  double *qlm_im        = steinhdt_pkg.qlm_im;

  double *sh_re         = steinhdt_pkg.sh_re;
  double *sh_im         = steinhdt_pkg.sh_im;


  /* START : Memory allocation */
  qlm_re = (double *) malloc ((l_ordr+1)*sizeof(double));
  qlm_im = (double *) malloc ((l_ordr+1)*sizeof(double));
  pre_sph2 = (double *) malloc ((l_ordr+1)*sizeof(double));
  sh_re = (double *) malloc ((l_ordr+1)*sizeof(double));
  sh_im = (double *) malloc ((l_ordr+1)*sizeof(double));
  cos_mphi = (double *) malloc ((l_ordr+1)*sizeof(double));
  sin_mphi = (double *) malloc ((l_ordr+1)*sizeof(double));
  dcos_mphi_x = (double *) malloc ((l_ordr+1)*sizeof(double));
  dsin_mphi_x = (double *) malloc ((l_ordr+1)*sizeof(double));
  dcos_mphi_y = (double *) malloc ((l_ordr+1)*sizeof(double));
  dsin_mphi_y = (double *) malloc ((l_ordr+1)*sizeof(double));
  alp = (double *) malloc ((l_ordr+1)*sizeof(double));
  alp_d = (double *) malloc ((l_ordr+1)*sizeof(double));
  /* END : Memory allocation */

  /* START : Initialize variables */
  for (i = 0; i < 3*Natoms; i++)
    Ql_force[i] = 0;

  for(i = 0; i <= l_ordr; i++)
    {
      qlm_re[i]=0.0; 
      qlm_im[i]=0.0; 
    } 
     
  for(i = 1; i <= l_ordr; i++)
    {
      pre_sph2[i] = 1.0;
      for(j = 1; j <= 2*i; j++)
	pre_sph2[i] *= (l_ordr-i+j);

      pre_sph2[i] = 1.0/pre_sph2[i];
    }
  pre_sph2[0] = 1.0*0.5;  /* divide by 2 for the purpose to time 2 later */
  /* END : Initialize variables */  


  /* START : Ql calculation */  
  for( i = 0; i <= (Natoms-2); i++)
    { /* the last atom has an empty neighbor list */
      for (Ni = 0; Ni < AtomsNeighborList[i].NeighbourMax; Ni++)
        {
          j = AtomsNeighborList[i].NList[Ni];
	  
	  /* DX is in Angstroms */
	  DX[0] = AtomsNeighborList[i].NDist[3*Ni+0];
	  DX[1] = AtomsNeighborList[i].NDist[3*Ni+1];
	  DX[2] = AtomsNeighborList[i].NDist[3*Ni+2];

	  dr = NormV3(DX);  /* dr in Angstroms */
	  
	  dx12t = -DX[0]*ANGS_BOHR;
	  dy12t = -DX[1]*ANGS_BOHR;
	  dz12t = -DX[2]*ANGS_BOHR;

	  r12t  = dr*ANGS_BOHR;
	  r12_inv = 1.0/r12t;
	  r12_inv3 = r12_inv*r12_inv*r12_inv;

	  r12xy = sqrt(dx12t*dx12t + dy12t*dy12t);
	  r12xy_inv = 1.0/r12xy;
	  r12xy_inv3 = r12xy_inv*r12xy_inv*r12xy_inv;

	  /*  switch function */
	  if(r12t <= rmin_ql)  /* r12t, rmin_ql in BOHR */
	    {
	      fc = 1.0;
	      dfc_x = 0.0;
	      dfc_y = 0.0;
	      dfc_z = 0.0;
	    }
	  else if(r12t > rmax_ql) /* r12t, rmax_ql in BOHR */
	    {
	      fc = 0.0;
	      dfc_x = 0.0;
	      dfc_y = 0.0;
	      dfc_z = 0.0;
	    }
	  else
	    {
	      arg1 = M_PI/(rmax_ql-rmin_ql);
	      fc = 0.5*cos((r12t-rmin_ql)*arg1)+0.5;
	      arg = -0.5*arg1*sin((r12t-rmin_ql)*arg1)*r12_inv;
	      dfc_x = -arg*dx12t;
	      dfc_y = -arg*dy12t;
	      dfc_z = -arg*dz12t;
	    }

	  /* Ql calculation */
	  if(r12t < rmax_ql)
	    {
	      cos_th    = dz12t*r12_inv;
	      
	      /* cos(mPhi), m = 0, 1 */
	      cos_mphi[0] = 1.0;
	      cos_mphi[1] = dx12t*r12xy_inv;

	      /* sin(mPhi), m = 0, 1 */
	      sin_mphi[0] = 0.0;
	      sin_mphi[1] = dy12t*r12xy_inv;

	      /* derivative cos(mPhi), m = 0, 1 */
	      dcos_mphi_x[0] = 0.0;
	      dcos_mphi_x[1] = dx12t*dx12t*r12xy_inv3-r12xy_inv;
	      dcos_mphi_y[0] = 0.0;
	      dcos_mphi_y[1] = dx12t*dy12t*r12xy_inv3;

	      /* derivative sin(mPhi), m = 0, 1 */
	      dsin_mphi_x[0] = 0.0;
	      dsin_mphi_x[1] = dcos_mphi_y[1];
	      dsin_mphi_y[0] = 0.0;
	      dsin_mphi_y[1] = dy12t*dy12t*r12xy_inv3-r12xy_inv;

	      /* obtain cos(mPhi) and sin(mPhi) for m >= 2 */
	      for(k = 2; k <= l_ordr; k++)
		{
		  cos_mphi[k] = 2.0*cos_mphi[1]*cos_mphi[k-1] - cos_mphi[k-2];
		  sin_mphi[k] = 2.0*cos_mphi[1]*sin_mphi[k-1] - sin_mphi[k-2];

		  /* derivatives */
		  dcos_mphi_x[k]  = 2.0*(dcos_mphi_x[1]*cos_mphi[k-1]
					 +cos_mphi[1]*dcos_mphi_x[k-1]) - dcos_mphi_x[k-2];
		  dcos_mphi_y[k]  = 2.0*(dcos_mphi_y[1]*cos_mphi[k-1]
					 +cos_mphi[1]*dcos_mphi_y[k-1]) - dcos_mphi_y[k-2];
		  dsin_mphi_x[k]  = 2.0*(dcos_mphi_x[1]*sin_mphi[k-1]
					 +cos_mphi[1]*dsin_mphi_x[k-1]) - dsin_mphi_x[k-2];
		  dsin_mphi_y[k]  = 2.0*(dcos_mphi_y[1]*sin_mphi[k-1]
					 +cos_mphi[1]*dsin_mphi_y[k-1]) - dsin_mphi_y[k-2];
		}
	      
	      /* use GSL to calculate the correct associated Legendre polynomials */
	      for(k = 0; k <= l_ordr; k++)
		{
		  gsl_sf_legendre_Plm_deriv_array(l_ordr, k, cos_th, tmpalp, tmpalp_d);
		  alp[k] = tmpalp[l_ordr-k];
		  alp_d[k] = tmpalp_d[l_ordr-k];
		}

	      for(k = 0; k <= l_ordr; k++)
		{
		  dcos_th_x = dz12t*dx12t*r12_inv3;
		  dcos_th_y = dz12t*dy12t*r12_inv3;
		  dcos_th_z = dz12t*dz12t*r12_inv3 - r12_inv;
		  
		  /* spherical Harmionics Real, Im */
		  sh_re[k] = cos_mphi[k]*alp[k];
		  sh_im[k] = sin_mphi[k]*alp[k];
		  
		  /* ql real, imaginary */
		  qlm_re[k] += fc*sh_re[k];  
		  qlm_im[k] += fc*sh_im[k];

		  /* derivatives w.r.t particle positions */
		  AtomsNeighborList[i].QlmDerivatives[(l_ordr+1)*Ni+k].Re[0] = sh_re[k]*dfc_x + fc*(cos_mphi[k]*dcos_th_x*alp_d[k]+alp[k]*dcos_mphi_x[k]);
		  AtomsNeighborList[i].QlmDerivatives[(l_ordr+1)*Ni+k].Re[1] = sh_re[k]*dfc_y + fc*(cos_mphi[k]*dcos_th_y*alp_d[k]+alp[k]*dcos_mphi_y[k]);
		  AtomsNeighborList[i].QlmDerivatives[(l_ordr+1)*Ni+k].Re[2] = sh_re[k]*dfc_z + fc*(cos_mphi[k]*dcos_th_z*alp_d[k]);

		  AtomsNeighborList[i].QlmDerivatives[(l_ordr+1)*Ni+k].Im[0] = sh_im[k]*dfc_x + fc*(sin_mphi[k]*dcos_th_x*alp_d[k]+alp[k]*dsin_mphi_x[k]);
		  AtomsNeighborList[i].QlmDerivatives[(l_ordr+1)*Ni+k].Im[1] = sh_im[k]*dfc_y + fc*(sin_mphi[k]*dcos_th_y*alp_d[k]+alp[k]*dsin_mphi_y[k]);
		  AtomsNeighborList[i].QlmDerivatives[(l_ordr+1)*Ni+k].Im[2] = sh_im[k]*dfc_z + fc*(sin_mphi[k]*dcos_th_z*alp_d[k]);
		}
	    } /* endif: exclude pair distance > rmax_ql*/
	} /* endfor Ni loop */
    } /* endfor i loop */
  /* END : Ql calculation */

  //  comm_fac  = 1.0/((double) (Natoms*Natoms));

  ql_v = 0.0;
  for(i = 0; i <= l_ordr; i++)
    {  
      ql_v += pre_sph2[i]*(qlm_re[i]*qlm_re[i] + qlm_im[i]*qlm_im[i]);
      //     printf("qlm %d %.15g %.15g\n",i,qlm_re[i],qlm_im[i]);
    }
  //  ql_v = sqrt(ql_v*2.0*comm_fac);  // using symmetry: |Qlm|^2 = |Ql,-m|^2
  ql_v = sqrt(ql_v*2.0);  
  ql_v = ql_v/((double) Natoms);  

  Ql->x_cv = ql_v;
  //  printf("The Q%d for the Cfg file is %f\n", l_ordr, ql_v);

  /* potential/force for extended variable */
  diff = ext_var.x_ev[Ql->ind_ev] - ql_v;   /* (ExtVar6-Q6), (ExtVar4-Q4), are dimensionless */
  ext_var.f_ev[Ql->ind_ev]  = -ext_var.kappa_ev[Ql->ind_ev]*diff;   /* force in Hartree */
  Energy_ExtVar = ext_var.kappa_ev[Ql->ind_ev]*diff*diff*0.5;   /* energy in Hartree */
  //  printf("The force and energy on ExtVar%d are %f, %f\n", l_ordr, ext_var.f_ev[Ql->ind_ev], Energy_ExtVar);

  /*=======================================================================*/
  /*  iii) get forces of each pair */
  //  f_pre = comm_fac/ql_v;
  f_pre = 2.00*(-ext_var.f_ev[Ql->ind_ev])/ql_v;
  f_pre = f_pre/((double) Natoms);

  InitV9(pvten_tmp);
  
  for( i = 0; i <= (Natoms-2); i++)
    {     /* the last atom has an empty neighbor list */
      for (Ni = 0; Ni < AtomsNeighborList[i].NeighbourMax; Ni++)
        {
          j = AtomsNeighborList[i].NList[Ni];

	  DX[0] = AtomsNeighborList[i].NDist[3*Ni+0];
	  DX[1] = AtomsNeighborList[i].NDist[3*Ni+1];
	  DX[2] = AtomsNeighborList[i].NDist[3*Ni+2];

          dr = (NormV3(DX))*ANGS_BOHR;  /* dr in Bohr */

	  DX[0] = -1*DX[0]*ANGS_BOHR;
	  DX[1] = -1*DX[1]*ANGS_BOHR;
	  DX[2] = -1*DX[2]*ANGS_BOHR;

	  ForceLoc[0] = 0.00;
	  ForceLoc[1] = 0.00;
	  ForceLoc[2] = 0.00;
	  
	  if(dr < rmax_ql)  /* dr, rmax_ql in Bohr */
	    {
	      for(k = 0; k <= l_ordr; k++)
		{
		  ForceLoc[0] += pre_sph2[k]*(qlm_re[k]*AtomsNeighborList[i].QlmDerivatives[(l_ordr+1)*Ni+k].Re[0]
					      + qlm_im[k]*AtomsNeighborList[i].QlmDerivatives[(l_ordr+1)*Ni+k].Im[0]);

		  ForceLoc[1] += pre_sph2[k]*(qlm_re[k]*AtomsNeighborList[i].QlmDerivatives[(l_ordr+1)*Ni+k].Re[1]
					      + qlm_im[k]*AtomsNeighborList[i].QlmDerivatives[(l_ordr+1)*Ni+k].Im[1]);

		  ForceLoc[2] += pre_sph2[k]*(qlm_re[k]*AtomsNeighborList[i].QlmDerivatives[(l_ordr+1)*Ni+k].Re[2]
					      + qlm_im[k]*AtomsNeighborList[i].QlmDerivatives[(l_ordr+1)*Ni+k].Im[2]);
		}  /* endfor k loop */
	      /* times 2 because of symmetry: ReQlm [ReQlm]' = ReQl,-m [ReQl,-m]' */
	      
	      ForceLoc[0] = ForceLoc[0]*f_pre/((double) Natoms);
	      ForceLoc[1] = ForceLoc[1]*f_pre/((double) Natoms);
	      ForceLoc[2] = ForceLoc[2]*f_pre/((double) Natoms);

	    } /* endif */

	  clatoms.force[3*i+0] += ForceLoc[0];
	  clatoms.force[3*i+1] += ForceLoc[1];
	  clatoms.force[3*i+2] += ForceLoc[2];

	  clatoms.force[3*j+0] += -ForceLoc[0];
	  clatoms.force[3*j+1] += -ForceLoc[1];
	  clatoms.force[3*j+2] += -ForceLoc[2];
	  
	  Ql_force[3*i+0] += ForceLoc[0]/(-ext_var.f_ev[Ql->ind_ev]);
          Ql_force[3*i+1] += ForceLoc[1]/(-ext_var.f_ev[Ql->ind_ev]);
          Ql_force[3*i+2] += ForceLoc[2]/(-ext_var.f_ev[Ql->ind_ev]);

          Ql_force[3*j+0] += -ForceLoc[0]/(-ext_var.f_ev[Ql->ind_ev]);
          Ql_force[3*j+1] += -ForceLoc[1]/(-ext_var.f_ev[Ql->ind_ev]);
          Ql_force[3*j+2] += -ForceLoc[2]/(-ext_var.f_ev[Ql->ind_ev]);

          DX[0] = -1*DX[0];
          DX[1] = -1*DX[1];
          DX[2] = -1*DX[2];
	  V3V3diadicV9addV9(DX, ForceLoc, pvten_tmp);

	} /* end for Ni */
    } /* end for i */

  /* pressure tensor calculations */
  V9symmetrizeV9addV9(pvten_tmp, &(ptens.pvten_Global[0]));

  free(qlm_re);
  free(qlm_im);
  free(pre_sph2);
  free(sh_re);
  free(sh_im);
  free(cos_mphi);
  free(sin_mphi);
  free(dcos_mphi_x);
  free(dsin_mphi_x);
  free(dcos_mphi_y);
  free(dsin_mphi_y);
  free(alp);
  free(alp_d);

  return Energy_ExtVar;
}
