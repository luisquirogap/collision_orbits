/////////////////////////////////
// ADITIONAL CODE //
////////////////////////////////

#include<../useful_routines/useful_routines.c>

////////////////////////////////////
// GLOBAL VARIABLES //
///////////////////////////////////


////////////////////////////////////
// ROUTINES //
///////////////////////////////////


//////////////
// MAIN //
/////////////

int main(int argc, char *argv[])
{
  size_t *indexHost, *indexSat;
  int *indHost, *indSat, *idsHost, *idsSat;
  int nHostDmStars, nSatDmStars;
  int *indexHostBounds, *indexSatBounds;
  int indexmin, indexmax, i, j, k, type, nSnapshots, nBounds, counter;
  int numberParGal[2][6], nTotalHost, nTotalSat;
  double *energyHost, *energySat, *energy, v2, r;
  char *infile, infiles[200], outfile[200];
  CM cmHost, cmSat, cmHostBounds, cmSatBounds;
  
  FILE *fNumbers, *fHostOrbit, *fSatOrbit, *fHostCenter;
  
  infile = argv[1];
  nSnapshots = atoi(argv[2]);
  
  sprintf(outfile,"%s_hostOrbit.output",infile);
  fHostOrbit = fopen(outfile,"w");
  sprintf(outfile,"%s_satOrbit.output",infile);
  fSatOrbit = fopen(outfile,"w");
  sprintf(outfile,"%s_HostCenter.output",infile);
  fHostCenter = fopen(outfile,"w");
  
  fNumbers = fopen("number_particles_matrix.input","r");

  // Reading number of particles per type per galaxy
  nTotalHost = nTotalSat = 0;
  for( type=0; type<6; type++ )
    {
      returnRead = fscanf(fNumbers,"%d %d",&numberParGal[0][type],&numberParGal[1][type]);
      nTotalHost = nTotalHost + numberParGal[0][type];
      nTotalSat = nTotalSat + numberParGal[1][type];
      printf("%d %d\n",numberParGal[0][type],numberParGal[1][type]);
    }
  
  printf("totals = %d %d\n",nTotalHost,nTotalSat);
  // Reading number of bound particles to find the galatic centers
  returnRead = fscanf(fNumbers,"%d",&nBounds);
  
  fclose(fNumbers);

  // Number of collisionless particles per galaxy 
 nHostDmStars = nTotalHost - numberParGal[0][0];
 nSatDmStars = nTotalSat - numberParGal[1][0];

 // Allocating memory
  energyHost = (double *)malloc((size_t)(nHostDmStars)*sizeof(double));
  if(energyHost == NULL){
    printf("Allocation of energyHost failed\n");
    exit(0);
  }
  
  indexHost = (size_t *)malloc((size_t)(nHostDmStars)*sizeof(size_t));
  if(indexHost == NULL){
    printf("Allocation of indexHost failed\n");
    exit(0);
  }
  
  indHost = (int *)malloc((size_t)(nHostDmStars)*sizeof(int));
  if(indHost == NULL){
    printf("Allocation of indHost failed\n");
    exit(0);
  }
  
  indexHostBounds = (int *)malloc((size_t)nBounds*sizeof(int));
  if(indexHostBounds == NULL){
    printf("Allocation of indexHostBounds failed\n");
    exit(0);
  }
  
  energySat = (double *)malloc((size_t)(nSatDmStars)*sizeof(double));
  if(energySat == NULL){
    printf("Allocation of energySat failed\n");
    exit(0);
  }
  
  indexSat = (size_t *)malloc((size_t)(nSatDmStars)*sizeof(size_t));
  if(indexSat == NULL){
    printf("Allocation of indexSat failed\n");
    exit(0);
  }
  
  indSat = (int *)malloc((size_t)(nSatDmStars)*sizeof(int));
  if(indSat == NULL){
    printf("Allocation of indSat failed\n");
    exit(0);
  }
  
  indexSatBounds = (int *)malloc((size_t)nBounds*sizeof(int));
  if(indexSatBounds == NULL){
    printf("Allocation of indexSatBounds failed\n");
    exit(0);
  }

  idsHost = (int *)malloc((size_t)nBounds*sizeof(int));
  if(idsHost == NULL){
    printf("Allocation of idsHost failed\n");
    exit(0);
  }

  idsSat = (int *)malloc((size_t)nBounds*sizeof(int));
  if(idsSat == NULL){
    printf("Allocation of idsSat failed\n");
    exit(0);
  }

  // Reading initial conditions
  sprintf(infiles,"%s_%.3d",infile,0);
  read_gadget1(infiles);

  
  energy = (double *)malloc((size_t)(N_part_total)*sizeof(double));
  if(energy == NULL){
    printf("Allocation of energy failed\n");
    exit(0);
  }

  // Particles index for each galaxy
  counter = 0;
  for( type=1; type<6; type++ )
    {      
      if( N_part[type]>0 )
	{
	  indexmin = 0;
	  for( i=0; i<type; i++ )
	    indexmin = indexmin + N_part[i];
	  indexmax =  indexmin + numberParGal[0][type];
	  
	  for( i=indexmin; i<indexmax; i++)
	    indHost[counter++] = i;
	}
    }
  
  counter = 0;
  for( type=1; type<6; type++ )
    {      
      if( N_part[type]>0 )
	{
	  printf("Computing energy for type %d\n",type);
	  indexmin = 0;
	  for( i=0; i<type; i++ )
	    indexmin = indexmin + N_part[i];
	  indexmin =  indexmin + numberParGal[0][type];
	  
	  indexmax = indexmin + numberParGal[1][type];
	  
	  for( i=indexmin; i<indexmax; i++)
	    indSat[counter++] = i;
	  
	}
    }

  // Initial center of mass for each galaxy
  centerMass(indHost, nHostDmStars, 0, 0, &cmHost);
  centerMass(indSat, nSatDmStars, 0, 0, &cmSat);
  
  printf("Host ini %d %.8lf %.8lf %.8lf %.8lf %.8lf %.8lf\n",0,
	 cmHost.cm[X],cmHost.cm[Y],cmHost.cm[Z],
	 cmHost.vcm[X],cmHost.vcm[Y],cmHost.vcm[Z] );
  printf("Sat ini %d %.8lf %.8lf %.8lf %.8lf %.8lf %.8lf\n",0,
	 cmSat.cm[X],cmSat.cm[Y],cmSat.cm[Z],
	 cmSat.vcm[X],cmSat.vcm[Y],cmSat.vcm[Z] );
  
  totalTranslationMinus(&cmHost, N_part_total);

  // Centers of mass in center of host
  centerMass(indHost, nHostDmStars, 0, 0, &cmHost);
  centerMass(indSat, nSatDmStars, 0, 0, &cmSat);
  
  printf("Host in Host %d %.8lf %.8lf %.8lf %.8lf %.8lf %.8lf\n",0,
	 cmHost.cm[X],cmHost.cm[Y],cmHost.cm[Z],
	 cmHost.vcm[X],cmHost.vcm[Y],cmHost.vcm[Z] );
  printf("Sat in Host %d %.8lf %.8lf %.8lf %.8lf %.8lf %.8lf\n",0,
	 cmSat.cm[X],cmSat.cm[Y],cmSat.cm[Z],
	 cmSat.vcm[X],cmSat.vcm[Y],cmSat.vcm[Z] );

  // Energy por all collisionless particles 
  FILE *pf1 = fopen("energia.dat","w");
  //for( i=N_part[0]; i<N_part_total; i++)
  for( i=N_part[0]; i<N_part_total; i++)
    {
      v2 =  particles[i].vel[X]*particles[i].vel[X]
	+ particles[i].vel[Y]*particles[i].vel[Y]
	+ particles[i].vel[Z]*particles[i].vel[Z];
      
      energy[i] = particles[i].mass*particles[i].pot + 0.5*particles[i].mass*v2;
      // energy[i-N_part[0]] = particles[i].mass*particles[i].pot + 0.5*particles[i].mass*v2;
      
      fprintf(pf1,"%d %.8lf %.8lf %.8lf\n",i,
	      particles[i].pos[X],particles[i].pos[Y],
	      energy[i]);
      if(energy[i]==0.0)
	printf("%lf\n",energy[i]);
    }
  fclose(pf1);

  // Energy for host galaxy in its center
  FILE *pf2 = fopen("energia_host.dat","w");
  for( i=0; i<nHostDmStars; i++)
    {
      fprintf(pf2,"%d %d %.8lf %.8lf %.8lf\n",i,indHost[i],
	      particles[indHost[i]].pos[X],particles[indHost[i]].pos[Y],
	      energy[indHost[i]]);
      //if(energy[indHost[i]]==0.0)
      //printf("%d %d\n",i,indHost[i]);
    }
  fclose(pf2);

  // Energy for satellite galaxy in host center
  FILE *pf3 = fopen("energia_sat.dat","w");
  for( i=0; i<nSatDmStars; i++)
    fprintf(pf3,"%d %d %.8lf %.8lf %.8lf\n",i,indSat[i],
	    particles[indSat[i]].pos[X],particles[indSat[i]].pos[Y],
	    energy[indSat[i]]); 
  fclose(pf3);

  // More bound particles for host in host center
  for( i=0; i<nHostDmStars; i++)
    energyHost[i] = energy[indHost[i]];
  
  gsl_sort_index(indexHost,energyHost,1,nHostDmStars);
  
  for( i=0; i<nBounds; i++ )
    indexHostBounds[i] = indHost[indexHost[i]];
  
  centerMass(indexHostBounds, nBounds, 0, 0, &cmHostBounds);
  printf("Host cm bounds in host %d %.8lf %.8lf %.8lf %.8lf %.8lf %.8lf\n",0,
	 cmHostBounds.cm[X],cmHostBounds.cm[Y],cmHostBounds.cm[Z],
	 cmHostBounds.vcm[X],cmHostBounds.vcm[Y],cmHostBounds.vcm[Z] );
  
  fprintf(fHostCenter,"%d %.8e %.8e %.8e %.8e %.8e %.8e\n",0,
	  cmHostBounds.cm[X],cmHostBounds.cm[Y],cmHostBounds.cm[Z],
	  cmHostBounds.vcm[X],cmHostBounds.vcm[Y],cmHostBounds.vcm[Z] );

  // More bounds particles for satellite in satellite center
  totalTranslationMinus(&cmSat, N_part_total);
  
  for( i=0; i<nSatDmStars; i++)
    {
      v2 =  particles[indSat[i]].vel[X]*particles[indSat[i]].vel[X]
	+ particles[indSat[i]].vel[Y]*particles[indSat[i]].vel[Y]
	+ particles[indSat[i]].vel[Z]*particles[indSat[i]].vel[Z];
      
      energySat[i] = particles[indSat[i]].mass*particles[indSat[i]].pot + 0.5*particles[indSat[i]].mass*v2;      
    }
  
  gsl_sort_index(indexSat,energySat,1,nSatDmStars);

  // More bound particles for satellite within 40 kpc of its center
  counter = 0;
  for( i=0; i<nSatDmStars; i++ )
    {
      r = sqrt( particles[indSat[indexSat[i]]].pos[X]*particles[indSat[indexSat[i]]].pos[X]
		+ particles[indSat[indexSat[i]]].pos[Y]*particles[indSat[indexSat[i]]].pos[Y]
		+ particles[indSat[indexSat[i]]].pos[Z]*particles[indSat[indexSat[i]]].pos[Z] );
      
      if( r < 5.0 )
	{
	  indexSatBounds[counter] = indSat[indexSat[i]];
	  counter++;
	  if(counter==nBounds)
	    break;
	}
    }

  // Translation to center of host
  centerMass(indHost, nHostDmStars, 0, 0, &cmHost);
  totalTranslationMinus(&cmHost, N_part_total);
  
  centerMass(indexSatBounds, nBounds, 0, 0, &cmSatBounds);
  printf("Sat cm Bounds in host %d %.8lf %.8lf %.8lf %.8lf %.8lf %.8lf\n",0,
	 cmSatBounds.cm[X],cmSatBounds.cm[Y],cmSatBounds.cm[Z],
	 cmSatBounds.vcm[X],cmSatBounds.vcm[Y],cmSatBounds.vcm[Z] );
  
    
  FILE *pf4 =  fopen("energia_ligadas_host.dat","w");
  FILE *pf5 =  fopen("energia_ligadas_sat.dat","w");
  
  // Index for nBounds particles more bound
  for( i=0; i<nBounds; i++ ){
    fprintf(pf4,"%d %.8lf %.8lf %.8lf %.8lf\n",
	    i,
	    particles[indexHostBounds[i]].pos[X],particles[indexHostBounds[i]].pos[Y],
	    energy[indexHostBounds[i]],energyHost[indexHost[i]]);
  }
  for( i=0; i<nBounds; i++ ){
    fprintf(pf5,"%d %.8lf %.8lf %.8lf %.8lf\n",
	    i,
	    particles[indexSatBounds[i]].pos[X],particles[indexSatBounds[i]].pos[Y],
	    energy[indexSatBounds[i]],energySat[indexSat[i]]);
  }
  
  fclose(pf4);
  fclose(pf5);

  free(energy);
  free(energyHost);
  free(energySat);
  free(indexHost);
  free(indexSat);
  free(indHost);
  free(indSat);

  // Translation to center of mass of more bonds for host 
  totalTranslationMinus(&cmHostBounds, N_part_total); 
  
  // Center for Host  and satellite in center of host
  centerMass(indexHostBounds, nBounds, 0, 0, &cmHost);
  centerMass(indexSatBounds, nBounds, 0, 0, &cmSat);
  
  printf("Host in host bounds %d %.8lf %.8lf %.8lf %.8lf %.8lf %.8lf\n",0,
	 cmHost.cm[X],cmHost.cm[Y],cmHost.cm[Z],
	 cmHost.vcm[X],cmHost.vcm[Y],cmHost.vcm[Z] );
  printf("Sat in host bounds %d %.8lf %.8lf %.8lf %.8lf %.8lf %.8lf\n",0,
	 cmSat.cm[X],cmSat.cm[Y],cmSat.cm[Z],
	 cmSat.vcm[X],cmSat.vcm[Y],cmSat.vcm[Z] );
   
  fprintf(fHostOrbit,"%d %.8lf %.8lf %.8lf %.8lf %.8lf %.8lf\n",0,
	  cmHost.cm[X], cmHost.cm[Y], cmHost.cm[Z],
	  cmHost.vcm[X], cmHost.vcm[Y], cmHost.vcm[Z]);
  
  fprintf(fSatOrbit,"%d %.8lf %.8lf %.8lf %.8lf %.8lf %.8lf\n",0,
	  cmSat.cm[X], cmSat.cm[Y], cmSat.cm[Z],
	  cmSat.vcm[X], cmSat.vcm[Y], cmSat.vcm[Z]);

  FILE *fLigadas;
  sprintf(infiles,"ligadas_sat_%.3d",0);
  fLigadas = fopen(infiles,"w");
  for( j=0; j<nBounds; j++)
    {
      fprintf(fLigadas,"%d %u %lf %lf %lf\n",j,
	      particles[indexSatBounds[j]].id,
	      particles[indexSatBounds[j]].pos[X],
	      particles[indexSatBounds[j]].pos[Y],
	      particles[indexSatBounds[j]].pos[Z]);
    }
  fclose(fLigadas);

  for( i=0; i<nBounds; i++ )
    {
      idsHost[i] = particles[indexHostBounds[i]].id;
      idsSat[i] = particles[indexSatBounds[i]].id;
    }

  // center of more bound particles in all simulation
  for( i=1; i<=nSnapshots; i++ )
    {
      sprintf(infiles,"%s_%.3d",infile,i);
      
      read_gadget1(infiles);
      indexmin = N_part[0];
      indexmax = N_part[0]+N_part[1]+N_part[2]+N_part[3];
      
      for( j=indexmin; j<indexmax; j++)
	{
	  for( k=0; k<nBounds; k++ )
	    {
	      if( particles[j].id == idsHost[k] )
		indexHostBounds[k] = j;
	      if( particles[j].id == idsSat[k] )
		indexSatBounds[k] = j;
	    }
	}
	
      // Center for Host
      centerMass(indexHostBounds, nBounds, 0, 0, &cmHostBounds);
      fprintf(fHostCenter,"%d %.8e %.8e %.8e %.8e %.8e %.8e\n",i,
	      cmHostBounds.cm[X],cmHostBounds.cm[Y],cmHostBounds.cm[Z],
	      cmHostBounds.vcm[X],cmHostBounds.vcm[Y],cmHostBounds.vcm[Z] );
      
      //Traslating galaxy to center of mass of the host
      totalTranslationMinus(&cmHostBounds, N_part_total);
      // Center for Host
      sprintf(infiles,"ligadas_sat_%.3d",i);
      fLigadas = fopen(infiles,"w");
      centerMass(indexHostBounds, nBounds, 0, 0, &cmHostBounds);
      centerMass(indexSatBounds, nBounds, 0, 0, &cmSatBounds);     
      for( j=0; j<nBounds; j++)
	{
	  fprintf(fLigadas,"%d %u %lf %lf %lf\n",j,
		  particles[indexSatBounds[j]].id,
		  particles[indexSatBounds[j]].pos[X],
		  particles[indexSatBounds[j]].pos[Y],
		  particles[indexSatBounds[j]].pos[Z]);
	}
      
      fprintf(fHostOrbit,"%d %.8lf %.8lf %.8lf %.8lf %.8lf %.8lf\n",i,
	      cmHostBounds.cm[X], cmHostBounds.cm[Y], cmHostBounds.cm[Z],
	      cmHostBounds.vcm[X], cmHostBounds.vcm[Y], cmHostBounds.vcm[Z]);
      
      fprintf(fSatOrbit,"%d %.8lf %.8lf %.8lf %.8lf %.8lf %.8lf\n",i,
	      cmSatBounds.cm[X], cmSatBounds.cm[Y], cmSatBounds.cm[Z],
	      cmSatBounds.vcm[X], cmSatBounds.vcm[Y], cmSatBounds.vcm[Z]);
      
      free(particles);
      free(gaspro);
      fclose(fLigadas);
    }
  
  fclose(fHostOrbit);
  fclose(fSatOrbit);
  fclose(fHostCenter);
  free(indexHostBounds);
  free(indexSatBounds);
  free(idsHost);
  free(idsSat);
  
  return 0;
}
