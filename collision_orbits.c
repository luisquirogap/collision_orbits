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
  int *indHost, nSnapshot, nTimes, snapToPrint;
  int indexmin, indexmax, i, type, nBounds, counter;
  int numberParGal[2][6], nTotalHost, nTotalSat;
  char *baseName, *infile, infiles[200], outfile[200];
  CM cmHost;
  
  FILE *fNumbers;
  
  baseName = argv[1];
  nSnapshot = atoi(argv[2]);
  
  fNumbers = fopen("number_particles_matrix.input","r");
  
  nTotalHost = nTotalSat = 0;
  for( type=0; type<6; type++ )
    {
      returnRead = fscanf(fNumbers,"%d %d",&numberParGal[0][type],&numberParGal[1][type]);
      nTotalHost = nTotalHost + numberParGal[0][type];
      nTotalSat = nTotalSat + numberParGal[1][type];
      printf("%d %d\n",numberParGal[0][type],numberParGal[1][type]);
    }
  
  printf("totals = %d %d\n",nTotalHost,nTotalSat);
  returnRead = fscanf(fNumbers,"%d",&nBounds);
  
  fclose(fNumbers);
  
  indHost = (int *)malloc((size_t)(nTotalHost-numberParGal[0][0])*sizeof(int));
  if(indHost == NULL){
    printf("Allocation of indHost failed\n");
    exit(0);
  }
 
  sprintf(infiles,"%s_%.3d",baseName,nSnapshot);
  read_gadget1(infiles);
  
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
  
  FILE *fHostCenter;
  infile =  argv[3];
  nTimes = counterLines(infile);
  fHostCenter = fopen(infile,"r");
  
  for( i=0; i<nTimes; i++ )
    {
      returnRead = fscanf(fHostCenter,"%d %lf %lf %lf %lf %lf %lf",
			  &snapToPrint,
			  &cmHost.cm[X],&cmHost.cm[Y],&cmHost.cm[Z],
			  &cmHost.vcm[X],&cmHost.vcm[Y],&cmHost.vcm[Z]);
      if( snapToPrint == nSnapshot )
	break;
    }
  
  printf("Host center %d %.8lf %.8lf %.8lf %.8lf %.8lf %.8lf\n",snapToPrint,
	 cmHost.cm[X],cmHost.cm[Y],cmHost.cm[Z],
	 cmHost.vcm[X],cmHost.vcm[Y],cmHost.vcm[Z] );
  
  totalTranslationMinus(&cmHost, N_part_total);
  
  
  //****************************************************
  // Writing data for all particles in ASCII 
  //**************************************************** 
  
  N_min = 0;
  
  for( type=0; type<6; type++)
    {
      
      N_max = N_min + N_part[type];	 
      
      if( N_part[type] > 0 )
	{
	  sprintf(outfile,"%s_centered.%d",infiles,type);
	  FILE *outfiles;
	  outfiles = fopen(outfile,"w");
	  if(outfiles==NULL) printf("No se pudo abrir %s\n",outfile);  
	  
	  for(i=N_min;i<N_max;i++)
	    {
	      
	      
#ifdef LONGIDS
	      fprintf(outfiles,"%lu",particles[i].id);
#else
	      fprintf(outfiles,"%u",particles[i].id);
#endif
	      
	      fprintf(outfiles," %f %f %f %f %f %f %f",
		      particles[i].pos[0],particles[i].pos[1],particles[i].pos[2],
		      particles[i].vel[0],particles[i].vel[1],particles[i].vel[2],
		      particles[i].mass);
	      
	      if( type == 0 &&  N_part[0] >0 )
		{
		  fprintf(outfiles," %f %f",gaspro[i].U,gaspro[i].rho);
#ifdef COOLING
		  fprintf(outfiles," %f %f",gaspro[i].Ne,gaspro[i].Nh);
#endif
		  fprintf(outfiles," %f",gaspro[i].h);
#ifdef SFR
		  fprintf(outfiles," %f",gaspro[i].sfr);
#endif
		}
	      
#ifdef OUTPUTPOTENTIAL
	      fprintf(outfiles," %f",particles[i].pot);
#endif
	      
#ifdef OUTPUTACCELERATION
	      fprintf(outfiles," %f",particles[i].acce);
#endif
	      
	      
#ifdef OUTPUTCHANGEOFENTROP
	      if( type == 0 && N_part[0] >0 )
		{
		  fprintf(outfiles," %f",gaspro[i].ecr);
		}
#endif 
	      
	      
#ifdef OUTPUTTIMESTEP
	      fprintf(outfiles," %f",particles[i].timestep);
#endif
	      
	      fprintf(outfiles,"\n");	
	      
	    }
	  fclose(outfiles);
	}

      N_min = N_max;
      
    }
  
  printf("Done.\n\n");
  
  
  return 0;
}
