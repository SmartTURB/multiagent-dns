//Initialize Belief field
void init_belief(double *p){
  double fac;
  double *a;
  fac=1./(double) NS;
  for(a=p;a<p+NS;++a)
      *a = fac;
}

//Initialize problem
void init_problem(){
  int iag,temp1,temp2,temp3,temp4,temp5,iread;
  FILE *fpara;
  string fname;

  //Parameters file
  fname="../input/settings.para.n"+to_string(run_index);
  fpara=fopen(fname.c_str(),"r"); 
  
  //Read and initialize grid parameters, max time and max number of hits
  NS=NX*NY*NZ;

  iread=fscanf(fpara,"%d",&DX);
  iread=fscanf(fpara,"%d",&tMax);
  iread=fscanf(fpara,"%d",&Nhits);
  
  fprintf(stdout, "###################################################\n");
  fprintf(stdout, "# Grid NXxNYxNZ=%d x %d x %d\n",NX,NY,NZ);
  fprintf(stdout, "# Grid spacing %d\n",DX);
  fprintf(stdout, "# Max time %d\n",tMax);
  if(Nhits==2)
      fprintf(stdout, "# Using BINARY DETECTION\n");
  else if(Nhits > 2)
      fprintf(stdout, "# Using NON-BINARY DETECTION\n");
  else{
      fprintf(stdout, "# Nhits must be > 1 \n");
      exit(2);
  }
  
  //Read source and first agent positions
  iread=fscanf(fpara,"%d %d %d %d %d %d",&given_rS.x,&given_rS.y,&given_rS.z,&given_rA.x,&given_rA.y,&given_rA.z);
  
  //Read if using DNS data or not and the corresponding particle threshold
  iread=fscanf(fpara,"%d %lf %d",&tOverlap,&dns_wind,&temp5);  
  hits_from_empLk=temp5;
  fprintf(stdout, "# Using DNS data with particle threshold c_thr=%d and wind=%lf\n",c_thr,dns_wind);
  
  //Read which quantities you want to save
  iread=fscanf(fpara,"%d %d %d %d",&temp1,&temp2,&temp3,&temp4);
  save_params = temp1;
  save_trajs = temp2;
  save_end = temp3;
  save_belief_info = temp4;
  
  //Read agents
  iread=fscanf(fpara,"%d",&Nagents);  
  fprintf(stdout, "# Num of Agents= %d\n",Nagents);
  fprintf(stdout, "# All agents use emp_lk as model\n");
  
  //Allocate agents
  agent=(agent_type *) malloc(Nagents*sizeof(agent_type));
  
  //Initialize agents weigths
  different_weights=-1;
  int j=0,n;
  double weight;
  while(j < Nagents){
      //Read agents' weights
      iread=fscanf(fpara,"%d %lf",&n,&weight);
      fprintf(stdout, "# %d agents have weight w= %g\n",n,weight);
          
      if(j+n>Nagents){
          fprintf(stderr, "Number of agents does not match Nag_true=%d Nag_false=%d \n",Nagents,j+n);
          exit(0);
      }
              
      for(iag=j;iag<j+n;iag++)
          agent[iag].weight = weight;
          
      j+=n;
      different_weights++;
  }
      
  if(j!=Nagents){
      fprintf(stderr, "Number of agents does not match Nag_true=%d Nag_false=%d \n",Nagents,j);
      exit(1);
  }

  //Set possible moves
  step=(vec3d *) malloc((Nactions+1)*sizeof(vec3d));
  step[0].x= DX; step[0].y=  0; step[0].z= 0;   //x-right
  step[1].x=-DX; step[1].y=  0; step[1].z= 0;   //x-left
  step[2].x=  0; step[2].y= DX; step[2].z= 0;   //y-right
  step[3].x=  0; step[3].y=-DX; step[3].z= 0;   //y-left
  step[4].x=  0; step[4].y=  0; step[4].z= DX;  //up
  step[5].x=  0; step[5].y=  0; step[5].z= -DX; //down
  step[6].x=  0; step[6].y=  0; step[6].z= 0;   //stay

  //Allocate belief and its moments if necessary
  belief=(double *) malloc(sizeof(double)*NS);
  bmoments = (double*) malloc(sizeof(double)*4);

  //Prepare folders
  fname="../output/Parameters/";
  iread=system(("mkdir -p "+fname).c_str());
  fname="../output/Agents_info/";
  iread=system(("mkdir -p "+fname).c_str());
  fname="../output/Belief_info/";
  iread=system(("mkdir -p "+fname).c_str());
  fname="../output/Time_statistics/";
  iread=system(("mkdir -p "+fname).c_str());
}

//Initializes hits from file, agents and source positions, and return file name of dns data
string init_agents_source_belief_hits(int iep, int times, int DeltaT,int *tStart){
  
  int min_hits, iag, iterPos;
  string fname, fname_dns;
  
  //Read hits from DNS data
  
  //Needed only if by chance we start from a (pathological) time where there are no hits in the arena
   newTime:
   iterPos=0;

   int iep_randTimes = iep%(times*5);
   fname_dns = select_dns_data(iep_randTimes,times);
   *tStart = (iep%times)*tOverlap + (int)(DeltaT*drand48()); //Start from random time
   cout << "TSTART " << *tStart << endl;
    
   //Read source position from input file
   place_source_fixed();
   
   //Print source position
   cout << "Source (" << source.rS.x << ", " << source.rS.y << ", " << source.rS.z << ")" << endl;
    
   //Read first agent position
   if(given_rA.x<0){
       import_dns_snapshot(fname_dns,*tStart);
       min_hits = shift_dns_snapshot(); //Shift DNS to source position and check there is at least one hit in the arena
       if(min_hits<c_thr)
           goto newTime;   
       else{
           newPos:
           agent[0].rA = place_from_hits(); //Place first agent on a random point where there is a hit
           if(manhattan(agent[0].rA.x,agent[0].rA.y,agent[0].rA.z,given_rS.x,given_rS.y,given_rS.z)<=5){ //Condition initial agents pos
               iterPos+=1;
               if(iterPos>10)
                   goto newTime;
               goto newPos;
           }    
       }
   }
   else
	   agent[0].rA = given_rA;
	      
  //Place the remaining agents in a spiral    
  place_multi_agents();
  
  //Print agents initial positions
  cout << "Agents:" << endl; 
  for(iag=0;iag<Nagents;iag++)
      cout << agent[iag].rA.x << "   " << agent[iag].rA.y << "   " << agent[iag].rA.z << endl;

  //Initialize and renormalize belief
  init_belief(belief);
  for(iag=0;iag<Nagents;iag++)
      renormalize_belief(belief,agent[iag].rA.x,agent[iag].rA.y,agent[iag].rA.z);
  
  return fname_dns;
}



//Deallocate all arrays
void close_problem(){
  free(belief);
  free(bmoments);
  free(agent);
  free(step);
}
