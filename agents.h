//Initialize the remaining agents positions
void place_multi_agents(){
  int iag=1;
  int layer=1,leg=0;

  //Place agents in an outward spiral (orthogonal to the mean wind) starting from first agent
  while(iag < Nagents){
  
      agent[iag].rA.x=agent[0].rA.x;
      
      if(leg == 0){
            agent[iag].rA.z = agent[iag-1].rA.z ;
            agent[iag].rA.y = agent[iag-1].rA.y + (int) DX;
            if(agent[iag].rA.y-agent[0].rA.y == layer)
                leg += 1;
            iag++;
            continue;
      }
      if(leg == 1){
            agent[iag].rA.z = agent[iag-1].rA.z- (int) DX;
            agent[iag].rA.y = agent[iag-1].rA.y ;
            if(agent[iag].rA.z-agent[0].rA.z == -layer)
                leg += 1;
            iag++;
            continue;
      }
      if(leg == 2){
            agent[iag].rA.z = agent[iag-1].rA.z ;
            agent[iag].rA.y = agent[iag-1].rA.y - (int) DX;
            if(agent[iag].rA.y-agent[0].rA.y == -layer)
                leg += 1;
            iag++;
            continue;
      }
      if(leg == 3){
            agent[iag].rA.z = agent[iag-1].rA.z + (int) DX;
            agent[iag].rA.y = agent[iag-1].rA.y ;
            if(agent[iag].rA.z-agent[0].rA.z == layer){
                leg = 0;
                layer += 1;
            }
            iag++;
            continue;
      }
  }
}

//Weighted Space Aware Infotaxis
#define TH 0.9999
int WSAI(int iag, double *bb){
 int i,nh,action,a;
 int overlap[Nactions], outside[Nactions];
 vec3d rnew;
 double *pp;
 double cum;
 double JJ[Nactions],minJ,probfound,probnh;
 double ES,SS,EMD,MD; 
 int memsize = NS*sizeof(double);

 pp=(double *) malloc(memsize);
    
 for(a=0;a<Nactions;a++){

   rnew.x=agent[iag].rA.x+step[a].x;
   rnew.y=agent[iag].rA.y+step[a].y;
   rnew.z=agent[iag].rA.z+step[a].z;
   
   overlap[a]=0;
   outside[a]=0;

   //Check if within arena
   if(rnew.x > NX-1 || rnew.y > NY-1 || rnew.z > NZ-1 || rnew.x < 0 || rnew.y < 0 || rnew.z < 0){
       JJ[a] = DBL_MAX;
       outside[a]=1;
       continue;
   }
   
   //Check if action makes agents overlap
   for(i=0;i<iag;i++){
       if(rnew.x==agent[i].rA.x && rnew.y==agent[i].rA.y && rnew.z==agent[i].rA.z){
           JJ[a] = DBL_MAX;
           overlap[a]=1;
	   //cout << "Agent " << iag << ", Action not to take " << a << endl;
           continue;
       }
   }
   
   //Compute probability of finding the source
   probfound = bb[(rnew.x*NY + rnew.y)*NZ + rnew.z];
   
   //If prob of finding the source is 1 and you do not overlap, then take that action
   if(fabs(probfound-1.)<DBL_EPSILON && overlap[a]==0){
       return a;
   }
       
   ES=EMD=0.0; 
   cum=0.0;
   nh=0;
   while(cum<TH && nh<Nhits-1){      
       memcpy(pp,bb,memsize);
       renormalize_belief(pp,rnew.x,rnew.y,rnew.z);
       probnh=BayesUpdate(rnew.x,rnew.y,rnew.z,pp,nh, &SS, &MD);       
       cum+=probnh;
       ES+=(1.-probfound)*probnh*SS;
       EMD+=(1.-probfound)*probnh*MD;
       nh++;
   }
   probnh=1.-cum;
   memcpy(pp,bb,memsize);
   renormalize_belief(pp,rnew.x,rnew.y,rnew.z);
   BayesUpdate(rnew.x,rnew.y,rnew.z,pp,nh, &SS, &MD);  
   ES+=(1.-probfound)*probnh*SS; 
   EMD+=(1.-probfound)*probnh*MD; 
   
   JJ[a] = (1.0-agent[iag].weight)*(pow(2,ES)-1.)*0.5 + agent[iag].weight*EMD;
   
   if(isnan(JJ[a]))
     JJ[a]=DBL_MAX;
 }

 //Select action according to policy
 minJ=JJ[0];
 action=0;
 for(a=1;a<Nactions;a++){
   if(JJ[a]<(minJ-DBL_EPSILON)){
     minJ=JJ[a];
     action=a;
   }
 }
 
 /*vectint act_array;
 int randomIndex;
 
 for(a=0;a<Nactions;a++){
   if(fabs(JJ[a]-JJ[action])<DBL_EPSILON)
     act_array.push_back(a);
 }
 
 randomIndex = (int) (drand48()*act_array.size());
 action = act_array[randomIndex];
 act_array.clear();
 
 */
 
 //If all actions make it overlap or exit the grid, then agent must stay in place
 if(overlap[action] || outside[action]){
     action=Nactions;
     if(minJ>DBL_MAX-1.)
         cout << "BELIEF IS ZERO FOR ALL ACTIONS: FAILED!" << endl;
 }
 
 free(pp);
 return action;
}
