/* ********************************* *\
       ODOR SEARCH CODE 3D
     Arbitrary number of Agents
\* ********************************* */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <ctime>
#include <vector>
#include <cfloat>
#include <set>
#include <algorithm>
#include <random>

using namespace std;

#define vect vector<double>
#define vectint vector<int>
#define NX 129
#define NY 99
#define NZ 99
#define NTIME 2670

#include "../headers/gvar.h"
#include "../headers/utility.h"
#include "../headers/io.h"
#include "../headers/source.h"
#include "../headers/dns_data.h"
#include "../headers/bayes.h"
#include "../headers/agents.h"
#include "../headers/initclosefuns.h"

int main(int argc,char *argv[]){
  int t,iep,iag,iagFound; //Useful indexes
  int nh,action,found,sumLost; //Agents' stuff
  int NEpisodes,NEpisodes_in; //Number of episodes
  int seed; //Rand seed
  int times,tStart=0,tEnd,DeltaT=100; //Useful stuff for DNS
  double S,MD; //Entropy and Manhattan distance of the belief
  vectint agent_info, obs;
  vect end_info, belief_info;
  string fname,fname_dns; //File names
  clock_t start_time, end_time; //Profiling
 
  //Start clock for profiling
  start_time = clock();
  
  //Read number of episodes and run index from command line
  if(argc!=4){
      fprintf(stderr,"usage %s <NEpisodes>\n",argv[0]);
      exit(1);
  }
  NEpisodes_in=atoi(argv[1]);
  run_index=atoi(argv[2]);
  c_thr=atoi(argv[3]);

  //Seed rng
  seed = run_index+1720790456;      
  srand48(seed);
  
  //Initialize grid, source, agents, moves and likelihood
  init_problem();
  
  //Print rng seed
  cout << "Rand seed: " << seed << endl;

  //# Episodes when using overlapping cuts of the DNS time series from the 5 sources
  times = max(1,(NTIME-tMax-DeltaT)/tOverlap+1); //Set tOverlap = tMax if you want non-overlapping cuts of the DNS
  NEpisodes = times*5*NEpisodes_in; //SET NEpisodes_in = 25 if you want 10^3 total episodes
  tEnd = min(tMax,NTIME);
  
  //Initialize DNS empirical likelihood
  fname_dns = select_dns_empLk();
  import_dns_empLk(fname_dns);
  shift_emp_lk();
    
  /*string output_file = "output_data.dat";
  export_dns_empLk(output_file);*/
      
  //Start loop over episodes
  for(iep=0;iep<NEpisodes;iep++){
      
      //Initializes hits from file, the agents and source positions, the belief(s), and returns dns data file 
      fname_dns = init_agents_source_belief_hits(iep,times,DeltaT,&tStart);
             
      //Read first hit for each agent to start the search
      for(iag=0;iag<Nagents;iag++){
          if(hits_from_empLk) //Generate hit from empLk 
              nh = rndUniform(agent[iag].rA.x,agent[iag].rA.y,agent[iag].rA.z);
          else //Take hit directly from dns
              nh = (dnsHits_snapshot[agent[iag].rA.x][agent[iag].rA.y][agent[iag].rA.z] >= c_thr) ? 1 : 0;
          
                  
          BayesUpdate(agent[iag].rA.x, agent[iag].rA.y, agent[iag].rA.z, belief, nh, &S, &MD);

          obs.push_back(nh);
      }
      
      //Compute initial Manhattan distance
      compute_iMD();
   
      //Reset initial values before time loop
      found=sumLost=0;
      iagFound=0;
      t=0;

      //Save source and agents parameters on a file
      if(save_params){
          fname="../output/Parameters/" + get_name_params_specific(NEpisodes,run_index) + ".dat";
          Save_params(fname,seed,iep,(iep==0?0:1));
      }
      
      //Store agent information at the beginning
      if(save_trajs){
          for(iag=0;iag<Nagents;iag++)
              Store_agent_info(agent_info,t,iag,obs[iag],iep);
      }
      
      //Store some stats about current belief
      if(save_belief_info){
          belief_moments(belief,bmoments,&S);
          Store_belief_info(belief_info,bmoments,S,t,iep);
      }
   
      //Save initial belief
      #ifdef PRINTBELIEF
          fname="../output/Belief_info/belief_" + get_name_params_specific_time(NEpisodes,run_index,iep,t) +".bin" ;
          Save_binary_double(fname,belief,NS);
      #endif
            
      //Start time loop
      while(found==0 && t<tEnd){
          
          //Import dns data
          import_dns_snapshot(fname_dns,t+tStart);
          shift_dns_snapshot();
          /*if(t==0 && iep==0){
                fname="dns_prova.bin";
                export_dns_snapshot(fname);
          }*/

          //Let agents take actions based on common belief and check if found or lost
          for(iag=0;iag<Nagents;iag++){
              if(!found){
              
                 //Determine action according to policy and common belief
                 action=WSAI(iag,belief);
                 
                 //Move agent
                 agent[iag].rA.x += step[action].x;
                 agent[iag].rA.x = min(max(agent[iag].rA.x,0),NX-1); //Ensures agents don't exit the boundaries
                 agent[iag].rA.y += step[action].y;     
                 agent[iag].rA.y = min(max(agent[iag].rA.y,0),NX-1); //Ensures agents don't exit the boundaries
                 agent[iag].rA.z += step[action].z;
                 agent[iag].rA.z = min(max(agent[iag].rA.z,0),NZ-1); //Ensures agents don't exit the boundaries
                  
                  //Check if agent found the source
                  if((pow(agent[iag].rA.x-source.rS.x,2)+pow(agent[iag].rA.y-source.rS.y,2)+pow(agent[iag].rA.z-source.rS.z,2))<0.9){
                      found=1;
                      iagFound=iag;
                      t++;
                   
                      //Store agent information at the end
                      if(save_trajs)
                          Store_agent_info(agent_info,t,iag,obs[iag],iep);
                  }
              }   
          }
         
          //Empty measurements vector
          obs.clear();
          
          //Put zeros where agents moved and renormalize belief
          for(iag=0;iag<Nagents;iag++)
             renormalize_belief(belief,agent[iag].rA.x,agent[iag].rA.y,agent[iag].rA.z);
     
          //Perform measurement and update common belief
          if(!found){
              t++; 
                  
              for(iag=0;iag<Nagents;iag++){  
              
                  //Store agent information at each timestep
                  if(save_trajs)
                      Store_agent_info(agent_info,t,iag,obs[iag],iep); 
                      
                  if(hits_from_empLk) //Generate hit from empLk 
                      nh = rndUniform(agent[iag].rA.x,agent[iag].rA.y,agent[iag].rA.z);
                  else //Take hit directly from dns
                      nh = (dnsHits_snapshot[agent[iag].rA.x][agent[iag].rA.y][agent[iag].rA.z] >= c_thr) ? 1 : 0;
                  
                  BayesUpdate(agent[iag].rA.x, agent[iag].rA.y, agent[iag].rA.z, belief, nh, &S, &MD);
                  
                  //Store measurements in a temp vector
                  obs.push_back(nh);
              }
          }
     
          //Store some stats about current belief
          if(save_belief_info){
              belief_moments(belief,bmoments,&S);
              Store_belief_info(belief_info,bmoments,S,t,iep);
          }
     
          //Save belief at time t
          #ifdef PRINTBELIEF
              fname="../output/Belief_info/belief_" + get_name_params_specific_time(NEpisodes,run_index,iep,t) +".bin";
              Save_binary_double(fname,belief,NS);
          #endif
      }
      //End loop over time
   
      //Consider all agents lost if episdode reaches max time allowed
      if(t==tEnd)
          sumLost=Nagents;
       
      //Output messages at the end of each episode
      if(found)
          cout << "In Ep = " << iep << " Agent = " << iagFound << " found the source in position (" << agent[iagFound].rA.x << ", " << agent[iagFound].rA.y << ", " << agent[iagFound].rA.z << ") at t = " << t << " while starting at dist = " << agent[iagFound].MD0 << endl;
      else
          cout << "In Ep = " << iep << " All agents got lost at t = " << t << endl;
   
      //Store and save useful info at the end of each episode
      if(save_end){
          end_info.push_back(iep);
          end_info.push_back(t);
          end_info.push_back(agent[iagFound].MD0);
          end_info.push_back(Nagents-sumLost);
          end_info.push_back(iagFound);
          end_info.push_back(agent[iagFound].weight);
          for(iag=0;iag<Nagents;iag++)
              end_info.push_back(manhattan(agent[iag].rA.x,agent[iag].rA.y,agent[iag].rA.z,source.rS.x,source.rS.y,source.rS.z));

          fname="../output/Time_statistics/" + get_name_params_specific(NEpisodes,run_index) + ".dat";
          Save_vector(fname,end_info,6+Nagents,(iep==0?0:1));
          end_info.clear();
      }
      
      //Save trajectories
      if(save_trajs){
          fname="../output/Agents_info/" + get_name_params_specific(NEpisodes,run_index) + ".bin";
          Save_vectorint_bin(fname,agent_info,8,(iep==0?0:1));
          agent_info.clear();
      }
   
      //Save belief moments and entropy at the end of each episode
      if(save_belief_info){
          fname="../output/Belief_info/" + get_name_params_specific(NEpisodes,run_index) + ".bin";
          Save_vector_bin(fname,belief_info,9,(iep==0?0:1));
          belief_info.clear();
      }
   
  }
  //End loop over episodes
 
  //Deallocate arrays
  close_problem();
 
  //Stop clock for profiling
  end_time = clock();
  print_time(start_time,end_time);
 
  exit(0);  
}
