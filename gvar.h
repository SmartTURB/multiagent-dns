/*  --------- Structures ---------- */

//Position on the 3D grid
struct vec3d{
  int x,y,z;
};

//Agent
struct agent_type{
  int MD0;                //Initial Manhattan distance
  double weight;          //Policy weight
  vec3d rA;               //Position
};

//Source
struct source_type{
  vec3d rS;               //Position
};

/* -------- Global Variables -------- */

int NS;                     //Lattice size for belief
int DX;                            //Grid spacing
int Nactions=6,Nagents;            //Number of actions and of agents
int Nhits,tMax;                    //Max possible hits and max possible time 
int run_index;                     //To avoid overwriting files
int different_weights;             //0=homogeneous, otherwise heterogenous policy weights
double *belief, *bmoments;         //Shared belief and its first 2 moments 
vec3d *step;                       //Moves
source_type source;                //Source
agent_type *agent;                 //Agent
vec3d given_rS,given_rA;           //Source position and agent position read from input file
bool save_params,save_trajs,save_end,save_belief_info; //Options to save info or not

//DNS simulations variables
bool hits_from_empLk;              //Option to draw hits from empirical likelihood (exact model)
double dns_wind;                   //Mean wind intensity
int c_thr,tOverlap;                //Threshold on the number of detectable particles and time shift in dns data         
double dnsHits_snapshot[NX][NY][NZ],temp3D[NX][NY][NZ]; //Original dns data and temporary array
double dns_empLk[NX][NY][NZ],dns_empLk_shifted[NX][NY][NZ]; //Empirical likelihood
