/* ---------------------------- Random routines ---------------------------- */

//Draw random variable from uniform distribution and compare with dns data
int rndUniform(int agent_x, int agent_y, int agent_z){
  double u=drand48();
  int hit;
  
  if(u>=dns_empLk_shifted[agent_x][agent_y][agent_z])
      hit = 0;
  else
      hit = 1;
  
  return hit;
}

/* ------------------------------------------------------------------------- */

/* ------------- Distance and belief moments routines ---------------------- */

int manhattan(int ix1,int iy1,int iz1,int ix2,int iy2,int iz2){
  return abs(ix1-ix2)+abs(iy1-iy2)+abs(iz1-iz2);
}

void compute_iMD(){
	int iag;
	for(iag=0;iag<Nagents;iag++)
        	agent[iag].MD0 = manhattan(agent[iag].rA.x,agent[iag].rA.y,agent[iag].rA.z,source.rS.x,source.rS.y,source.rS.z);
}

double dist(double x1,double y1,double z1,double x2,double y2,double z2){
  double dx,dy,dz;
  dx=(x1-x2)*DX;
  dy=(y1-y2)*DX;
  dz=(z1-z2)*DX;
  return sqrt(dx*dx+dy*dy+dz*dz);
}  

void belief_moments(double* p, double* moments, double* S){
	int i,j,l,k,n;

	*S=0.;
	for(n=0;n<6;n++)
        	moments[n] = 0.;

	for(i=0;i<NX;i++){
		for(j=0;j<NY;j++){
    		for(l=0;l<NZ;l++){
    			k = (NY*i+j)*NZ+l;

    			if(p[k]>0.)
    				*S -= p[k]*log(p[k]);
				
        		//*MD += (abs(i-agent[0].rA.x)+abs(j-agent[0].rA.y))*p[k];
			
    			for(n=0;n<3;n++){
    				moments[2*n]+=p[k]*pow(i,n+1);
    				moments[2*n+1]+=p[k]*pow(j,n+1);
    			}
    		}
    	}
    }
    
	*S/=log(2.);

	return;
}

/* ------------------------------------------------------------------------- */

/* ---------------------- Routines to store data --------------------------- */

void Store_agent_info(vectint& agent_info, const int t, const int iag, const int nh, const int iep){
  agent_info.push_back(iep);
  agent_info.push_back(t);
  agent_info.push_back(iag);
  agent_info.push_back(agent[iag].rA.x);
  agent_info.push_back(agent[iag].rA.y);
  agent_info.push_back(agent[iag].rA.z);
  agent_info.push_back(nh);
  agent_info.push_back(manhattan(agent[iag].rA.x,agent[iag].rA.y,agent[iag].rA.z,source.rS.x,source.rS.y,source.rS.z));
}

void Store_belief_info(vect& belief_info, double* m, const double S, const int t, const int iep){
  belief_info.push_back(iep);
  belief_info.push_back(t);
  //mean X Y Z
  belief_info.push_back(m[0]);
  belief_info.push_back(m[1]);
  belief_info.push_back(m[2]);
  //second moment X Y Z
  belief_info.push_back(m[3]);
  belief_info.push_back(m[4]);
  belief_info.push_back(m[5]);
  //entropy
  belief_info.push_back(S);
}

/* ------------------------------------------------------------------------- */

/* ------------------------ Routines for output ---------------------------- */

static inline string int_to_string(const int a){
  ostringstream str;
  str << a;
  return str.str();
}

static inline string double_to_string(const double a, const int precision){
  ostringstream str;
  str << fixed << setprecision(precision) << a;
  return str.str();
}

string get_name_params_specific(const int Nepis, const int n){
    string a;

    a = "dnsData_cThr" + int_to_string(c_thr) + "_windDns" + double_to_string(dns_wind,1) + "_empLk";

    if(hits_from_empLk)
        a += "_hitsFromEmpLk";
    
    a += "_Nag" + int_to_string(Nagents) + "_NEpis" + int_to_string(Nepis) + "_weights";
    
    for(int iag=0;iag<Nagents;iag++)
        a += double_to_string(agent[iag].weight,2) + "_";
        
    if(Nhits==2)
        a += "binary_"; 
    
    if(given_rS.x>0.)
        a += "source" + int_to_string(given_rS.x) + "," + int_to_string(given_rS.y) + "," + int_to_string(given_rS.z) + "_";     
        
    if(given_rA.x>0.)
        a += "agent" + int_to_string(given_rA.x) + "," + int_to_string(given_rA.y) + "," + int_to_string(given_rA.z) + "_";  
      
    a += "n" + int_to_string(n); 
    
    return a;
}

string get_name_params_specific_time(const int Nepis, const int n, const int iep, const int t){
    string a;
    
    a = "dnsData_cThr" + int_to_string(c_thr) + "_windDns" + double_to_string(dns_wind,1) + "_empLk";

    if(hits_from_empLk)
        a += "_hitsFromEmpLk";
    
    a += "_Nag" + int_to_string(Nagents) + "_NEpis" + int_to_string(Nepis) + "_Episode" + int_to_string(iep) + "_t" + int_to_string(t) + "_weights";
    
    for(int iag=0;iag<Nagents;iag++)
        a += double_to_string(agent[iag].weight,2) + "_";
    
    if(Nhits==2)
        a += "binary_";
               
    if(given_rS.x>0.)
        a += "source" + int_to_string(given_rS.x) + "," + int_to_string(given_rS.y) + "," + int_to_string(given_rS.z) + "_";     
        
    if(given_rA.x>0.)
        a += "agent" + int_to_string(given_rA.x) + "," + int_to_string(given_rA.y) + "," + int_to_string(given_rA.z) + "_";  
    
    a += "n" + int_to_string(n);
    
    return a;
}

/* ------------------------------------------------------------------------- */

/* -------- Routines to place agents and source from the likelihood -------- */

// Function to print a randomly selected entry exceeding the threshold
vec3d place_from_hits() {
    vectint exceedingX,exceedingY,exceedingZ;
    vec3d position;
    int size_exceed,randomIndex;

    for (int x = 5; x < NX-5; x++) {
        for (int y = 5; y < NY-5; y++) {
            for (int z = 5; z < NZ-5; z++) {
                if (dnsHits_snapshot[x][y][z] >= c_thr) {
                    exceedingX.push_back(x);
                    exceedingY.push_back(y);
                    exceedingZ.push_back(y);
                }
            }    
        }
    }
    
    size_exceed = exceedingX.size();

    randomIndex = (int) (drand48()*size_exceed);
    randomIndex = min(randomIndex, size_exceed-1);

    position.x = exceedingX[randomIndex];
    position.y = exceedingY[randomIndex];
    position.z = exceedingZ[randomIndex];
    
    exceedingX.clear();
    exceedingY.clear();
    exceedingZ.clear();

    return position;
}

/* ------------------------------------------------------------------------- */

void print_time(const clock_t start_time, const clock_t end_time)
{
  double length = (double)(end_time - start_time)/CLOCKS_PER_SEC;
  
  printf("Total running time = ");
  
  int Nhours = floor(length)/3600;
  
  if(Nhours > 0)
    {
      int Nmin = floor(length - 3600*Nhours)/60;
      printf("%dh%dm%ds\n",Nhours,Nmin,(int)(length - 3600*Nhours - 60*Nmin));
    }
  else
    {
      int Nmin = floor(length)/60;
      if(Nmin > 0)
	printf("%dm%ds\n",Nmin,(int)(length - 60*Nmin));
      else
        {
	  if(length >= 1)
	    printf("%1.2lfs\n",length);
	  else
	    printf("%1.2es\n",length);
        }
    }
}
