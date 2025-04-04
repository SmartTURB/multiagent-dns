//Renormalize belief with hole in ix,iy,iz
void renormalize_belief(double *p,int ix,int iy,int iz){
  double norm,appo;
  double *a;
  int i,k;
  
  //Adaptable to other found conditions -- now good when found at dist < 1
  appo=0.;
  for(i=6;i<7;i++){
      k = (min(max(ix+step[i].x,0),NX-1)*NY + min(max(iy+step[i].y, 0),NY-1))*NZ + min(max(iz+step[i].z, 0),NZ-1);
      appo+=p[k];
      p[k]=0.0;
  }
  
  norm=1.0/(1.0-appo);
  for(a=p;a<p+NS;++a)
      *a *= norm;
}

//Update belief and compute useful quantities
double BayesUpdate(int ax, int ay, int az, double *p, int nh, double *entropy, double *dist){
  double norm,appo,pp,empiric_lk;
  int i,j,l;

  //Source position in the DNS file
  int sourceDns_x = 128;
  int sourceDns_y = 49;
  int sourceDns_z = 49;
      
  //Update belief
  norm=0.0;
  for(i=0;i<NX;i++){
      for(j=0;j<NY;j++){
          for(l=0;l<NZ;l++){
              if(ax + sourceDns_x - i >= NX || ax + sourceDns_x - i < 0 || ay + sourceDns_y - j >= NY || ay + sourceDns_y - j < 0 || az + sourceDns_z - l >= NZ || az + sourceDns_z - l < 0)
                  empiric_lk = (nh==1) ? 0. : 1.;
              else
                  empiric_lk = (nh==1) ? dns_empLk[ax+sourceDns_x-i][ay+sourceDns_y-j][az+sourceDns_z-l] : 1. - dns_empLk[ax+sourceDns_x-i][ay+sourceDns_y-j][az+sourceDns_z-l]; 
       
              p[(i*NY+j)*NZ+l] *= empiric_lk;  
              norm+=p[(i*NY+j)*NZ+l];
          }
      }
  }
  
  //Normalize belief and compute its entropy (for infotaxis) and Manhattan distance (for greedy policy)
  appo=1./norm;
  *entropy=0.0;
  *dist=0.0;
  for(i=0;i<NX;i++){
      for(j=0;j<NY;j++){
          for(l=0;l<NZ;l++){  
              pp=p[(i*NY+j)*NZ+l]*=appo;      
              if(pp>0){
            	  *entropy-=pp*log(pp);
            	  *dist+=(abs(i-ax)+abs(j-ay)+abs(l-az))*pp;
              }
          }
      }
  }
    
  *entropy/=log(2);
  
  return norm;
}
