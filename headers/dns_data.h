/* ------------ Routines to manipulate dns full time series ---------------- */

string select_dns_data(const int iep, const int times){

    string dns_file;
    
    if(abs(dns_wind-3.)<DBL_EPSILON){
        if(iep<times)  
          dns_file = "../../dns_data/time_series/time_series_data_source1_wind3_x_extent_1548_y_extent_1188_z_extent_1188_nx_129_ny_99_nz_99.bin";
        else if(iep<2*times)
          dns_file = "../../dns_data/time_series/time_series_data_source2_wind3_x_extent_1548_y_extent_1188_z_extent_1188_nx_129_ny_99_nz_99.bin";
        else if(iep<3*times)
          dns_file = "../../dns_data/time_series/time_series_data_source3_wind3_x_extent_1548_y_extent_1188_z_extent_1188_nx_129_ny_99_nz_99.bin";
        else if(iep<4*times)
          dns_file = "../../dns_data/time_series/time_series_data_source4_wind3_x_extent_1548_y_extent_1188_z_extent_1188_nx_129_ny_99_nz_99.bin";
        else if(iep<5*times)
          dns_file = "../../dns_data/time_series/time_series_data_source5_wind3_x_extent_1548_y_extent_1188_z_extent_1188_nx_129_ny_99_nz_99.bin";
    }
    else
        cout << "No such dns time series file in the dns_data folder. Please correct." << endl;
        
    return dns_file;
}

//Import one snapshot of the time series only
void import_dns_snapshot(const string& dns_file, const int t){

  // Calculate the offset based on the size of each matrix and the timestep
  streampos offset = static_cast<streampos>(t) * sizeof(double) * NX * NY * NZ;

  // Read data from the binary file
  ifstream input_file(dns_file, std::ios::binary);
  
  if(!input_file.is_open())
    cerr << "Error opening input file: " << dns_file << endl;

  // Move the file pointer to the desired position
  input_file.seekg(offset, std::ios::beg);
  
  //for(int ix=0;ix<NX;ix++)
  	//input_file.read(reinterpret_cast<char*>(dnsHits_snapshot[ix]), sizeof(int) * NY);
  
  // Read the chunk of data directly into the pre-allocated array  
  input_file.read(reinterpret_cast<char*>(dnsHits_snapshot), sizeof(double) * NX * NY * NZ);

  if(!input_file)
    cerr << "Error reading from input file: " << dns_file << endl;
    
  // Close the file
  input_file.close();
}

//Export one snapshot of the time series only
void export_dns_snapshot(const string& dns_file){

  ofstream output_file(dns_file, std::ios::out | std::ios::binary);
  
  if(!output_file.is_open())
      cerr << "Error opening output file: " << dns_file << endl;
  
  output_file.write(reinterpret_cast<char*>(dnsHits_snapshot), sizeof(double) * NX * NY * NZ);    

  if(!output_file)
      cerr << "Error writing on output file: " << dns_file << endl;
    
  // Close the file
  output_file.close();
}

//Shift snapshot of the time series to adjust to the actual source position
int shift_dns_snapshot() {

    //Shift data depending on where we put the source
    int shiftX = given_rS.x-128;
    int shiftY = given_rS.y-49; 
    int shiftZ = given_rS.z-49; 
    int min_hits = 0;

    //Initialize the temporary array to zeros
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int l = 0; l < NZ; l++)
                temp3D[i][j][l] = 0.;
        }
    }

    // Apply the shift operation to the temporary array
    for (int x = 0; x < NX; x++) {
        for (int y = 0; y < NY; y++) {
            for (int z = 0; z < NZ; z++) {
                int newX = x + shiftX;
                int newY = y + shiftY;
                int newZ = z + shiftZ;

                // Check if the new coordinates are within bounds
                if (newX >= 0 && newX < NX && newY >= 0 && newY < NY && newZ >= 0 && newZ < NZ)
                    temp3D[newX][newY][newZ] = dnsHits_snapshot[x][y][z];
            }
        }   
    }   
    
    for (int x = 0; x < NX; x++) {
        for (int y = 0; y < NY; y++) {
            for (int z = 0; z < NZ; z++) {
                dnsHits_snapshot[x][y][z] = temp3D[x][y][z];
                min_hits = max(min_hits,(int) dnsHits_snapshot[x][y][z]);
            }
        }
    }
    
    return min_hits;
}

//Shift empirical likelihood data 
void shift_emp_lk(){

  //Shift data depending on where we put the source
  int shiftX = given_rS.x-128;  
  int shiftY = given_rS.y-49; 
  int shiftZ = given_rS.z-49; 
      
  //Initialize the temporary array to zeros
  for (int i = 0; i < NX; i++) {
      for (int j = 0; j < NY; j++) {
          for (int l = 0; l < NZ; l++)
              temp3D[i][j][l] = 0.;
      }
  }

  // Read data from the file into the array
  for (int x = 0; x < NX; x++) {
      for (int y = 0; y < NY; y++) {
          for (int z = 0; z < NZ; z++) {
              // Calculate the shifted indices
              int newX = x + shiftX;
              int newY = y + shiftY;
              int newZ = z + shiftZ;
              
              // Check if the shifted indices are within bounds
              if (newX >= 0 && newX < NX && newY >= 0 && newY < NY && newZ >= 0 && newZ < NZ) {
                  // Populate the matrix with the shifted value
                  temp3D[newX][newY][newZ] = dns_empLk[x][y][z];
              }
          } 
      }
  }
  
  // Copy the temporary array back to the original array
  for (int x = 0; x < NX; x++) {
    for (int y = 0; y < NY; y++) {
        for (int z = 0; z < NZ; z++)
            dns_empLk_shifted[x][y][z] = temp3D[x][y][z];
    }
  }  
  
}

/* ------------------------------------------------------------------------- */


/* ------------ Routines to manipulate empirical likelihood ---------------- */

string select_dns_empLk(){
    
    string dns_stationary_file;
    
    if(fabs(dns_wind-3.)<DBL_EPSILON){
        if(fabs(c_thr-5.)<DBL_EPSILON)  
          dns_stationary_file = "../../dns_data/likelihood/emp_likelihood_5sources_wind3_thr5_nx_129_ny_99_nz_99.dat";
        else if(fabs(c_thr-10.)<DBL_EPSILON)  
          dns_stationary_file = "../../dns_data/likelihood/emp_likelihood_5sources_wind3_thr10_nx_129_ny_99_nz_99.dat";
        else if(fabs(c_thr-20.)<DBL_EPSILON)
          dns_stationary_file = "../../dns_data/likelihood/emp_likelihood_5sources_wind3_thr20_nx_129_ny_99_nz_99.dat";
        else{
          cout << "Invalid threshold. Using default threshold for DNS data with wind=3, i.e. c_thr = 10" << endl;
          dns_stationary_file = "../../dns_data/likelihood/emp_likelihood_5sources_wind3_thr10_nx_129_ny_99_nz_99.dat";
        }
    }
    else
        cout << "No such dns empirical likelihood file in the dns_data folder. Please correct." << endl;
         
    return dns_stationary_file;
}

//Import empirical likelihood data
void import_dns_empLk(const string& dns_stationary_file){

  // Open the file for reading.
  ifstream inputFile(dns_stationary_file);
  
  // Check if the file opened successfully
  if (!inputFile.is_open())
    cerr << "Error: Could not open file: " << dns_stationary_file << endl;

  // Read data from the file into the array
  for (int x = 0; x < NX; x++){
      for (int y = 0; y < NY; y++){
          for (int z = 0; z < NZ; z++){
              double value;
              if (inputFile >> value)
                  dns_empLk[x][y][z] = value;
              else {
                  std::cerr << "Error reading data from file: " << dns_stationary_file << endl;
                  inputFile.close();
              }
          }
      }
  }
  
  // Close the file
  inputFile.close();
}

void export_dns_empLk(const std::string& output_file) {

    // Open the file for writing.
    std::ofstream outputFile(output_file);

    // Check if the file opened successfully
    if (!outputFile.is_open()) {
        std::cerr << "Error: Could not open file: " << output_file << std::endl;
        return;
    }

    // Write data from the array to the file as a single list
    for (int x = 0; x < NX; x++) {
        for (int y = 0; y < NY; y++) {
            for (int z = 0; z < NZ; z++) {
                outputFile << dns_empLk[x][y][z] << endl;
            }
        }
    }

    // Close the file
    outputFile.close();

    // Optional: Notify that the export is complete
    std::cout << "Data successfully exported to " << output_file << std::endl;
}

/* ------------------------------------------------------------------------- */
