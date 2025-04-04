/* -------------------------- Routines for input-output ---------------------------*/
void Save_binary_double(const string filename, const double array[], const int size){
    // Open a binary file for writing.
    std::ofstream outfile(filename, std::ios::out | std::ios::binary);

    if (!outfile) {
        std::cerr << "Failed to open the file for writing." << std::endl;
        return;
    }

    // Write the array to the binary file.
    outfile.write(reinterpret_cast<const char*>(array), sizeof(double) * size);

    // Close the file.
    outfile.close();
}

void Save_params(const string name, const int seed, const int iep, const bool append){
  ofstream out;
  if(append)
      out.open(name.c_str(),ios::out|ios::app);
  else{
      out.open(name.c_str(),ios::out);
      out << "# Random seed, source position and agents parameters" << endl;
      out << seed << endl;
      out << source.rS.x << "  " << source.rS.y << "  " << source.rS.z << endl;
  }
  
  if(out.is_open()){
      for(int iag=0;iag<Nagents;iag++)
          out << iep << "  " << iag << "  " << agent[iag].rA.x << "  " << agent[iag].rA.y << "  " << agent[iag].rA.z << "  " << agent[iag].weight << endl;
      out.close();
  }
  else{
      cout << "ERROR: Unable to open the file: " << name << endl;
  }
}


void Save_vector(const string name, const vect vect_variable, const int col, const bool append)
{
  ofstream out;
  
  if(append)
      out.open(name.c_str(),ios::out|ios::app);
  else
      out.open(name.c_str(),ios::out);

  if(out.is_open()){
    int Npoints = vect_variable.size()/col;
    
    for(int i=0;i<Npoints;i++){
        for(int j=0;j<col;j++)
          out << vect_variable[col*i+j] << "   ";
        out << endl;
    }
    out.close();
  }
  else{
    cout << "ERROR: Unable to open the file: " << name << endl;
  }
}

void Save_vector_bin(const string name, const vect v, const int col, const bool append)
{
    // Determine the file open mode based on the specified WriteMode
    std::ios_base::openmode open_mode = std::ios::binary | std::ios::out;
    
    if(append)
        open_mode |= std::ios::app;

    // Open the file with the specified mode
    std::fstream out(name, open_mode);
	//fstream out (name, ios::binary|ios::out);
	out.write((char*)&v[0], v.size()*sizeof(double));
	out.close();

}

void Save_vectorint_bin(const string name, const vectint v, const int col, const bool append)
{
    // Determine the file open mode based on the specified WriteMode
    std::ios_base::openmode open_mode = std::ios::binary | std::ios::out;
    
    if(append)
        open_mode |= std::ios::app;

    // Open the file with the specified mode
    std::fstream out(name, open_mode);
	//fstream out (name, ios::binary|ios::out);
	out.write((char*)&v[0], v.size()*sizeof(int));
	out.close();
}
