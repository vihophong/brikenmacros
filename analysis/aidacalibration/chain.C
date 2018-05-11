void chain(char* listfile){
  ch = new TChain("aida");
  std::ifstream ifs(listfile);
  string filelist[1000];

  Int_t nfiles=0;
  while (!ifs.eof()){
      ifs>>filelist[nfiles];
      cout<<filelist[nfiles]<<endl;
      nfiles++;
  }
  nfiles=nfiles-1;
  cout<<"There are "<<nfiles<<" files in total!"<<endl;

  for (Int_t i=0;i<nfiles;i++){
      char tempchar[1000];
      sprintf(tempchar,"%s/aida",filelist[i].c_str());
      ch->Add(tempchar);
  }
}
