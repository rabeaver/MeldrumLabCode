function [fheader]=readheader(f,n)

     %f: file pointer to read from
   %n: number of points at the end of fid to use for baseline correction
   %(-1 = no correction)
   %fheader: returned fid file header information

   fheader=struct('nblocks',{0},'ntraces',{0},'np',{0},'ebytes',{0},'tbytes',{0},'bbytes',{0},'vers_id',{0},'status',{0},'nb_headers',{0});
   bheader=struct('scale',{0},'status',{0},'index',{0},'mode',{0},'ctcount',{0},'lpval',{0},'rpval',{0},'lvl',{0},'tlt',{0});

   %Read vnmr file

   fheader.nblocks=fread(f,1,'long');
   fheader.ntraces=fread(f,1,'long');
   fheader.np=fread(f,1,'long');
   fheader.ebytes=fread(f,1,'long');
   fheader.tbytes=fread(f,1,'long');
   fheader.bbytes=fread(f,1,'long');
   fheader.vers_id=fread(f,1,'short');
   fheader.status=fread(f,1,'short');
   fheader.nb_headers=fread(f,1,'long');
 
end
