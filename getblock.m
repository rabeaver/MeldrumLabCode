function [dta,bheader]=getblock(file,np,n)

   %file: file to read block from
   %dta: returned block
   %bheader: returned fid block header information for the last block that
   %was read

   bheader=struct('scale',{0},'status',{0},'index',{0},'mode',{0},'ctcount',{0},'lpval',{0},'rpval',{0},'lvl',{0},'tlt',{0});

   %Read vnmr file
   dta=[];

   bheader.scale=fread(file,1,'short');
   bheader.status=fread(file,1,'short');
   bheader.index=fread(file,1,'short');
   bheader.mode=fread(file,1,'short');
   bheader.ctcount=fread(file,1,'long');
     bheader.lpval=fread(file,1,'float');
     bheader.rpval=fread(file,1,'float');
     bheader.lvl=fread(file,1,'float');
     bheader.tlt=fread(file,1,'float');
     d=fread(file,np,'int32');

     d=transpose(reshape(d,2,size(d,1)/2));
     d=complex(d(:,1),d(:,2));
     if n>=0
        d=d-mean(d(end-n:end));
        dta=[dta,d];
     end

end
