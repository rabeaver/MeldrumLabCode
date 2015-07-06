function p=readpar(filename,parname)

   %reads a parameter from vnmr parameter file
   %filename= path to vnmr procpar file
   %parname= name of parameter to read
   %p= returned parameter value
   
   f=fopen(filename,'r','b');
   line=fgetl(f)
   while (~strncmp(line,parname,length(parname)))&&(~feof(f))
       line=fgetl(f)
   end
   line=fgetl(f)  
   fclose(f);
   
   if line == -1
       p=1;
   else
       param = str2num(line);
       param = param(2:(size(param,2)));
       p=param';
   end
   
end