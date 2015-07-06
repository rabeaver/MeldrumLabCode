function p=readpar(filename,parname)

   %reads a parameter from vnmr parameter file
   %filename= path to vnmr procpar file
   %parname= name of parameter to read
   %p= returned parameter value
   
   parname = strcat('##$',parname);
   
   f=fopen(filename,'r','b');
   line=fgetl(f);
   while (~strncmp(line,parname,length(parname)))&&(~feof(f))
       line=fgetl(f);
   end
   fclose(f);
   
   if line == -1
       p=1;
   else
   
   param = line;
   param = param(size(parname,2)+2:size(param,2));
   param = str2double(param);
   p=param;

   end

end