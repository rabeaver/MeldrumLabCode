function p=readpar_Kea(filename,parname)

   %reads a parameter from vnmr parameter file
   %filename= path to vnmr procpar file
   %parname= name of parameter to read
   %p= returned parameter value
   
   f=fopen(filename,'r','b');
   line=fgetl(f);
   while (~strncmp(line,parname,length(parname)))&&(~feof(f))
       line=fgetl(f);
   end
   fclose(f);

   
   if line == -1
       p=1;
   else
       q = line(length(parname) + 4:length(line));
       p = str2double(q);
   end
   
   if strcmp(q,'"yes"')
       p = 1;
   elseif strcmp(q,'"no"')
       p = 0;
   end
   
end