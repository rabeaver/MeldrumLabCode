function data4 = loaddata(filename,startint,endint)
    data=load(filename,'-ascii');
    data2=data(:,3);
    np = max(size(sattime));
    pts = max(size(data2(:,1)))/np;
    data3=reshape(data2,pts,np);
    data4=squeeze(sum(data3(startint:endint,:)));
    data4 = data4';
end