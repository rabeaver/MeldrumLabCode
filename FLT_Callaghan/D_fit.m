    load D_data.dat
    
    smallD = 0.0005;
    bigD = 0.002;
    gamma = 4.006e+03;
    
    for i = 1:1:32;
    k(i,1) = ((2*pi*gamma*D_data(i,2)*100*smallD)^2)*(bigD - smallD/3);
    end;
    
    Expt_D = D_data(:,3);
    func_D = inline('y(1)*exp(-k*y(2))','y','k');
    y = [1.0 1.0e-6];
    [y,r1,j1] = nlinfit(k,Expt_D,func_D,y);
    p = y(1);
    D = y(2);
    Cal_D = y(1)*exp(-k*y(2));
    
    plot(k,Expt_D,'*k');
    hold;
    plot(k,Cal_D);