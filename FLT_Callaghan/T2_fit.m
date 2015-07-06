    load T2_data.dat
    
    t = T2_data(:,1)*0.00015;
    Expt_t2 = T2_data(:,5);
    func_t2 = inline('y(3) + y(1)*exp(-t*y(2))','y','t');
    y = [1.0 0.005 1.0];
    [y,r1,j1] = nlinfit(t,Expt_t2,func_t2,y);
    I0 = y(3);
    A  = y(1);
    T2 = y(2);
    Cal_t2 = y(3) + y(1)*exp(-t*y(2));
    
    plot(t,Expt_t2,'*k')
    hold
    plot(t,Cal_t2)