    load T1_data.dat
    
    t = T1_data(:,1)*0.6;
    Expt_t1 = T1_data(:,5);
    func_t1 = inline('y(3) + y(1)*exp(-t*y(2))','y','t');
    y = [-1.0 1.0 1.0];
    [y,r1,j1] = nlinfit(t,Expt_t1,func_t1,y);
    I0 = y(3);
    A  = y(1);
    T1 = y(2);
    Cal_t1 = y(3) + y(1)*exp(-t*y(2));
    
    plot(t,Expt_t1,'*k')
    hold
    plot(t,Cal_t1)
  