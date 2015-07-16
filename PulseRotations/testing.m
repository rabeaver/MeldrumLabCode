M0 = [0;0;1];
T2 = 100e3; %us
T1 = 600e3; %us

pulseSeq = [ 5   1    0     pi/2;   %delay (us)   rotx   roty   flip
             50  0    0     0;
             5   0    1     pi;
             200 0    0     0];

         
totTime = sum(pulseSeq(:,1));
taxis = linspace(0,totTime,totTime*10+1);
M = zeros(3,length(taxis));
M(:,1) = M0;

nrSeg = size(pulseSeq,2);

for i = 1:nrSeg
    p(i,:) = [pulseSeq(i,1)*10 pulseSeq(i,2) pulseSeq(i,3) pulseSeq(i,4)/10];
end

for tt = 2:length(totTime)

    if pulseSeq(seg,2) == 1;
        for i = 2:pulseSeq(seg,1)*10+1;
            M(:,i) = xrot(pulseSeq(seg,4)/pulseSeq(seg,1)/10)*M(:,i-1);
        end
    elseif pulseSeq(seg,3) == 1;
        for i = 2:pulseSeq(seg,1)*10+1;
            M(:,i) = yrot(pulseSeq(seg,4)/pulseSeq(seg,1)/10)*M(:,i-1);
        end
%     elseif pulseSeq(seg,2)+pulseSeq(seg,3) == 0;
%         freeprec
    end
        
end  

% ===== Plot the Results ======

time = 1:52;
plot(time,M(1,1:52),'b-',time,M(2,1:52),'r--',time,M(3,1:52),'g-.');
legend('M_x','M_y','M_z');
xlabel('Time (ms)');
ylabel('Magnetization');
axis([min(time) max(time) -1 1]);
grid on;