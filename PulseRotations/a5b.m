

% Bloch Equation Simulation
% -----------------------------------------
% 

M0 = [0;0;1];
a = pi/2;

M_p1 = xrot(a)*M0;

%

dT = 1;		% 1ms time step.
T = 100;	% 100ms total duration
N = ceil(T/dT)+1; % number of time steps.
df = 20;	% Hz offset.
T1 = 200;	% ms.
T2 = 30;	% ms.

% ===== Get the Propagation Matrix ======

[A,B] = freeprecess(dT,T1,T2,df);


% ===== Simulate the Decay ======

M = zeros(3,N);	% Keep track of magnetization at all time points.
M(:,1)=M_p1;	% Starting magnetization.

for k=2:N
	M(:,k) = A*M(:,k-1)+B;
end;



% ===== Plot the Results ======

time = [0:N-1]*dT;
plot(time,M(1,:),'b-',time,M(2,:),'r--',time,M(3,:),'g-.');
legend('M_x','M_y','M_z');
xlabel('Time (ms)');
ylabel('Magnetization');
axis([min(time) max(time) -1 1]);
grid on;



