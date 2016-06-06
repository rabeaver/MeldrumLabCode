clear all;
pwd1 = pwd;

% load of the two distribution files
[first,c1] = uigetfile('*.out','                                                              first data');
cd(c1)
firstdis = load(first,'-ASCII');
[last,c1] = uigetfile('*.out','                                                              last data');
lastdis = load(last)

% load of the two horizontal and vertical "time"files.
f2 = uigetfile('*.out','                                                          horizontal "time"');
horiz = load(f2)
f3 = uigetfile('*.out','                                                          vertical "time"');
vert = load(f3)
cd(pwd1)

[sh,nimportequoi]=size(horiz);
[sv,nimportequoi]=size(vert);

% weightfirst = sum(sum(firstdis))
% weightlast = sum(sum(lastdis))

% transformation of the matrices (Laplace shift and flip)
q02 = inputdlg('enter q0_square value', 'q0_square shift information',1,{'0.5'})
q02 = str2num(q02{1})
for j = 1:sv
    lastdis(j,:) = lastdis(j,:)*exp(-q02*horiz(j));
end
lastdis = fliplr(lastdis);
firstdis = fliplr(firstdis);
horiz = 1./horiz;
horiz = flipud(horiz);
horizl = log10(horiz);
vertl = log10(vert);




% distribution = weightfirst.*firstdis + weightlast.*lastdis
distribution = firstdis + lastdis;

% figure of the first distribution
figure
surf(horizl, vertl, firstdis)
axis([horizl(1),horizl(sh),vertl(1),vertl(sv)]);
shading interp;
set(gcf,'Renderer','zbuffer');
title('first distribution')

% figure of the last distribution
figure
surf(horizl, vertl, lastdis)
axis([horizl(1),horizl(sh),vertl(1),vertl(sv)]);
shading interp;
set(gcf,'Renderer','zbuffer');
title('last distribution')

% figure of the finale distribution
figure
surf(horizl, vertl, distribution)
axis([horizl(1),horizl(sh),vertl(1),vertl(sv)]);
shading interp;
set(gcf,'Renderer','zbuffer');
title('finale distribution')

distributionl = log(distribution+1)
figure
surf(horizl, vertl, distributionl)
axis([horizl(1),horizl(sh),vertl(1),vertl(sv)]);
shading interp;
set(gcf,'Renderer','zbuffer');
title('finale distribution')

% saving of the new distribution and vertical and horizontal files
cd(c1);
[nimportequoi,tailletexte] = size(first);
filename = first;
filename(tailletexte-8:tailletexte) = [];
savefile = [filename 'average.out'];
save(savefile,'distribution','-ASCII');

[nimportequoi,tailletexte] = size(f2);
filename = f2;
filename(tailletexte-3:tailletexte) = [];
savefile = [filename 'average.out'];
save(savefile,'horiz','-ASCII');

[nimportequoi,tailletexte] = size(f3);
filename = f3;
filename(tailletexte-3:tailletexte) = [];
savefile = [filename 'average.out'];
save(savefile,'vert','-ASCII');

filename = first;
filename(tailletexte-8:tailletexte) = [];

savefile = [filename 'average.fig'];
saveas(gcf,savefile);
savefile = [filename 'average.jpg'];
saveas(gcf,savefile);


