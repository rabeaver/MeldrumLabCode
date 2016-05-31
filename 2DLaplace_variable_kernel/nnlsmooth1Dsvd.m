function [spectrum, chisq,compte]=nnlsmooth1Dsvd(spectrum,weighting, tau, time, data, kernel1)


% conditions to end the program
repeat = 1;
fini = 0;

% definition of the size of the data and the result
[~,tailledata] = size(data);
[~,Xdim] = size(tau);

% step #1 of the L&H algorithm: initialization of the vectors

% Make Designmat and systemat
% matrice of the exp(-time/tau). from Tailledata + 1 to tailledata+Xdim-2: includes the weighting, includes data
% this matrix will be diagonalized again and again 
% matrice of the exp(-time/tau). from Tailledata + 1 to tailledata+Xdim-2: includes the weighting

systemmat = zeros(tailledata,Xdim);
% designmat = zeros(tailledata,Xdim);

kernel1=vectorize(strtrim(strrep(strrep(strrep(kernel1,'h','time'''),'D','(1/tau(a))'),'T','tau(a)')));
func1=inline(kernel1,'a','time','tau');

for a = 1 : Xdim
    systemmat (:,a) = func1(a,time,tau);
end

% ssyb = size(systemmat);
% svd decomposition to reduce the size of the data

[u,~,~] = svd(systemmat,0);

systemmat = u'*systemmat;
data = u'*data';

[~,tailledata] = size(systemmat);
conct = zeros(Xdim-2,Xdim);
systemmat = [systemmat;conct];

for i = 1 : Xdim-2
    systemmat(i+tailledata,i) = weighting(i+1);
    systemmat(i+tailledata,i+1) = -2*weighting(i+1);
    systemmat(i+tailledata,i+2) = weighting(i+1);
end

designmat = systemmat;

conct = zeros(Xdim-2,1);
datac = [data;conct];

% including data into systemmat
systemmat = [systemmat datac];

% ssya = size(systemmat)

designmatT=designmat';      % tranpose of the matrix designmat for the step #2 of the L&H algorithm


chivector = zeros(tailledata,1);
% chivector: at the end: product between designmatT and spectrum
% (chivectro-dat)^2 = chisq;  same dimension than data
answ = zeros(Xdim + tailledata -2,1);
% answ = f-Ex, dimension Xdim + tailledata -2
% oldsetp = zeros(1,Xdim);
setp = zeros(1,Xdim);
% set of P indices (at the beginning "empty")
soln = zeros(1,Xdim);
% soln: solution of the linear eq with the diagonalized systemmat
% W = zeros(Xdim,1);
% W = E'(f-Ex), same dimension as sepctrum (Xdim), needed for step #2
% spectrum = zeros(1,Xdim); %distribution (result of the problem)
setz = zeros(1,Xdim);
setz = setz + 1; 
% 1 indicates position filled, -1 indicates position empty


% step # 2 L&H
% first stage: spectrum is all zero
% answ = f-Ex
for i = 1 : tailledata
    answ(i) = data(i) - answ(i);
end
W = designmatT*answ;

% step #4 L&H
% Find an index t in setz such that W(t) is the max of W
% for the first step: no constraint
[~,indice]=max(W);
compte = 0;
while fini == 0
    compte = compte + 1;
    % step #5 of the L&H algorithm
    % move the index "indice" from the set setz to the set setp
    oldsetp = setp;
    ind = find(setp == 0);
    setp(ind(1))=indice;
    setz(indice)=-1;
    [systemmat] = qrdecmp(systemmat,setp,oldsetp,Xdim+1,tailledata+Xdim-2);
    [setp,soln] = lineareqsolve(systemmat,setp,soln,Xdim+1);
    
    if soln(indice) < 0
        warning('correction made for rounoff error');
    end
    
    % only if soln(indice) < 0, find new indice
    while soln(indice) < 0
        W(indice) = 0;

        biggest = 0;
        newindice = 0;
        for i = 1 : Xdim
            if W(i) > biggest && setz(i) == 1
                biggest = W(i);
                indice = i;
            end
        end        
        if newindice == 0
            fini = 1;
            break
        end
        setz(indice)=1;
        setz(newindice)=-1;
        indice = newindice;
        ind = find(setp == 0);
        setp(ind(1)-1)=indice;     
        [systemmat] = qrdecmp(systemmat,setp,oldsetp,Xdim+1,tailledata+Xdim-2);
        [setp,soln] = lineareqsolve(systemmat,setp,soln,Xdim+1);
    end
    
    if fini == 1
        break;
    end
    repeat = 1;
    
    % step #6 L&H
    % step #7: if Zj >0 for all j in setp, set x=z (spectrum = soln) and go to step #2, else ...
    while repeat == 1
        [lambda,indextozero] = findlambda(spectrum,soln,setp,Xdim);
        % step #8 #9. find indextozero=q in setp, such that xq/(xq-zq)=min.
        % If there is a negative value in soln, lambda won't be zero
        if indextozero ~= 0
            %if there is a negative value in soln, we have to keep removingindices from setp until all elt in soln are positive
            if lambda == 0      % not necessary
                warning('algorithm fault - infinite loop');
                repeat = 0;
                fini = 1;
                break;
            end
            % step #10 set x = x + lambda(z-x)
            spectrum = spectrum + lambda.*(soln-spectrum);
            % steps #11. Move from setp to setz all indices j in setp for which xj=0. go to step #6
            oldsetp = setp;
            [setp,setz]=deletefromset(setp,setz,indextozero,Xdim);
            spectrum(indextozero)=0;
            soln(indextozero)=0;
            %step # 6
            [systemmat] = qrdecmp(systemmat,setp,oldsetp,Xdim+1,tailledata+Xdim-2);     % Get to step 6 in the algorithm
            [setp,soln]=lineareqsolve(systemmat,setp,soln,Xdim+1);                      % until no more soln is negative
        else %Step #7 : x=z, go to step #2
            repeat = 0;
            spectrum = soln;
        end
    end
    if fini == 1
        break;
    end
    
    % step #2 of the L&H algorithm. first step: spectrum is all zero
    answ = designmat*spectrum';
    for i = 1 : tailledata
        answ(i)=data(i)-answ(i);
    end
    % w = E'(f-Ex)
    for i = tailledata+1 : tailledata+Xdim-2
        answ(i)=-answ(i);
    end
    W = designmatT*answ;
    
    % step #4 of the L&H algorithm
    % find index of bigger number in W, whose index is still in setz
%     indice = 0;
    biggest = 0;
    indice = 0;
    for i = 1 : Xdim
        if W(i) > biggest && setz(i) == 1
            biggest = W(i);
            indice = i;
        end
    end
    
    if indice == 0
        fini = 1;
    end
    if fini ~= 1 && setp(Xdim) ~= 0
        fini = 1;
    end
    
end

for i = 1 : tailledata
        chivector(i) = designmat(i,:)*spectrum';
end

chisq = sum((chivector-data).^2);
