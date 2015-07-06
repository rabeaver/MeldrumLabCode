% BAYESIAN_ESTIMATE Estimate components of a T2 distribution
%   [T2 W FLIP PHASE] = BAYESIAN_ESTIMATE(Y,OPTS) Estimates the
%   T2 relaxation times, weights, flip angle and phase of the complex echo
%   data in Y. The Extended Phase Graph (EPG) algorithm is used to jointly
%   estimate the distribution components and the flip angle.
%
%   [T2 W FLIP PHASE COV] = BAYESIAN_ESTIMATE(Y,OPTS) Also returns the
%   covariance of the posterior approximation as a measure of uncertainty.
%
% Author: Kelvin Layton (klayton@unimelb.edu.au)
% Date: June 2012
%
function [T2, w, D, flip, phase, covariance]= bayesian_estimate(y,opts)

% Process inputs parameters
%
sigma2 = opts.sigma2;
tau = opts.tau;
gammas = opts.correctionSchedule;
nCorrections = length(gammas);
priorMu = opts.priorMean;
priorPsi = opts.priorCov;
D = opts.D;
G = opts.G;
decayIn = opts.decayIn;
DELTA = opts.DELTA;
deltaVec = opts.deltaVec;
nDPnts = opts.nDPnts;
% gamma = opts.gamma;	% Gamma, Rad/s/T.

nEchoes = opts.nEchoes;               % Pass nEchoes as a variable of opts
N = 2*opts.nEchoes*nDPnts;     % Number of measurements (real & imag)
nComponents = length(D);    % Number of components
nParams = length(priorMu);


% Check validity of prior
%
% if ~check_prior(priorMu,priorPsi)
%     return;
% end

% Set initial posterior to prior
%
mu_prev = priorMu;
Psi_prev = priorPsi;

% Split measurements into real and imaginary parts
%
%ytilde = [real(y); imag(y)];   

% Iterative over corrections
tic
for j = 1:nCorrections
    
    % Split mean vector into individual parameters
    %
%     f_prev = mu_prev(1:nComponents);                % Weights
%     g_prev = mu_prev(nComponents+1:2*nComponents);  % T2 Guesses
%     d_prev = mu_prev(2*nComponents+1:end-2);        % D Guesses
%     a_prev = mu_prev(end-1);                        % Alpha
%     p_prev = mu_prev(end);                          % Phase
    
      f_prev = (mu_prev(1:nComponents));                % Weights
      g_prev = log10(mu_prev(nComponents+1:2*nComponents));  % T2 Guesses
      d_prev = log10(mu_prev(2*nComponents+1:end-2));        % D Guesses
      a_prev = wrapToPi(mu_prev(end-1));                        % Alpha--FlipAngle
      p_prev = wrapToPi(mu_prev(end));                          % Phase
    
%     decayIn = exp(-gamma.^2 .* G.^2 .* deltaVec.^2 .* (d_prev) .* (DELTA + (2*deltaVec./3)));
    % Get EPG signal (real and imaginary components)
    for l = 1:nDPnts
        Hx(:,l) = epg_diff(nEchoes,tau,1,10^g_prev,a_prev, (10^d_prev), G, deltaVec(l), DELTA) * exp(1j*p_prev);
    end
%     H = Hbase * ones(size(decayIn)); %decayIn;
    H = reshape(Hx, N/2,1);
    s_prev = [real(H*f_prev); imag(H*f_prev)];
    s_prev = reshape(s_prev, N,1);
    
    % Compute derivatives of each component
    %
    [Jg, Jd, Ja] = epg_diff_derivatives(nEchoes,tau,1,10^g_prev,a_prev, 10^d_prev, G, deltaVec(l), DELTA);  % Needs to return Jd
    
    % Calculate Jacobian by combining components and adding phase term
    %
    Jf = H;
    Jd = (exp(1j*p_prev).*bsxfun(@times,Jd, f_prev')) * ones(size(decayIn)); % * decayIn;
    Jd = reshape(Jd, N/2, 1);
    Jg = (exp(1j*p_prev).*bsxfun(@times,Jg, f_prev')) * ones(size(decayIn));
    Jg = reshape(Jg, N/2, 1);
    Ja = (exp(1j*p_prev).*(Ja*f_prev)) * ones(size(decayIn));
    Ja = reshape(Ja, N/2, 1);
    Jp = 1j.*(H*f_prev); 
    
    J = [real([Jf Jg Jd Ja Jp]); imag([Jf Jg Jd Ja Jp])];

%     figure(2)
%     subplot(2,1,1)
%     plot(real(Jf))
%     subplot(2,1,2)
%     plot(real(Jd))
    
    % Compute matrices for this correction step
    % ScaledSigma = sigma2/gammas(j)*eye(N);
    %
%     ScaledSigma_inv = gammas(j)/sigma2*eye(N);
%     Psi_prev_inv = inv(Psi_prev);
%     
%     % We want inverse as: Psi_inv = inv(J*Psi_prev*J' + ScaledSigma);
%     % Since ScaledSigma is diagonal, we use the Matrix Inversion Lemma
%     % otherwise known as Woodbury matrix identity. This is much faster.
%     %
%     B = Psi_prev_inv + gammas(j)/sigma2*(J'*J);
%     C = ScaledSigma_inv*J;
%     D = J'*ScaledSigma_inv;
%     A = C/B*D;
%     % THIS IS THE SLOW STEP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% %     Psi_inv = ScaledSigma_inv - ScaledSigma_inv*J*inv(Psi_prev_inv + gammas(j)/sigma2*(J'*J))*J'*ScaledSigma_inv;
%     Psi_inv = ScaledSigma_inv - A;
    


    % Perform posterior update
%     %
%     figure(5)
%     subplot(5,1,1)
%     hold off
%     plot(J(:,1)); legend('w'); 
%     hold on
%     plot(diff(J(:,1)));
%     subplot(5,1,2)
%     hold off
%     plot(J(:,2)); legend('T2');
%     hold on
%     plot(diff(J(:,2)));
%     subplot(5,1,3)
%     hold off
%     plot(J(:,3)); legend('D');
%     hold on
%     plot(diff(J(:,3)));
%     subplot(5,1,4)
%     hold off
%     plot(J(:,4)); legend('flip');
%     hold on
%     plot(diff(J(:,4)));
%     subplot(5,1,5)
%     hold off
%     plot(J(:,5)); legend('phase');
%     hold on
%     plot(diff(J(:,5)));
%     
    ScaledSigma = sigma2*eye(N);
    Kf = Psi_prev*J'/((J*Psi_prev*J'+ ScaledSigma)/gammas(j));
    L = Kf*(y-s_prev); %amount by which the guesses will change for each step in the corrections
%     M = y-s_prev; %difference between data and current guess in corrections schedule
    mu = mu_prev + L;
    Psi = (eye(nParams) - Kf*J)*Psi_prev;

    % Save current posterior for next iteration
    %
    Psi_prev = Psi;
    mu_prev = mu;
end
toc
% Extract individual estimates form state vector
%
T2 = 1./mu(nComponents+1:2*nComponents);    % Convert rates to times

if nargout>1
    w = mu(1:nComponents);  % Component weights
end
if nargout>2
    D = mu(2*nComponents+1:end-2);  % Component weights
end
if nargout>3
    flip = mu(end-1); % Flip angle
end
if nargout>4
    phase = mu(end);   % Phase
end
if nargout>5
    covariance = Psi;   % Covariance matrix 
end

%%
% Helper function to check the validity of the user specified prior
%
% In addition to basic checks, certain conditions are problematic for the
% estimation algorithm and this function serves to identify these
% situations. These problematic conditions are easy to avoid in practice.
%
function ok = check_prior(priorMu,priorPsi)
    nPar = length(priorMu);    % Number of parameters to estimate
    nComp = (nPar-2)/2;        % Number of components
    
    % Check size of mean and covariance
    %
    assert(mod(nComp,1)==0,'Prior must specify weight and location of each mode followed by flip angle and phase');
    assert(length(priorMu)==size(priorPsi,1),'Prior mean and covariance sizes must match');

    % Extract prior means
    %
    priorRates = priorMu(nComp+1:2*nComp);
    priorAlpha = priorMu(2*nComp+1);
    priorPhase = priorMu(end);
    
    % The prior mean for the flip angle should be away from 180 degrees. 
    % At 180 degrees the derivative of the signal is zero and the posterior
    % approximation has trouble updating.
    %
    assert(priorAlpha < pi,'Prior mean for flip angle should be less than 180 degrees');
    
    % Phase only needs to be between -pi and pi to represent the full range
    % of angles.
    %
    assert(priorPhase>=-pi && priorPhase <= pi,'Prior mean for phase should be between -180 and +180 degrees');   
    
    % Prior covariance must not be singular otherwise initial inversions
    % fail.
    %
    assert(cond(priorPsi)<1e9,'Prior covariance is close to singular');

    % Check the prior doesn't have identical means for the relaxation rates 
    % of different components. This will not work well due to the 'label
    % switching' problem of mixture estimation.
    %
    perms=nchoosek(1:nComp,1);
    rateDiff = 1;
    assert(all(rateDiff>0.01),'Prior means for rates must be sufficiently separated');
    
    ok=true;
end

end