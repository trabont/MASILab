function [x,chi2,chi2p,res,resn,out,exit] = fit_SSPulseMT_yarnykh_Full_Fit(x0,Mdat,pwMT,ts,TR,R1obs,theta,w1e,delta,lineshape)
%% Function to Calculate qMT using Yarnykh Full Fit Model
 
% Fit data - x0 = [M0b kba T2a T2b]
opt = optimset('Display','off','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',2000);
[x,resn,res,exit,out,lam,J] = ...
    lsqnonlin(@fitModel,x0,[1e-3 1e-3 1e-4 1e-9],[1 50 1 1]...
    ,opt,Mdat,pwMT,ts,TR,R1obs,theta,w1e,delta,lineshape);

% Calculate 95% CIs
% ci = nlparci(x,res,'jacobian',J);

chi2 = res.^2./Mdat;
chi2 = sum(chi2(:));
chi2p = chi2cdf(chi2,2*length(delta)-1);

%--------------------------------------------------------------------------
function res = fitModel(x,Mdat,pwMT,ts,TR,R1obs,thetaEX,w1e,delta,lineshape)

% Set model parameters
% Richard's Code has M0(1) as 1, not as a difference, this changes the code
% slightly and makes our answers differ slightly, but they are equivalent
% when changed to reflect his code using my model

% [(1-x(1) x(1)] - MPF, [1 x(1)] - PSR
M0 = [1 x(1)]; 
kba = x(2); 
T2(1) = x(3);
T2(2) = x(4);


% Solve for R1a and R1b
% R1(2) = 1;
% % Get R1a based upon Robs and current model parameters
% R1(1) = R1obs - ((kab*(R1(2)-R1obs))./(R1(2)-R1obs+kab*M0(1)/M0(2)));

% In Preparation for Move to 7T, Changed to R1a=R1b=R1obs, as in Yarnykh,
% 2012
R1(1) = R1obs;
R1(2) = R1(1);


% Calculate model

[Mact,Mnorm] = yarnykh_pulseMT(M0,R1,T2,TR,kba,pwMT,ts,thetaEX,delta,w1e,lineshape);
Mact = double(Mact);
Mnorm = double(Mnorm);

if size(Mnorm,2) ~= size(Mdat,2)
    Mnorm = Mnorm';
end

% Calculate residuals
res = Mnorm - Mdat;

% % Used to Check Fits
% figure(10);
% semilogx(delta,Mdat,'o'), hold on
% semilogx(delta,Mnorm), drawnow, hold off

