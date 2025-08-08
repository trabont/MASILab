function [PSR,R1obs,chi2,chi2p,res,resn] = Analysis_Yarnykh_1pt(p)
%% Analysis_Yarnykh_1pt - Finds T1 and Runs the Fitting function for the SP Analysis
% Finds T1 using the MFA method from Haacke, then prepares the data to run
% the fitting function for finding the SP qMT
%
% Syntax:  [PSR,R1obs,chi2,chi2p,res,resn] = Analysis_Yarnykh_1pt(p)
%
% Inputs:
%   p - Structure which needs the following inputs:
%       p.B1 - B1 correction data (fractional)
%       p.B0 - B0 correction data
%       p.MFA - Number of dynamics in the MFA image
%       p.T1flip - FA from the Examcard in degrees (calculates rest of FA in code)
%       p.Ernst - T1 Data (in order from largest FA to smallest FA)
%       p.T1TR - TR from T1 data (in ms)
%       p.pwMT - pulse length of MT pulse (in s)
%       p.MT_flip - Flip angle of MT saturation pulse (in degrees)
%       p.qMTflip - Flip angle of Excitation pulse (in degrees)
%       p.M - MT data
%       p.TR - TR of MT data (in ms)
% 
% Outputs:
%   PSR - PSR of the MT data
%   R1obs - R1obs (in 1/s)
%   chi2 - chi2 goodness of fit
%   chi2p - p-value of chi2 goodness of fit
%   res - residuals
%   resn - normalized residual from lsqnonlin function
%




%% MFA - R1obs - (Ke's Code)
corrB1 = double(p.B1/100);


T1flip = double(p.T1flip);
MFA = double(p.MFA);
Ernst = double(p.Ernst);
T1TR = double(p.T1TR);

del_flip = T1flip/MFA;
thetaT1 = T1flip:-del_flip:del_flip; % deg
yval = Ernst./sind(thetaT1*corrB1);
xval = Ernst./tand(thetaT1*corrB1);
pfit=polyfit( xval ,yval, 1);

E1 = pfit(1);

T1obs = real(-T1TR/log(E1));
R1obs = 1/T1obs*1000; % s

% Set R1obs for ON Grant
% R1obs = 1;

if R1obs < 0 || isnan(R1obs) || isinf(R1obs) || isinf(T1obs)
    PSR = 0;
    R1obs = 0;
    chi2 = 0;
    chi2p = 0;
    res = 0;
    resn = 0;
    return
end

%% qMT Data Prep

pwMT = double(p.pwMT);
MT_flip = double(p.MT_flip);
qMTflip = double(p.qMTflip);
deltaMT = double(p.deltaMT);
M = double(p.M);
TR = double(p.TR);

gamma = 42.58*2*pi; % Larmor - rad/s-uT


B1ampMT = MT_flip./(gamma*pwMT); % uT
[B1MT,tMT] = philipsRFpulse_FA(MT_flip,pwMT,'am_sg_100_100_0');
B1eMT = CWEqMTPulse(B1MT*corrB1,tMT,pwMT);

thetaEX = ([qMTflip qMTflip]*pi/180)*corrB1;
ts = 1e-3; % s
corrB0 = double(p.B0);
TR = [TR*10^-3, TR*10^-3];

%% 1 Parm Fit

% Initial estimate
%       PSR   kba   T2b    T2aR1a
% p0 = [0.15  10   10e-6   10e-3];

p0 = 0.15;
p.kba = 10;
p.T2b = 10e-6;
p.T2aR1a =  10e-3;

[PSR,chi2,chi2p,res,resn] = fit_SSPulseMT_yarnykh_1pt(p0,M,pwMT,ts,TR,R1obs,thetaEX,B1eMT,deltaMT+corrB0,'super-lorentzian',p.kba,p.T2b,p.T2aR1a);




