function dimensionlessparameters=define_parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define parameters in dimensional form %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The parameter values set here are those used to produce the data
% point for extraction corresponding to the cafe recipe documented 
% in the paper with grounds processed at a 1.5 setting

R0=29.2e-3;
L=18.7e-3;                       

a1=49.75e-6;
a2=335.41e-6;

rhoout=997;
rhogrounds=330;

csat=2.124e2;
cs0=1.18e2;

phis=0.8272;

Mout=0.04;

P=5;

Ds_star=0.625e-9;
Deff=0.01;

k=6e-7;

tshot=33.9;

% The vol. fracs of fines and boulders (of the grounds)
f1=0.1624;
f2=1-f1;

% Calculate vol. fracs. of fines and boulders (of baskets)
phis1=phis*f1;
phis2=phis*f2;

bet1_star=3*phis1/a1;
bet2_star=3*phis2/a2;
bet0_star=bet1_star/2+bet2_star/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert to dimensionless parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bet1=bet1_star/bet0_star;
bet2=bet2_star/bet0_star;
Ds1=Ds_star*tshot/a1^2;
Ds2=Ds_star*tshot/a2^2;
Q1=1/a1/bet0_star;
Q2=1/a2/bet0_star;
beta=csat/cs0;
q=Mout/pi/R0^2/rhoout/L;
alpha=cs0/rhogrounds;

K=k*cs0^2*tshot*bet0_star;

dimensionlessparameters=[bet1 bet2 Ds1 Ds2 Q1 Q2 beta q Deff K alpha phis];