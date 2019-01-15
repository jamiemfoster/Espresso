clear all
clf
close all
clc

tic

plot_flag=1; % Uses figures 1-4
checkmass_flag=1; % Uses figure 5

N=40;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the dimenionless parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dimensionlessparameters=define_parameters;
bet1=dimensionlessparameters(1);
bet2=dimensionlessparameters(2);
Ds1=dimensionlessparameters(3);
Ds2=dimensionlessparameters(4);
Q1=dimensionlessparameters(5);
Q2=dimensionlessparameters(6);
beta=dimensionlessparameters(7);
q=dimensionlessparameters(8);
Deff=dimensionlessparameters(9);
K=dimensionlessparameters(10);
alpha=dimensionlessparameters(11);
phis=dimensionlessparameters(12);

dx=1/(N-1);
x=linspace(0,1,N);

T_end=10;
tout=linspace(0,T_end,1e4);

cl0=zeros(N,1);
u10=ones(N,1);
u20=ones(N,1);

disp(['The value of D_eff is ' num2str(Deff)])
disp(['The value of D_s1 is ' num2str(Ds1)])
disp(['The value of D_s2 is ' num2str(Ds2)])
disp(['The value of Q_1 is ' num2str(Q1)])
disp(['The value of Q_2 is ' num2str(Q2)])
disp(['The value of b_et1 is ' num2str(bet1)])
disp(['The value of b_et2 is ' num2str(bet2)])
disp(['The value of K is ' num2str(K)])
disp(['The value of beta is ' num2str(beta)])
disp(['The value of nu is ' num2str(q)])

% Concatenate the initial data
u0=[cl0; u10];
for i=1:N-1
    u0=[u0; u10];
end
for i=1:N
    u0=[u0; u20];
end

% Build the mass matrix
M=build_mass(N,dx,x);

% Solve the problem
options=odeset('Mass',M);
params=[N dx Deff Ds1 Ds2 bet1 bet2 K Q1 Q2 beta phis q];
[t,u]=ode15s(@(t,u) RHS(t,u,params,x),tout,u0,options);

% Pull out the concentration at the exit
c_exit=zeros(length(tout),1);
for i=1:length(tout)
    c_exit(i)=u(i,N);
end

% Do the post-processing to find the EY
extract=trapz(tout,c_exit);
disp('----------------------------')
disp([num2str(1e2*q*beta*extract/phis) '% of the coffee was extracted'])
disp(['The EY was ' num2str(1e2*alpha*q*beta*extract/phis) '%'])
comp_time=toc;
disp(['The computation took ' num2str(comp_time) ' seconds'])
disp('----------------------------')

if plot_flag==1
    plotresults(tout,u,N,x)
end
if checkmass_flag==1
    checkmass(tout,u,x,params)
end