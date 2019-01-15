function du=RHS(t,u,params,x)
	
N=      params(1);
dx=     params(2);
Deff=   params(3);
Ds1=    params(4);
Ds2=    params(5);
bet1=   params(6);
bet2=   params(7);
K=      params(8);
Q1=     params(9);
Q2=     params(10);
beta=   params(11);
phis=   params(12);
q=      params(13);

du=zeros((2*N+1)*N,1);

% Extract the surface concentrations
csurf1=zeros(N,1);
for i=1:N
    csurf1(i)=u((i+1)*N);
end
csurf2=zeros(N,1);
for i=1:N
    csurf2(i)=u(N*N+(i+1)*N);
end

% Compute the reaction rate
G1=K*(1-u(1:N)).*csurf1.*(csurf1-beta.*u(1:N));
G2=K*(1-u(1:N)).*csurf2.*(csurf2-beta.*u(1:N));

% Eqn for c_l (concentration in the liquid)
gamma=1/(1-phis);

du(1)=-Deff/dx*(-3*u(1)/2+2*u(2)-u(3)/2)+q.*u(1);
for i=2:N-1
    du(i)=  -q.*gamma.*(u(i+1)-u(i-1))/2/dx ...
            +gamma.*Deff*(u(i+1)-2*u(i)+u(i-1))/dx^2 ...
            +gamma.*bet1*G1(i) ...
            +gamma.*bet2*G2(i);
end
du(N)=u(N-2)/2 - 2*u(N-1) + 3*u(N)/2;

% Eqns for cs1
% j is the index for the particle number (with the j=1 particle at z=0)
% i is the index for the grid points in r
for j=1:N
    du(j*N+1)=  4*pi*(x(1)/2+x(2)/2)^2 * Ds1 * (u(j*N+2)-u(j*N+1)) / dx;    
    for i=2:N-1
        du(j*N+i)=  4*pi*(x(i+1)/2+x(i)/2)^2 * Ds1 * (u(j*N+i+1)-u(j*N+i)) / dx ...
                   -4*pi*(x(i)/2+x(i-1)/2)^2 * Ds1 * (u(j*N+i)-u(j*N+i-1)) / dx;
    end
    du(j*N+N)=  4*pi*(x(N-1)/2+x(N)/2)^2 * Ds1 * (u(j*N+N-1)-u(j*N+N)) / dx - 4*pi*beta*Q1*G1(j);
end

% Eqns for cs2
% j is the index for the particle number (with the j=1 particle at z=0)
% i is the index for the grid points in r
for j=N+1:2*N
    du(j*N+1)=  4*pi*(x(1)/2+x(2)/2)^2 * Ds2 * (u(j*N+2)-u(j*N+1)) / dx;    
    for i=2:N-1
        du(j*N+i)=  4*pi*(x(i+1)/2+x(i)/2)^2 * Ds2 * (u(j*N+i+1)-u(j*N+i)) / dx ...
                   -4*pi*(x(i)/2+x(i-1)/2)^2 * Ds2 * (u(j*N+i)-u(j*N+i-1)) / dx;
    end
    du(j*N+N)=  4*pi*(x(N-1)/2+x(N)/2)^2 * Ds2 * (u(j*N+N-1)-u(j*N+N)) / dx - 4*pi*beta*Q2*G2(j-N);
end