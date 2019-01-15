function M=build_mass(N,dr,r)

M1=eye(N,N);
M2=eye(N,N);
M1(1,1)=4/3*pi*(dr/2)^3;
for i=2:N-1
    V_dum=4*pi*(r(i)^2*dr+dr^3/12);
    M1(i,i)=V_dum;
end
M1(N,N)=4/3*pi-4/3*pi*(1-dr/2)^3;

M2(1,1)=3/4;
M2(1,2)=1/8;
M2(2,1)=1/4;
M2(2,2)=6/8;
M2(2,3)=1/8;
for i=3:N-2
    M2(i,i)=6/8;
    M2(i,i-1)=1/8;
    M2(i,i+1)=1/8;
end
M2(N-1,N-2)=1/8;
M2(N-1,N-1)=6/8;
M2(N-1,N)=1/4;
M2(N,N-1)=1/8;
M2(N,N)=3/4;
M_sub=M2*M1;

% M_sub for each spherical diffusion equation and ones for
% those related to macroscopic transport
M=eye((2*N+1)*N,(2*N+1)*N);
for i=1:2*N
    M(i*N+1:i*N+N,i*N+1:i*N+N)=M_sub;
end

% For the algebraic boundary conditions
M(1,1)=0;
M(N,N)=0;