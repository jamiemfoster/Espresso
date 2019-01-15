function checkmass(tout,u,x,params)

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

% Check that solubles are conserved in the liquid
Ml1=zeros(length(tout),1);
for i=1:length(tout)
    Ml1(i)=(1-phis)*trapz(x,u(i,1:N));
end

Ml2=zeros(length(tout),1);
c_exit=zeros(length(tout),1);
R1=zeros(length(tout),1);
R2=zeros(length(tout),1);
for i=1:length(tout)
    c_exit(i)=u(i,N);
    csurf1=zeros(N,1);
    for j=1:N
        csurf1(j)=u(i,(j+1)*N);
    end
    csurf2=zeros(N,1);
    for j=1:N
        csurf2(j)=u(i,N*N+(j+1)*N);
    end
    R1(i)=trapz(x,K*(1-u(i,1:N))'.*csurf1.*(csurf1-beta.*u(i,1:N)'));
    R2(i)=trapz(x,K*(1-u(i,1:N))'.*csurf2.*(csurf2-beta.*u(i,1:N)'));
end
for i=1:length(tout)
    Ml2(i)=Ml1(1)+trapz(tout(1:i),-q*c_exit(1:i)+bet1.*R1(1:i)+bet2.*R2(1:i),1);
end


figure(5)
subplot(2,2,1)
hold on
plot(tout,Ml1,'ok')
plot(tout,Ml2,'.r')
grid on
box on
xlabel('t')
ylabel('Ml')

% Check that solubles are conserved in the fines
V=eye(N,N);
V(1,1)=4/3*pi*(dx/2)^3;
for i=2:N-1
    V_dum=4*pi*(x(i)^2*dx+dx^3/12);
    V(i,i)=V_dum;
end
V(N,N)=4/3*pi-4/3*pi*(1-dx/2)^3;

Ms1_dum=zeros(length(tout),N);
for i=1:length(tout)
    for j=1:N
        Ms1_dum(i,j)=sum(V*u(i,j*N+1:(j+1)*N)')/4/pi;
    end
end

Ms1_1=zeros(length(tout),1);
for i=1:length(tout)
    Ms1_1(i)=bet1*trapz(x,Ms1_dum(i,:))/beta/Q1;
end

Ms1_2=zeros(length(tout),1);
for i=1:length(tout)
    Ms1_2(i)=Ms1_1(1)-bet1*trapz(tout(1:i),R1(1:i),1);
end

Ms1_init=bet1/beta/Q1/3;
figure(5)
subplot(2,2,2)
hold on
plot(tout,Ms1_1,'ok')
plot(tout,Ms1_2,'.r')
plot(0,Ms1_init,'.g','MarkerSize',15)
grid on
box on
xlabel('t')
ylabel('Ms1')

% Check that solubles are conserved in the boulders
Ms2_dum=zeros(length(tout),N);
for i=1:length(tout)
    for j=1:N
        Ms2_dum(i,j)=sum(V*u(i,N*N+j*N+1:N*N+(j+1)*N)')/4/pi;
    end
end

Ms2_1=zeros(length(tout),1);
for i=1:length(tout)
    Ms2_1(i)=bet2*trapz(x,Ms2_dum(i,:))/beta/Q2; 
end

Ms2_2=zeros(length(tout),1);
for i=1:length(tout)
    Ms2_2(i)=Ms2_1(1)-bet2*trapz(tout(1:i),R2(1:i),1);
end

Ms2_init=bet2/beta/Q2/3;
figure(5)
subplot(2,2,3)
hold on
plot(tout,Ms2_1,'ok')
plot(tout,Ms2_2,'.r')
plot(0,Ms2_init,'.g','MarkerSize',15)
grid on
box on
xlabel('t')
ylabel('Ms2')

% Check that solubles are conserved globally
M_tot1=Ml1+Ms1_1+Ms2_1;
M_tot2=zeros(length(tout),1);
for i=1:length(tout)
    M_tot2(i)=M_tot1(1)+trapz(tout(1:i),-q.*c_exit(1:i),1);
end
figure(5)
subplot(2,2,4)
hold on
plot(tout,M_tot1,'ok')
plot(tout,M_tot2,'.r')
plot(0,Ms1_init+Ms2_init,'.g','MarkerSize',15)
grid on
box on
xlabel('t')
ylabel('Mtot')