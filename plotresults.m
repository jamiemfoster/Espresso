function plotresults(tout,u,N,x)

c_exit=zeros(length(tout),1);
for i=1:length(tout)
    c_exit(i)=u(i,N);
end

h=figure(1);
hold on
plot(tout,c_exit,'k','linewidth',3);
grid on
box on
xlabel('$t$','Interpreter','latex','FontSize',18)
ylabel('$c_{exit}$','Interpreter','latex','FontSize',18)

figure(2);
colormap(flipud(hot))
shading interp
hold on
cs1=u(:,N+1:N*(N+1));
cs1_dum=cs1(1,:);
cs1_dum=reshape(cs1_dum,[N,N]);
h=surf(x,x,cs1_dum);
set(h,'edgecolor','none')
for j=1:9
    cs1_dum=cs1(j*length(tout)/10,:);
    cs1_dum=reshape(cs1_dum,[N,N]);
    h=surf(x,x,cs1_dum);
    set(h,'edgecolor','none')
end
cs1_dum=cs1(end,:);
cs1_dum=reshape(cs1_dum,[N,N]);
h=surf(x,x,cs1_dum);
set(h,'edgecolor','none')
xlabel('$z$','Interpreter','latex','FontSize',18)
ylabel('$r$','Interpreter','latex','FontSize',18)
zlabel('$c_{s1}$','Interpreter','latex','FontSize',18)
campos([5.5 6.5 4])
set(gca,'xlim',[0 1],'ylim',[0 1],'zlim',[0 1]);

figure(3);
colormap(flipud(hot))
shading interp
hold on
cs2=u(:,N*(N+1)+1:end);
cs2_dum=cs2(1,:);
cs2_dum=reshape(cs2_dum,[N,N]);
h=surf(x,x,cs2_dum);
set(h,'edgecolor','none')
for j=1:9
    cs2_dum=cs2(j*length(tout)/10,:);
    cs2_dum=reshape(cs2_dum,[N,N]);
    h=surf(x,x,cs2_dum);
    set(h,'edgecolor','none')
end
cs2_dum=cs2(end,:);
cs2_dum=reshape(cs2_dum,[N,N]);
h=surf(x,x,cs2_dum);
set(h,'edgecolor','none')
xlabel('$z$','Interpreter','latex','FontSize',18)
ylabel('$r$','Interpreter','latex','FontSize',18)
zlabel('$c_{s2}$','Interpreter','latex','FontSize',18)
campos([5.5 6.5 4])
set(gca,'xlim',[0 1],'ylim',[0 1],'zlim',[0 1]);

figure(4)
hold on
plot(x,u(1,1:N),'.k','MarkerSize',10)
plot(x,u(1,1:N),'k')
for i=1:9
    plot(x,u(i*length(tout)/10,1:N),'.','MarkerSize',10,'color',[0 0 i/10])
    plot(x,u(i*length(tout)/10,1:N),'color',[0 0 (i-1)/10])
end
plot(x,u(end,1:N),'.b','MarkerSize',10)
plot(x,u(end,1:N),'b')
grid on
box on
xlabel('$z$','Interpreter','latex','FontSize',18)
ylabel('$c_l$','Interpreter','latex','FontSize',18)