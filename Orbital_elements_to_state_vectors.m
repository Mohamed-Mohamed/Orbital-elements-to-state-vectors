%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg




% Get orbit from orbital element 
close all; clear all; clc;
%% constants
G=6.67428E-11; % gravitational constant
%% inputs
muo=398600; % Gravitational Parameter
M=[5.972*10^24,419600];       % [M1,M2]
R1=[0;0;0];                 % position of M(1)
R2=[3207;5459*0;2714];                 % position of M(2)
I=[0,0,0];              % location of initial axis
V1=[0;0;0];                 % velocity of M(1)
V2=[-6.532;0.7835;6.142];                 % velocity of M(2)
%% get orbital element
[ h1, h, i, Omega, e1, e, w, theta ] = OrbitalElements ( R2,V2,muo );
%% RK4 parameter of orbital element 
order=6;
X0=[h;e;i;Omega;w;theta];
B=[0;0;0;0;0;0;];
sol(1:6,1)=X0;
dt=1000;
t_initial=0;
t_final=1e6;
%% solution of orbital element 
for n=1:length(t_initial:dt:t_final)
    A=[0,0,0,0,0,0;...
      0,0,0,0,0,0;...
      0,0,0,0,0,0;...
      0,0,0,0,0,0;...
      0,0,0,0,0,0;...
      muo^2/sol(1,n)^4*(1+sol(2,n)*cosd(sol(6,n)))^2,0,0,0,0,0];    
    [ XX, t ] = RK4( A,B,sol(1:6,n),dt,n*dt,(n+1)*dt,order );
    sol(1:6,n+1)=XX(1:6,2);
end
for m=1:length(sol(1,:))
    [ r, v ] = OrbitalElements2rvGeo( sol(1,m), muo, sol(2,m), sol(3,m), sol(4,m), sol(5,m), sol(6,m) );
    r_XYZ(1:3,m)=r;
    v_XYZ(1:3,m)=v;
end
%% RK4 parameter of orbital element with oblateness
order1=6;
X01=[h;e;i;Omega;w;theta];
B1=[0;0;0;0;0;0;];
sol1(1:6,1)=X0;
Re=6378;
J2=1.08263e-3;
%% solution of orbital element with oblateness
for n=1:length(t_initial:dt:t_final)
    A1=[0,0,0,0,0,0;...
      0,0,0,0,0,0;...
      0,0,0,0,0,0;...
      -(3/2/sol1(1,n)^8*muo^4*J2*Re^2*(1-sol1(2,n)^2)^(3/2))*cosd(sol1(3,n)),0,0,0,0,0;...
      -(3/2/sol1(1,n)^8*muo^4*J2*Re^2*(1-sol1(2,n)^2)^(3/2))*(5/2*(sind(sol1(3,n)))^2-2),0,0,0,0,0;...
      muo^2/sol1(1,n)^4*(1+sol1(2,n)*cosd(sol1(6,n)))^2,0,0,0,0,0];    
    [ XX1, t1 ] = RK4( A1,B1,sol1(1:6,n),dt,n*dt,(n+1)*dt,order1 );
    sol1(1:6,n+1)=XX1(1:6,2);
end
for m=1:length(sol1(1,:))
    [ r1, v1 ] = OrbitalElements2rvGeo( sol1(1,m), muo, sol1(2,m), sol1(3,m), sol1(4,m), sol1(5,m), sol1(6,m) );
    r_XYZ1(1:3,m)=r1;
    v_XYZ1(1:3,m)=v1;
end
%% RK45 parameter
% X02=[R1;R2;V1;V2];
% B2=[0;0;0;0;0;0;0;0;0;0;0;0];
% sol2(1:12,1)=X02;
% order2=12;
%% RK4
% for n=1:length(t_initial:dt:t_final)
%     b=G*M(2)/(norm(sol2(1:3,n)-sol2(4:6,n)))^3;
%     c=-G*M(1)/(norm(sol2(1:3,n)-sol2(4:6,n)))^3;
%     A2=[0,0,0,0,0,0,1,0,0,0,0,0; ...
%         0,0,0,0,0,0,0,1,0,0,0,0; ...
%         0,0,0,0,0,0,0,0,1,0,0,0; ...
%         0,0,0,0,0,0,0,0,0,1,0,0; ...
%         0,0,0,0,0,0,0,0,0,0,1,0; ...
%         0,0,0,0,0,0,0,0,0,0,0,1;...
%         -b,0,0,b,0,0,0,0,0,0,0,0; ...
%         0,-b,0,0,b,0,0,0,0,0,0,0; ...
%         0,0,-b,0,0,b,0,0,0,0,0,0; ...
%         -c,0,0,c,0,0,0,0,0,0,0,0; ...
%         0,-c,0,0,c,0,0,0,0,0,0,0; ...
%         0,0,-c,0,0,c,0,0,0,0,0,0 ];
%     [ XX2 ] = RK4( A2,B2,sol2(1:12,n),dt,n*dt,(n+1)*dt,order2 );
%     sol2(1:12,n+1)=XX2(1:12,2);
% end
% R1_x=sol2(1,:);
% R1_y=sol2(2,:);
% R1_z=sol2(3,:);
% R2_x=sol2(4,:);
% R2_y=sol2(5,:);
% R2_z=sol2(6,:);
% V1_x=sol2(7,:);
% V1_y=sol2(8,:);
% V1_z=sol2(9,:);
% V2_x=sol2(10,:);
% V2_y=sol2(11,:);
% V2_z=sol2(12,:);
%% plotting
figure(1);
view(3);
hold all;
set(gcf,'color','w');
subplot(1,2,1)
plot3(r_XYZ(1,:),r_XYZ(2,:),r_XYZ(3,:),'linewidth',2);
grid on;
xlabel('X','fontsize',18);
ylabel('Y','fontsize',18);
zlabel('Z','fontsize',18);
title('Solution without oblateness','fontsize',18);
subplot(1,2,2)
plot3(r_XYZ1(1,:),r_XYZ1(2,:),r_XYZ1(3,:),'linewidth',2);
grid on;
xlabel('X','fontsize',18);
ylabel('Y','fontsize',18);
zlabel('Z','fontsize',18);
title('Solution with oblateness','fontsize',18);
% error
figure(2);
set(gcf,'color','w');
% plot3(I(1)*ones(1,length(R1_x))+R2_x-R1_x,I(2)*ones(1,length(R1_y))+R2_y-R1_y,I(3)*ones(1,length(R1_z))+R2_z-R1_z,'g','LineWidth',2);
%% error 
for k=1:length(r_XYZ1(2,:))
    E1(1,k)=norm([r_XYZ(1:3,k)]);
    E2(1,k)=norm([r_XYZ1(1:3,k)]);
end
[ E,Max_e,std_e, mean_e, RMS_e ] = ERROR ( E1(1,:),E2 );
figure(2);
plot(t_initial:dt:t_final+dt,E)
xlim([t_initial, t_final+dt])
xlabel('Time (sec)','fontsize',18);
ylabel('Sol_w_i_t_h_o_u_t _o_b_l-Sol_w_i_t_h _o_b_l','fontsize',18);
title('Error','fontsize',18);
legend(['Max.(Error) = ' num2str(Max_e) ' Std(Error) = ' num2str(std_e) ' mean(Error) = ' num2str(mean_e) ' RMS(Error) = ' num2str(RMS_e) ]);
grid on;