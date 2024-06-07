clear
clc
close all

% load('result_rho_is_0point5.mat')
% load('result_rho_is_0point5_extendedto_10000.mat')
load('result_rho_is_0point5_extendedto_40000.mat')
s=zeros(iter-1,len+1+22);
Bbs=zeros(T+1,iter-1);
J1s=zeros(1,iter-1);
J2s=zeros(1,iter-1);
STL=zeros(1,NN);
for j=1:iter-1
    s(j,:)=[param{j} lam_param{j} beta_param{j}];
    Bbs(:,j)=Bs{j}';
    J1s(1,j)=J1{j};
    J2s(1,j)=J2{j};
end


figure(1)
for j=1:len+1+22
    plot(s(:,j))
    hold on
end


figure(2)
plot(J1s,'blue')
hold on
plot(J2s,'green')

figure(3)
for j=1:T+1
    plot(Bbs(j,:))
    hold on
end
plot(J2s,'-green')


Param=param{iter};
Lam=lam_param{iter};
Bet=beta_param{iter};


[JP, J_STL_h, J_STL_s] = expectations_Unicycle([0.6 1.4],[0.6 1.4],[pi/2-pi/10 , pi/2+pi/10],[0.99, 1.01], 1, M, P, Param, Lam, Bet, T);

plotTraj_rand([0.6 1.4],[0.6 1.4],[pi/2-pi/10 , pi/2+pi/10],[0.99, 1.01], 550, M, P, Param,  T, 4, 5);
for fignum=4:5
figure(fignum)
hold on
NNN=10000;
x1=zeros(2,NNN);
x2=zeros(2,NNN);
x3=zeros(2,NNN);
x4=zeros(2,NNN);
for j=1:NNN
    thet=2*pi*j/NNN;
    r=[cos(thet);sin(thet)];
    x1(:,j)=[5;5]+sqrt(1.5)*r;
    x2(:,j)=[2;8]+sqrt(1.5)*r;
    x3(:,j)=[8;2]+sqrt(1.5)*r;
    x4(:,j)=[8;8]+sqrt(1.5)*r;
end

plot(x1(1,:),x1(2,:), '-r');
hold on
plot(x2(1,:),x2(2,:), '-b');
hold on
plot(x3(1,:),x3(2,:), '-b');
hold on
end
plot(8,8, '*', 'color', 'black','Markersize', 20, 'Linewidth', 2 )
% the following line produces trajectories for the trained controller with
% initial condition and model uncertainty
% plotTraj([0.6 1.4],[0.6 1.4],[pi/2-pi/10 , pi/2+pi/10],[0.99, 1.01], [5, 5, 4, 3], M, Ma, P, Param,  T, 4, 5);



figure(6)
Param=param{1};
for ii=1:100
    q= 0.99+rand*(0.02);
    s0 = [0.6;0.6;2*pi/5]+rand(3,1).*[0.8;  0.8;  pi/5];
    a0 = Actor(s0,0,M,P,Param);
    xx = zeros(3,T);
    xx(:,1)=Dynamics(s0,a0,q);
    for ij=1:T
        a(:,ij)=Actor(xx(:,ij),ij,M,P,Param);
        xx(:,ij+1)=Dynamics(xx(:,ij),a(:,ij),q);
    end
    plot([s0(1,1) xx(1,:)],[s0(2,1) xx(2,:)], 'Linewidth', 1.5);
    hold on
end
xlim([-5 15]);
ylim([-5 15]);
hold on
NNN=10000;
x1=zeros(2,NNN);
x2=zeros(2,NNN);
x3=zeros(2,NNN);
x4=zeros(2,NNN);
for j=1:NNN
    thet=2*pi*j/NNN;
    r=[cos(thet);sin(thet)];
    x1(:,j)=[5;5]+sqrt(1.5)*r;
    x2(:,j)=[2;8]+sqrt(1.5)*r;
    x3(:,j)=[8;2]+sqrt(1.5)*r;
end

plot(x1(1,:),x1(2,:), '-r');
hold on
plot(x2(1,:),x2(2,:), '-b');
hold on
plot(x3(1,:),x3(2,:), '-b');
hold on
plot(8,8, '*', 'color', 'black','Markersize', 20, 'Linewidth', 2 )