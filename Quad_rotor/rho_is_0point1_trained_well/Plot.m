clear
clc
close all
% load('Result_rho_is_0point1.mat')
load('Result_rho_is_0point1_extendedto_47000.mat')
s=zeros(iter-1,len+lenSTL);
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
for j=1:len
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


[JP, J_STL_h, J_STL_s] = expectations_Quadrotor([0,1],[0,2*pi],[-pi/2, pi/2], [-0.01, 0.01], As, B, 10, M, Param, Lam, Bet, T);


for fignum=4:5
    figure(fignum)
    
    [X,Y,Z]=sphere;
    r  = 0.5 *FF;
    X2 = X * r;
    Y2 = Y * r;
    Z2 = Z * r;
    
    surface=surf(X2+FF*2,Y2+FF*2,Z2,'FaceAlpha',0.2);
    surface.EdgeColor = 'none';
    
    
    
    hold on
    
    
    [X,Y,Z]=sphere;
    r=sqrt(1.5)*FF;
    X2 = X * r;
    Y2 = Y * r;
    Z2 = Z * r;
    
    surface=surf(X2+FF*5,Y2+FF*5,Z2,'FaceAlpha',0.2,'Facecolor',[0.8500 0.3250 0.0980], 'EdgeColor', [1,0,0]);
%     surface.EdgeColor = 'none';
    
    
    hold on
    [X,Y,Z]=sphere;
    r=sqrt(1.5)*FF;
    X2 = X * r;
    Y2 = Y * r;
    Z2 = Z * r;
    
    surface=surf(X2+FF*2,Y2+FF*8,Z2,'FaceAlpha',0.2,'Facecolor',[0,1,1], 'EdgeColor', [0,0,1]);
%     surface.EdgeColor = 'none';
    
    
    
    hold on
    [X,Y,Z]=sphere;
    r=sqrt(1.5)*FF;
    X2 = X * r;
    Y2 = Y * r;
    Z2 = Z * r;
    
    surface=surf(X2+FF*8,Y2+FF*2,Z2,'FaceAlpha',0.2, 'Facecolor',[0,1,1], 'EdgeColor', [0,0,1]);
%     surface.EdgeColor = 'none';
        
    
    axis equal
    
end


hold on



rng(0)

rrs     = [0,1];
pphis   = [0,2*pi];
tthetas = [-pi/2, pi/2];
deltas  = [-0.01, 0.01];


l=20.5;
u=28;%27.5;;
lowest=0.2;

for ii=1:750
    rr=rrs(1)+rand*(rrs(2)-rrs(1));
    pphi=pphis(1)+rand*(pphis(2)-pphis(1));
    ttheta=tthetas(1)+rand*(tthetas(2)-tthetas(1));
    delta_model= deltas(1)+rand*(deltas(2)-deltas(1));
    initxx = 2+ 0.5*rr*cos(ttheta)*cos(pphi);
    inityy = 2+ 0.5*rr*cos(ttheta)*sin(pphi);
    initzz =    0.5*rr*sin(ttheta);
    s0 = FF* [initxx;inityy;initzz;zeros(3,1)];
    A=As+delta_model*eye(6);
    a0 = Actor(s0,0,M,P,Param);
    aa=a0;
    uu=[tan(0.1*tanh(aa(1,1)/10)); tan(0.1*tanh(aa(2,1)/10));   2*tanh(aa(3,1)/10) ];
    xx = zeros(6,T);
    xx(:,1)=A*s0+B*uu;
    for ij=1:T
        a(:,ij)=Actor(xx(:,ij),ij,M,P,Param);
        aa=a(:,ij);
        uu=[tan(0.1*tanh(aa(1,1)/10)); tan(0.1*tanh(aa(2,1)/10));   2*tanh(aa(3,1)/10) ];
        xx(:,ij+1)=A*xx(:,ij)+B*uu;
    end
    
    rob = robustness_Quadrotor(s0,Param, M, T, A, B);
    
    if (rob < 0)
        figure(4)
%         title('Violating Trajectories')
        color = '-r.' ;
        plot3([s0(1,1) xx(1,:)],[s0(2,1) xx(2,:)], color, 'Linewidth', 0.75);
        
    else
        figure(5)
%         title('Satisfactory Trajectories')
        reward=Perf_obj(s0, Param, M, T, A, B)
        color=[  reward/(l-u)-u/(l-u),   (1-lowest)*reward/(l-u)+(lowest*l-u)/(l-u), reward/(l-u)-u/(l-u)];
        plot3([s0(1,1) xx(1,:)],[s0(2,1) xx(2,:)],[s0(3,1) xx(3,:)], 'color',color, 'Linewidth', 1);
        
    end
    hold on
end
plot3(FF*8, FF*8, -FF*3, '*', 'color', 'black','Markersize', 60, 'Linewidth', 3)



Param=param{1};
T=10;
figure(6)
[X,Y,Z]=sphere;
r  = 0.5 *FF;
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;
surface=surf(X2+FF*2,Y2+FF*2,Z2,'FaceAlpha',0.2);
surface.EdgeColor = 'none';
hold on
[X,Y,Z]=sphere;
r=sqrt(1.5)*FF;
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;
surface=surf(X2+FF*5,Y2+FF*5,Z2,'FaceAlpha',0.2,'Facecolor',[0.8500 0.3250 0.0980], 'EdgeColor', [1,0,0]);
%     surface.EdgeColor = 'none';
hold on
[X,Y,Z]=sphere;
r=sqrt(1.5)*FF;
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;
surface=surf(X2+FF*2,Y2+FF*8,Z2,'FaceAlpha',0.2,'Facecolor',[0,1,1], 'EdgeColor', [0,0,1]);
%     surface.EdgeColor = 'none';
hold on
[X,Y,Z]=sphere;
r=sqrt(1.5)*FF;
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;
surface=surf(X2+FF*8,Y2+FF*2,Z2,'FaceAlpha',0.2, 'Facecolor',[0,1,1], 'EdgeColor', [0,0,1]);
%     surface.EdgeColor = 'none';
axis equal
hold on

rng(0)

rrs     = [0,1];
pphis   = [0,2*pi];
tthetas = [-pi/2, pi/2];
deltas  = [-0.01, 0.01];
for ii=1:100
    rr=rrs(1)+rand*(rrs(2)-rrs(1));
    pphi=pphis(1)+rand*(pphis(2)-pphis(1));
    ttheta=tthetas(1)+rand*(tthetas(2)-tthetas(1));
    delta_model= deltas(1)+rand*(deltas(2)-deltas(1));
    initxx = 2+ 0.5*rr*cos(ttheta)*cos(pphi);
    inityy = 2+ 0.5*rr*cos(ttheta)*sin(pphi);
    initzz =    0.5*rr*sin(ttheta);
    s0 = FF* [initxx;inityy;initzz;zeros(3,1)];
    A=As+delta_model*eye(6);
    a0 = Actor(s0,0,M,P,Param);
    aa=a0;
    uu=[tan(0.1*tanh(aa(1,1)/10)); tan(0.1*tanh(aa(2,1)/10));   2*tanh(aa(3,1)/10) ];
    xx = zeros(6,T);
    xx(:,1)=A*s0+B*uu;
    for ij=1:T
        a(:,ij)=Actor(xx(:,ij),ij,M,P,Param);
        aa=a(:,ij);
        uu=[tan(0.1*tanh(aa(1,1)/10)); tan(0.1*tanh(aa(2,1)/10));   2*tanh(aa(3,1)/10) ];
        xx(:,ij+1)=A*xx(:,ij)+B*uu;
    end
        plot3([s0(1,1) xx(1,:)],[s0(2,1) xx(2,:)],[s0(3,1) xx(3,:)], 'Linewidth', 1);
       
    hold on
end
plot3(FF*8, FF*8, -FF*3, '*', 'color', 'black','Markersize', 20, 'Linewidth', 3)


figure(7)
V_sum=sum(exp(s(:,len+1+10+1:len+1+10+10)),2);
for j=len+1+10+1:len+1+10+10
    V=exp(s(:,j));
    V=V./V_sum;
    plot(V,'Linewidth', 2)
    hold on
end    
V_sum=sum(exp(s(:,len+1+21:len+1+22)),2);
V=exp(s(:,len+1+21));
V=V./V_sum;
plot(V,'--blue','Linewidth', 2)
hold on
V=exp(s(:,len+1+22));
V=V./V_sum;
plot(V,'--red','Linewidth', 2)





% rrs     = linspace(0,1,10);
% pphi    = linspace(0,2*pi,10);
% ttheta  = linspace(-pi/2,pi/2,8);
% deltas  = linspace(-0.01,0.01,6);
% 
% 
% for ii=1:length(rrs) 
%     for jj=1:length(pphi)
%         for kk=1:length(ttheta)
%             for ll=1:length(deltas)
%                 A=As+deltas(ll)*eye(6);
%                 initxx = 2+ 0.5*rrs(ii)*cos(ttheta(kk))*cos(pphi(jj));
%                 inityy = 2+ 0.5*rrs(ii)*cos(ttheta(kk))*sin(pphi(jj));
%                 initzz =    0.5*rrs(ii)*sin(ttheta(kk));
%                 s0 = FF* [initxx;inityy;initzz;zeros(3,1)];
%                 a0 = Actor(s0,0,M,P,Param);
%                 aa=a0;
%                 uu=[tan(0.1*tanh(aa(1,1)/10)); tan(0.1*tanh(aa(2,1)/10));   2*tanh(aa(3,1)/10) ];
%                 xx = zeros(6,T);
%                 xx(:,1)=A*s0+B*uu;
%                 for ij=1:T 
%                     a(:,ij)=Actor(xx(:,ij),ij,M,P,Param);
%                     aa=a(:,ij);
%                     uu=[tan(0.1*tanh(aa(1,1)/10)); tan(0.1*tanh(aa(2,1)/10));   2*tanh(aa(3,1)/10) ];
%                     xx(:,ij+1)=A*xx(:,ij)+B*uu;
%                 end      
%             
%                 rob = robustness_Quadrotor(s0,Param, M, T, A, B);
%             
%                 if (rob < 0)
%                     figure(4)
%                     title('Violating Trajectories')
%                     color = '-r.' ;
%                     plot3([s0(1,1) xx(1,:)],[s0(2,1) xx(2,:)], color, 'Linewidth', 0.75);
%                 
%                 else
%                     figure(5)
%                     title('Satisfactory Trajectories')
%                     
%                     reward=Perf_obj(s0, Param, M, T, A, B);
%                     plot3([s0(1,1) xx(1,:)],[s0(2,1) xx(2,:)],[s0(3,1) xx(3,:)], 'color',0.9*(slope*reward+intersect)*[0,1,0], 'Linewidth', 1);
%                     
%                 end
%             
%             end
%         end
%     end
% end
% 
% hold off
% 
