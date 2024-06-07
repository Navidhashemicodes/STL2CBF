clear all
clc
close all
n=6;
P=3;
M=15;
T=35;
delta_t=0.05;
A_c=[zeros(3,3) eye(3);zeros(3,6)];
g=9.81;
B_c=[zeros(3,3); g 0 0; 0 -g 0; 0 0 1];
%%% x(k+1) =A*x(k)+B*u(k) %%%%%%%%
A=expm(A_c*delta_t);
B= G_maker(delta_t,A_c,B_c);

Param=2*rand(1, M*(n+1)+M+P*M+P)-1;

%
s0=[0; 0; 0; zeros(3,1)];
a0 = Actor(s0,0,M,P,Param);
aa=a0;
uu=[tan(0.1*tanh(aa(1,1)/50)); tan(0.1*tanh(aa(2,1)/50));   2*tanh(aa(3,1)/50) ];
xx = zeros(6,T);
% xode=zeros(6,T);
xx(:,1)=A*s0+B*uu;
% [~,in_out] =  ode45(@(t,x)eq1f(x,aa),[0 delta_t],s0);
% xode(:,1)=in_out(end,:)';
for i=1:T
%     a(:,i)=Actor(xx(:,i),i,M,P,Param);
%     aa=a(:,i);
%     uu=[tan(0.1*tanh(aa(1,1)/50)); tan(0.1*tanh(aa(2,1)/50));   2*tanh(aa(3,1)/50) ];
    uu=[0.1;0.1;2];
    xx(:,i+1)=A*xx(:,i)+B*uu;
%     [~,in_out] =  ode45(@(t,x)eq1f(x,aa),[0 delta_t],xode(:,i));
%     xode(:,i+1)=in_out(end,:)';
end
figure
plot(xx(1,:))
hold on
% plot(xode(1,:))


figure
plot(xx(2,:))
hold on
% plot(xode(2,:))


figure
plot(xx(3,:))
hold on
% plot(xode(3,:))


% figure
% plot(xx(4,:))
% hold on
% plot(xode(4,:))
% 
% figure
% plot(xx(5,:))
% hold on
% plot(xode(5,:))
% 
% 
% figure
% plot(xx(6,:))
% hold on
% plot(xode(6,:))



function a=Actor(s,k, M,P, Param)
n=size(s,1);
index=0;
Weights{1}=reshape(Param(1,index+1:index+M*(n+1)), [n+1, M])';  
index=index+M*(n+1);
Biases{1}= Param(1,index+1:index+M)';
index=index+M;
Weights{2}=reshape(Param(1,index+1:index+M*P), [M, P])';
index=index+M*P;
Biases{2}=Param(1,index+1:index+P)';
Sab=Weights{1}*[s;k]+Biases{1};
Sa=tanh(Sab);
a=Weights{2}*Sa+Biases{2};
end