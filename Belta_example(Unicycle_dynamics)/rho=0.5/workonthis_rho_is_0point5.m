clear
clc
close all
n=3;
P=2;
M=5;
Weights{1}=2*rand(M,n+1)-1;
Weights{2}=2*rand(P,M)-1;
Biases{1}=2*rand(M,1)-1;
Biases{2}=2*rand(P,1)-1;
q_mean=1;
epsi=0.01;


discr=100;
Beta1=0.95;
Beta2=0.95;
eps=0.0000001;
tol=10^(-10);


T=20;
NN=40000;
Batch_num=1000;
len=M*(n+1+1)+P*(M+1);
param=cell(1,NN);
lam_param=cell(1,NN);
beta_param=cell(1,NN);

load('candidate2_initial');
param{1}=Param;
% param{1}=2*rand(1,len)-1;
beta_param{1}=ones(1,22);
lam_param{1}=2*rand-1;

J1=cell(1,NN+1);
J2=cell(1,NN+1);
Bs=cell(1,NN+1);
alpha=0.1;
i=0;

s0=[0.5+0.1*rand(2,1); pi/2]; 
q=q_mean+epsi*(2*rand-1);

J1{1}=Perf_obj(s0,param{1},q,M,P,T);
[J2{1}, Bs{1}]=STL_obj(s0,param{1}, lam_param{1}, beta_param{1}, q,M,P,T);
rho=0.5;
rng(0)

KK=60;



% s0s=[0.5+rand(2,KK); pi/2+0.1*pi*(2*rand(1,KK)-1)];
tic
for iter=1:NN
    s0s=[0.5+rand(2,KK); pi/2+0.1*pi*(2*rand(1,KK)-1)];
    q=q_mean+epsi*(2*rand-1);
    
    if mod(iter, Batch_num)==1
        Vd=cell(1,Batch_num);Vd{1}=zeros(1,len+23);
        Sd=cell(1,Batch_num);Sd{1}=zeros(1,len+23);
        Vd1=cell(1,Batch_num);Vd1{1}=zeros(1,len+23);
        Sd1=cell(1,Batch_num);Sd1{1}=zeros(1,len+23);
        Vd2=cell(1,Batch_num);Vd2{1}=zeros(1,len+23);
        Sd2=cell(1,Batch_num);Sd2{1}=zeros(1,len+23);
        index=0;
    end
    index=index+1;
    
    
    parfor initNum=1:KK
        s0 = s0s(:,initNum);
        Ggrad1(initNum,:)=[Back_prop_Perf(s0, param{iter}, q,M, P, T) 0 zeros(1,22)]; %#ok<*PFBNS>
        Ggrad2(initNum,:)=Back_prop_STL(s0, param{iter}, lam_param{iter}, beta_param{iter}, q,M, P, T);
    end
    if ~isnan(sum(sum(Ggrad1)))
        Norms1=norms(Ggrad1',2);
        i1=find(Norms1==max(Norms1)); i1=i1(1);
        grad1 = Ggrad1(i1,:);
    else 
        grad1=Grad(iter-1,:);
    end
    if ~isnan(sum(sum(Ggrad2)))
        Norms2=norms(Ggrad2',2);
        i2=find(Norms2==max(Norms2)); i2=i2(1);
        grad2 = Ggrad2(i2,:);

    else
        grad2=Grad(iter-1,:);
    end
    
    
    fprintf('Iteration %d: G1 = %s. G2 = %s.\n\n', iter, mat2str(grad1,1), mat2str(grad2,1));
    

    Vd1{index+1} = Beta1*Vd{index}+(1-Beta1)*grad1   ;
    Sd1{index+1} = Beta2*Sd{index}+(1-Beta2)*grad1.^2;
    Vd2{index+1} = Beta1*Vd{index}+(1-Beta1)*grad2   ;
    Sd2{index+1} = Beta2*Sd{index}+(1-Beta2)*grad2.^2;
    
    
    S1=[param{iter} , lam_param{iter}, beta_param{iter}] + alpha*Vd1{index+1}./sqrt(Sd1{index+1}+eps);
    S1p=[param{iter} , lam_param{iter}, beta_param{iter}] + (1/discr) * alpha*Vd1{index+1}./sqrt(Sd1{index+1}+eps);
    S2=[param{iter} , lam_param{iter}, beta_param{iter}] + alpha*Vd2{index+1}./sqrt(Sd2{index+1}+eps);
    
    
    
    paramnext1=S1(1,1:len);
    paramnext1p=S1p(1,1:len);
    paramnext2=S2(1,1:len);
    lam_paramnext1=S1(1,len+1);
    lam_paramnext1p=S1p(1,len+1);
    lam_paramnext2=S2(1,len+1);
    beta_paramnext1=S1(1,len+2:end);
    beta_paramnext1p=S1p(1,len+2:end);
    beta_paramnext2=S2(1,len+2:end);
    
    s0=[0.5+rand(2,1); pi/2+0.1*pi*(2*rand-1)];
    q=q_mean+epsi*(2*rand-1);
    
    J2_discuss= STL_obj(s0,param{iter}, lam_param{iter}, beta_param{iter}, q,M,P,T);
    J1next1=Perf_obj(s0,paramnext1,q,M,P,T);
    J1next1p=Perf_obj(s0,paramnext1p,q,M,P,T);
    J1next2=Perf_obj(s0,paramnext2,q,M,P,T);
    [J2next1, Bsnext1]=STL_obj(s0,paramnext1, lam_paramnext1, beta_paramnext1, q,M,P,T);
    [J2next1p, Bsnext1p]=STL_obj(s0,paramnext1p, lam_paramnext1p, beta_paramnext1p, q,M,P,T);
    [J2next2, Bsnext2]=STL_obj(s0,paramnext2, lam_paramnext2, beta_paramnext2, q,M,P,T);
    if J2_discuss<rho
        if J2next1>J2_discuss
            param{iter+1}=paramnext1;
            lam_param{iter+1}=lam_paramnext1;
            beta_param{iter+1}=beta_paramnext1;
            J1{iter+1}=J1next1;
            J2{iter+1}=J2next1; Bs{iter+1}=Bsnext1;
            Vd{index+1}=Vd1{index+1}; Sd{index+1}=Sd1{index+1};
            Grad(iter,:)  =grad1;
        else
            param{iter+1}=paramnext2;
            lam_param{iter+1}=lam_paramnext2;
            beta_param{iter+1}=beta_paramnext2;
            J1{iter+1}=J1next2;
            J2{iter+1}=J2next2; Bs{iter+1}=Bsnext2;
            Vd{index+1}=Vd2{index+1}; Sd{index+1}=Sd2{index+1};
            Grad(iter,:)  =grad2;
        end
    else
        param{iter+1}=paramnext1p;
        lam_param{iter+1}=lam_paramnext1p;
        beta_param{iter+1}=beta_paramnext1p;
        J1{iter+1}=J1next1p;
        J2{iter+1}=J2next1p; Bs{iter+1}=Bsnext1p;
        Vd{index+1}=(1/discr)*Vd1{index+1}; Sd{index+1}=Sd1{index+1};
        Grad(iter,:)  =(1/discr)*grad1;
    end
end
Run_time=toc;



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
plot(x4(1,:),x4(2,:), '-b');
hold on
end
% the following line produces trajectories for the trained controller with
% initial condition and model uncertainty
Ma=[min(J1s(floor(0.9*iter):end)),max(J1s(floor(0.9*iter):end))];
plotTraj([0.6 1.4],[0.6 1.4],[pi/2-pi/10 , pi/2+pi/10],[0.99, 1.01], [5, 5, 4, 3], M, Ma, P, Param,  T, 4, 5);



Param=param{1};

for fignum=6:7
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
plot(x4(1,:),x4(2,:), '-b');
hold on
end
% the following line produces trajectories for the trained controller with
% initial condition and model uncertainty
Ma=[0 10000];
plotTraj([0.6 1.4],[0.6 1.4],[pi/2-pi/10 , pi/2+pi/10],[0.99, 1.01], [2, 2, 2, 2], M, Ma, P, Param,  T, 6, 7);

