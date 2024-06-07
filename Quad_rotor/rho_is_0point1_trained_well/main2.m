% clear all
% clc
% close all
% n=6;
% P=3;
% M=10;
% delta_t=0.05;
% A_c=[zeros(3,3) eye(3);zeros(3,6)];
% g=9.81;
% B_c=[zeros(3,3); g 0 0; 0 -g 0; 0 0 1];
% %%% x(k+1) =A*x(k)+B*u(k) %%%%%%%%
% As=expm(A_c*delta_t);
% B= G_maker(delta_t,A_c,B_c);
% 
% epsi=0.01;
% 
% Weights{1}=2*rand(M,n+1)-1;
% Weights{2}=2*rand(P,M)-1;
% Biases{1}=2*rand(M,1)-1;
% Biases{2}=2*rand(P,1)-1;
% 
% 
% discr=50000;
% NN=10000;
% lambda=rand;
% T=20;
% Beta1_Perf=0.99;
% Beta2_Perf=0.99;
% Beta1_STL=0.5;
% Beta2_STL=0.5;
% eps=0.0000001;
% 
% 
% alpha=0.01;
% Batch_num=400;
% len=M*(n+1+1)+P*(M+1);
% lenSTL=23;
% param=cell(1,NN);
% lam_param=cell(1,NN);
% beta_param=cell(1,NN);
% 
% load('init.mat')
% 
% % param{1}=[reshape(Weights{1}',[1, M*(n+1)]) , Biases{1}', reshape(Weights{2}',[1, P*M]) , Biases{2}'];
% param{1}=the_Param;
% lam_param{1}=lambda;
% % beta_param{1}=ones(1,lenSTL-1);
% beta_param{1}=the_beta;
% J1=cell(1,NN+1);
% J2=cell(1,NN+1);
% Bs=cell(1,NN+1);
% 
% FF=0.1*0.5*0.25;
% 
% Rand3=2*rand(3,1)-1;
% s0 = FF*[  [2;2;0]+0.5*Rand3/norm(Rand3)  ; zeros(3,1)   ];
% A=As+(epsi*(2*rand-1))*eye(6);
% 
% J1{1}=Perf_obj(s0,param{1},M, T,A,B);
% [J2{1}, Bs{1}]=STL_obj(s0,param{1}, lam_param{1}, beta_param{1}, M, T,A,B);
% rho=0.1;
% rng(0)
J1(NN:6*NN)=cell(1,5*NN+1);
J2(NN:6*NN)=cell(1,5*NN+1);
param(NN:6*NN)=cell(1,5*NN+1);
beta_param(NN:6*NN)=cell(1,5*NN+1);
lam_param(NN:6*NN)=cell(1,5*NN+1);
KK=20;

% s0s= FF  *  [ [2;2;0]+(rand(3,KK)-0.5)   ;    zeros(3,KK)   ];
tic
for iter=NN-1:6*NN
    
    A=As+(epsi*(2*rand-1))*eye(6);
    Rand3=2*rand(3,KK)-1;
    s0s= FF  *  [ [2;2;0]+0.5*(Rand3./norms(Rand3,2))   ;    zeros(3,KK)   ];
    
    if mod(iter, Batch_num)==1
        Vd=cell(1,Batch_num);Vd{1}=zeros(1,len+lenSTL);
        Sd=cell(1,Batch_num);Sd{1}=zeros(1,len+lenSTL);
        Vd1=cell(1,Batch_num);Vd1{1}=zeros(1,len+lenSTL);
        Sd1=cell(1,Batch_num);Sd1{1}=zeros(1,len+lenSTL);
        Vd2=cell(1,Batch_num);Vd2{1}=zeros(1,len+lenSTL);
        Sd2=cell(1,Batch_num);Sd2{1}=zeros(1,len+lenSTL);
        index=0;
    end
    index=index+1;
    
    
    parfor initNum=1:KK
        s0 = s0s(:,initNum);
        Ggrad1(initNum,:)=[Back_prop_Perf(s0, param{iter}, M, T,A,B) 0 zeros(1,lenSTL-1)];
        Ggrad2(initNum,:)=Back_prop_STL(s0, param{iter}, lam_param{iter}, beta_param{iter}, M, T,A,B);
    end
    
    
    
    jj=[];
    for i=1:KK
        if isnan(sum(Ggrad1(i,:)))
           jj=[jj,i];
        end
    end
    Ggrad1(jj,:)=[];
    
    if size(Ggrad1,1)>=1
        Norms1=norms(Ggrad1',2);
        i1=find(Norms1==max(Norms1)); i1=i1(1);
        grad1 = Ggrad1(i1,:);
    else 
        error('Every thing is NaN');
    end
    
    
    jj=[];
    for i=1:KK
        if isnan(sum(Ggrad2(i,:)))
           jj=[jj,i];
        end
    end
    Ggrad2(jj,:)=[];
    if size(Ggrad2,1)>=1
        Norms2=norms(Ggrad2',2);
        i2=find(Norms2==max(Norms2)); i2=i2(1);
        grad2 = Ggrad2(i2,:);
    else
        error('Every thing is NaN');
    end
    
    
    fprintf('Iteration %d: G1 = %s. G2 = %s.\n\n', iter, mat2str(grad1,1), mat2str(grad2,1));
    
    
    Vd1{index+1} = Beta1_Perf*Vd{index}+(1-Beta1_Perf)*grad1   ;
    Sd1{index+1} = Beta2_Perf*Sd{index}+(1-Beta2_Perf)*grad1.^2;
    Vd2{index+1} = Beta1_STL*Vd{index}+(1-Beta1_STL)*grad2   ;
    Sd2{index+1} = Beta2_STL*Sd{index}+(1-Beta2_STL)*grad2.^2;
    
    
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
    
    Rand3=2*rand(3,1)-1;
    s0 = FF*[  [2;2;0]+0.5*Rand3/norm(Rand3)  ; zeros(3,1)   ]; 
    A=As+(epsi*(2*rand-1))*eye(6);
    
    
    J2_discuss= STL_obj(s0,param{iter}, lam_param{iter}, beta_param{iter}, M, T,A,B);
    J1next1=Perf_obj(s0,paramnext1,M, T,A,B);
    J1next1p=Perf_obj(s0,paramnext1p,M, T,A,B);
    J1next2=Perf_obj(s0,paramnext2,M, T,A,B);
    [J2next1, Bsnext1]=STL_obj(s0,paramnext1, lam_paramnext1, beta_paramnext1, M, T,A,B);
    [J2next1p, Bsnext1p]=STL_obj(s0,paramnext1p, lam_paramnext1p, beta_paramnext1p, M, T,A,B);
    [J2next2, Bsnext2]=STL_obj(s0,paramnext2, lam_paramnext2, beta_paramnext2, M, T,A,B);
    if J2_discuss<rho
        if J2next1>J2_discuss
            param{iter+1}=paramnext1;
            lam_param{iter+1}=lam_paramnext1;
            beta_param{iter+1}=beta_paramnext1;
            J1{iter+1}=J1next1;
            J2{iter+1}=J2next1; Bs{iter+1}=Bsnext1;
%             Grad(iter,:)  =grad1;
            Vd{index+1}=Vd1{index+1}; Sd{index+1}=Sd1{index+1};
        
        else
            param{iter+1}=paramnext2;
            lam_param{iter+1}=lam_paramnext2;
            beta_param{iter+1}=beta_paramnext2;
            J1{iter+1}=J1next2;
            J2{iter+1}=J2next2; Bs{iter+1}=Bsnext2;
%             Grad(iter,:)  =grad2;
            Vd{index+1}=Vd2{index+1}; Sd{index+1}=Sd2{index+1};
        end
    else
        param{iter+1}=paramnext1p;
        lam_param{iter+1}=lam_paramnext1p;
        beta_param{iter+1}=beta_paramnext1p;
        J1{iter+1}=J1next1p;
        J2{iter+1}=J2next1p; Bs{iter+1}=Bsnext1p;
        Vd{index+1}=Vd1{index+1}; Sd{index+1}=Sd1{index+1};
%         Grad(iter,:)=grad1;
    end
end
Run_time=toc;

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
for j=1:len+lenSTL
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


Ma=[min(J1s(floor(0.9*iter):end)),max(J1s(floor(0.9*iter):end))];
slope=1/(2*(Ma(1)-Ma(2)));
intersect = 1- Ma(1)*slope;


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
    
    surface=surf(X2+FF*5,Y2+FF*5,Z2,'FaceAlpha',0.2);
    surface.EdgeColor = 'none';
    
    
    hold on
    [X,Y,Z]=sphere;
    r=sqrt(1.5)*FF;
    X2 = X * r;
    Y2 = Y * r;
    Z2 = Z * r;
    
    surface=surf(X2+FF*2,Y2+FF*8,Z2,'FaceAlpha',0.2);
    surface.EdgeColor = 'none';
    
    
    
    hold on
    [X,Y,Z]=sphere;
    r=sqrt(1.5)*FF;
    X2 = X * r;
    Y2 = Y * r;
    Z2 = Z * r;
    
    surface=surf(X2+FF*8,Y2+FF*2,Z2,'FaceAlpha',0.2);
    surface.EdgeColor = 'none';
    
    
    hold on
    [X,Y,Z]=sphere;
    r=sqrt(1.5)*FF;
    X2 = X * r;
    Y2 = Y * r;
    Z2 = Z * r;
    
    surface=surf(X2+FF*8,Y2+FF*8,Z2-FF*3,'FaceAlpha',0.2);
    surface.EdgeColor = 'none';
    
    
    
    axis equal
    
end


hold on




Param=param{iter};

rrs     = linspace(0,1,10);
pphi    = linspace(0,2*pi,10);
ttheta  = linspace(-pi/2,pi/2,8);
deltas  = linspace(-0.01,0.01,6);


for ii=1:length(rrs) 
    for jj=1:length(pphi)
        for kk=1:length(ttheta)
            for ll=1:length(deltas)
                A=As+deltas(ll)*eye(6);
                initxx = 2+ 0.5*rrs(ii)*cos(ttheta(kk))*cos(pphi(jj));
                inityy = 2+ 0.5*rrs(ii)*cos(ttheta(kk))*sin(pphi(jj));
                initzz =    0.5*rrs(ii)*sin(ttheta(kk));
                s0 = FF* [initxx;inityy;initzz;zeros(3,1)];
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
                    title('Violating Trajectories')
                    color = '-r.' ;
                    plot3([s0(1,1) xx(1,:)],[s0(2,1) xx(2,:)], color, 'Linewidth', 0.75);
                
                else
                    figure(5)
                    title('Satisfactory Trajectories')
                    
                    reward=Perf_obj(s0, Param, M, T, A, B);
                    plot3([s0(1,1) xx(1,:)],[s0(2,1) xx(2,:)],[s0(3,1) xx(3,:)], 'color',0.9*(slope*reward+intersect)*[0,1,0], 'Linewidth', 1);
                    
                end
            
            end
        end
    end
end

hold off

