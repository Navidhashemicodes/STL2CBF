function gradient_STL = Back_prop_STL(s0, param, lam_param, beta, q, M, P, T)

n=size(s0,1);

index=0;
Weights{1}=reshape(param(1,index+1:index+M*(n+1)), [n+1, M])';  
index=index+M*(n+1);
Biases{1}= param(1,index+1:index+M)';
index=index+M;
Weights{2}=reshape(param(1,index+1:index+M*P), [M, P])';
index=index+M*P;
Biases{2}=param(1,index+1:index+P)';



eta=lam_param^2+1;


Q1=1.5*eye(2);
P1=inv(Q1);
center1=[5;5];

Q2=1.5*eye(2);
P2=inv(Q2);
center2=[2;8];

Q3=1.5*eye(2);
P3=inv(Q3);
center3=[8;2];




Sab{1}=Weights{1}*[s0;0]+Biases{1};
Sa{1}=tanh(Sab{1});

a{1}=Weights{2}*Sa{1}+Biases{2};
s{1}=Dynamics(s0,a{1},q);

D_f_a{1}=Dyn_a(s0,a{1});
D_f_s{1}=[];
D_S_step{1}=0;
exp_sum= exp(-eta*CBF1(s{1}));
exp_bsum= CBF1(s{1})*exp(-eta*CBF1(s{1}));
D_J_S{1}=zeros(1,n);    %%% this is only an initialization. This value is to be counted with back propagation

F1=(s{1}(1:2,1)-center1)'*P1*(s{1}(1:2,1)-center1);
Omega1{1}=[2*(s{1}(1:2,1)-center1)'*P1*exp(1-F1) 0];

for i=2:T
    Sab{i}=Weights{1}*[s{i-1};i-1]+Biases{1};
    Sa{i}=tanh(Sab{i});
    a{i}=Weights{2}*Sa{i}+Biases{2};
    s{i}=Dynamics(s{i-1},a{i},q);
    D_f_a{i}=Dyn_a(s{i-1},a{i});
    D_f_s{i}=Dyn_s(s{i-1},a{i},q);
    D_S_step{i}=D_f_s{i}+D_f_a{i}*Weights{2}*diag(1-Sa{i}.^2)*Weights{1}(:,1:n);
    exp_sum=exp_sum + exp(-eta*CBF1(s{i}));
    exp_bsum= exp_bsum + CBF1(s{i})*exp(-eta*CBF1(s{i}));
    D_J_S{i}=zeros(1,n);   %%% this is only an initialization, This value is to be counted with back propagation
    F1=(s{i}(1:2,1)-center1)'*P1*(s{i}(1:2,1)-center1);
    Omega1{i}=[2*(s{i}(1:2,1)-center1)'*P1*exp(1-F1) 0];
end


sq_sum1=0;
beta_sum1=0;
for i=1:10
    sq_sum1=sq_sum1 + beta(i)^2*CBF2(s{i});
    beta_sum1=beta_sum1+beta(i)^2;
end
w_bar1=sq_sum1/beta_sum1;


for i=1:10
    F2=(s{i}(1:2,1)-center2)'*P2*(s{i}(1:2,1)-center2);
    Omega2{i}= ((beta(21)^2)/(beta(21)^2+beta(22)^2))*(beta(i)^2/beta_sum1)*[-2*(s{i}(1:2,1)-center2)'*P2   0];
end


sq_sum2=0;
beta_sum2=0;
for i=1:10
    sq_sum2=sq_sum2 + beta(i+10)^2*CBF3(s{i});
    beta_sum2=beta_sum2+beta(i+10)^2;
end
w_bar2=sq_sum2/beta_sum2;


for i=1:10
    F3=(s{i}(1:2,1)-center3)'*P3*(s{i}(1:2,1)-center3);
    Omega3{i}=((beta(22)^2)/(beta(21)^2+beta(22)^2))*(beta(i+10)^2/beta_sum2)*[-2*(s{i}(1:2,1)-center3)'*P3   0];
end

sq_sum=(beta(21)^2)*w_bar1  +  (beta(22)^2)*w_bar2;
beta_sum=beta(21)^2  +  beta(22)^2;
w_bar= sq_sum   /   beta_sum  ;
exp_sum=exp_sum+exp(-eta*w_bar);
exp_bsum=exp_bsum+w_bar*exp(-eta*w_bar);






for i=1:T
    E1{i}=exp(-eta*CBF1(s{i}))/exp_sum;
end

for i=1:10
    E2{i}=exp(-eta*w_bar)/exp_sum;
end





%%%%%%% Start of Back Propagation
D_J_S{T}=E1{T}*Omega1{T};
for i=T-1:-1:1
    EO=E1{i}*Omega1{i};
    if i<=10
        EO=EO+E2{i}*(Omega2{i}+Omega3{i});
    end
    D_J_S{i}=EO+D_J_S{i+1}*D_S_step{i+1};
end
%%%%%%  End   of Back Propagation

%%%%%% First Layer

I=eye(M);
for m=1:M
    D_Sa_W1{1}{m}=kron(I(:,m), ((1-(Sa{1}(m,1))^2)*[s0;0])' );  %%% NN is designed for the time steps not the real times in seconds.
    D_Sa_B1{1}{m}=(1-(Sa{1}(m,1))^2)*I(:,m);
    D_S_W1{1}{m}=D_f_a{1}*Weights{2}*D_Sa_W1{1}{m};
    D_S_B1{1}{m}=D_f_a{1}*Weights{2}*D_Sa_B1{1}{m};
    for i=2:T
        D_Sa_W1{i}{m}=kron(I(:,m), ((1-(Sa{i}(m,1))^2)*[s{i-1};i-1])' );
        D_Sa_B1{i}{m}=(1-(Sa{i}(m,1))^2)*I(:,m);
        D_S_W1{i}{m}=D_f_a{i}*Weights{2}*D_Sa_W1{i}{m}+D_f_s{i}*D_S_W1{i-1}{m};
        D_S_B1{i}{m}=D_f_a{i}*Weights{2}*D_Sa_B1{i}{m}+D_f_s{i}*D_S_B1{i-1}{m};
    end
end

D_J_Weight1=[];
D_J_Biases1=[];
for m=1:M
    D_J_W1{m}=zeros(1,n+1);
    D_J_B1{m}=zeros(1,1);
    for i=1:T
        D_J_W1{m}=D_J_W1{m}+D_J_S{i}*D_S_W1{i}{m};
        D_J_B1{m}=D_J_B1{m}+D_J_S{i}*D_S_B1{i}{m};
    end
    D_J_Weight1=[D_J_Weight1,D_J_W1{m}];
    D_J_Biases1=[D_J_Biases1,D_J_B1{m}];
end
D_J_param1=[D_J_Weight1,D_J_Biases1];





%%%%% Second Layer

I=eye(P);
for p=1:P
    D_S_W2{1}{p}=D_f_a{1}*kron(I(:,p) , Sa{1}');
    D_S_B2{1}{p}=D_f_a{1}*I(:,p);
    for i=2:T
        D_S_W2{i}{p}=D_f_a{i}*kron(I(:,p) , Sa{i}')+D_f_s{i}*D_S_W2{i-1}{p};
        D_S_B2{i}{p}=D_f_a{i}*I(:,p)+D_f_s{i}*D_S_B2{i-1}{p};
    end
end


D_J_Weight2=[];
D_J_Biases2=[];
for p=1:P
    D_J_W2{p}=zeros(1,M);
    D_J_B2{p}=zeros(1,1);
    for i=1:T
        D_J_W2{p}=D_J_W2{p}+D_J_S{i}*D_S_W2{i}{p};
        D_J_B2{p}=D_J_B2{p}+D_J_S{i}*D_S_B2{i}{p};
    end
    D_J_Weight2=[D_J_Weight2,D_J_W2{p}];
    D_J_Biases2=[D_J_Biases2,D_J_B2{p}];
end

D_J_param2=[D_J_Weight2,D_J_Biases2];


D_J_lam= 2*lam_param*(  exp_bsum/(eta*exp_sum)  +  log(exp_sum)/(eta^2)   );
for i=1:10
    D_J_beta(i)=(exp(-eta*w_bar)/exp_sum)*((beta(21)^2)/(beta(21)^2+beta(22)^2))*(  (2*beta(i)*CBF2(s{i})/beta_sum1)  -  (2*beta(i)*sq_sum1/(beta_sum1^2))  );
end

for i=1:10
    D_J_beta(i+10)=(exp(-eta*w_bar)/exp_sum)*((beta(22)^2)/(beta(21)^2+beta(22)^2))*(  (2*beta(i+10)*CBF3(s{i})/beta_sum2)  -  (2*beta(i+10)*sq_sum2/(beta_sum2^2))  );
end

D_J_beta(21)=(exp(-eta*w_bar)/exp_sum)*(  (2*beta(21)*w_bar1/beta_sum)  -  (2*beta(21)*sq_sum/(beta_sum^2))  );
D_J_beta(22)=(exp(-eta*w_bar)/exp_sum)*(  (2*beta(22)*w_bar2/beta_sum)  -  (2*beta(22)*sq_sum/(beta_sum^2))  );

gradient_STL=[D_J_param1,D_J_param2, D_J_lam, D_J_beta];

end



function z=CBF1(s)
s=s(1:2,1);
Q1=1.5*eye(2);
P1=inv(Q1);
center1=[5;5];
F1=(s-center1)'*P1*(s-center1);
z=1-exp(1-F1);
end

function z=CBF2(s)
s=s(1:2,1);
Q2=1.5*eye(2);
P2=inv(Q2);
center2=[2;8];
F2=(s-center2)'*P2*(s-center2);
z=1-F2;
end

function z=CBF3(s)
s=s(1:2,1);
Q3=1.5*eye(2);
P3=inv(Q3);
center3=[8;2];
F3=(s-center3)'*P3*(s-center3);
z=1-F3;
end