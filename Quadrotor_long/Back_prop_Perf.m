function gradient_Perf = Back_prop_Perf(s0, param, M, T,A,B,FF)

n=size(s0,1);
P=size(B,2);
index=0;
Weights{1}=reshape(param(1,index+1:index+M*(n+1)), [n+1, M])';  
index=index+M*(n+1);
Biases{1}= param(1,index+1:index+M)';
index=index+M;
Weights{2}=reshape(param(1,index+1:index+M*P), [M, P])';
index=index+M*P;
Biases{2}=param(1,index+1:index+P)';



discount=0.9;
dest=FF*[8;8;-3];

Sab{1}=Weights{1}*[s0;0]+Biases{1};
Sa{1}=tanh(Sab{1});
a{1}=Weights{2}*Sa{1}+Biases{2};
aa=a{1};
uu=[tan(0.1*tanh(aa(1,1)/10)); tan(0.1*tanh(aa(2,1)/10));   2*tanh(aa(3,1)/10) ];
s{1}=A*s0+B*uu;
D_f_a{1}=[B(:,1)*((0.1)*(0.1)*(1-tanh(aa(1,1)/10)^2)/(1+tan(0.1*tanh(aa(1,1)/10))^2))...
          B(:,2)*((0.1)*(0.1)*(1-tanh(aa(2,1)/10)^2)/(1+tan(0.1*tanh(aa(2,1)/10))^2))...
          B(:,3)*((2)*(0.1)*(1-tanh(aa(3,1)/10)^2))  ];
D_S_step{1}= 0;
D_J_S{1}=zeros(1,n);    %%% this is only an initialization. This value is to be counted with back propagation
Omega{1}=[-2*(s{1}(1:3,1)-dest)'*reward(s{1},FF)/((FF^2)*36)  zeros(1,3)] ;
for i=2:T
    Sab{i}=Weights{1}*[s{i-1};i-1]+Biases{1};
    Sa{i}=tanh(Sab{i});
    a{i}=Weights{2}*Sa{i}+Biases{2};
    aa=a{i};
    uu=[tan(0.1*tanh(aa(1,1)/10)); tan(0.1*tanh(aa(2,1)/10));   2*tanh(aa(3,1)/10) ];
    s{i}=A*s{i-1}+B*uu;
    D_f_a{i}=[B(:,1)*((0.1)*(0.1)*(1-tanh(aa(1,1)/10)^2)/(1+tan(0.1*tanh(aa(1,1)/10))^2))...
              B(:,2)*((0.1)*(0.1)*(1-tanh(aa(2,1)/10)^2)/(1+tan(0.1*tanh(aa(2,1)/10))^2))...
              B(:,3)*((2)*(0.1)*(1-tanh(aa(3,1)/10)^2))  ];
    D_f_s=A;
    D_S_step{i}=D_f_s+D_f_a{i}*Weights{2}*diag(1-Sa{i}.^2)*Weights{1}(:,1:n);
    D_J_S{i}=zeros(1,n);   %%% this is only an initialization, This value is to be counted with back propagation
    Omega{i}=[-2*(s{i}(1:3,1)-dest)'*reward(s{i},FF)/((FF^2)*36)   zeros(1,3)];
end


for i=1:T
    E{i}=discount^(i);
end

%%%%%%% Start of Back Propagation
D_J_S{T}=E{T}*Omega{T};
for i=T-1:-1:1
    D_J_S{i}=E{i}*Omega{i}+D_J_S{i+1}*D_S_step{i+1};
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
        D_S_W1{i}{m}=D_f_a{i}*Weights{2}*D_Sa_W1{i}{m}+A*D_S_W1{i-1}{m};
        D_S_B1{i}{m}=D_f_a{i}*Weights{2}*D_Sa_B1{i}{m}+A*D_S_B1{i-1}{m};
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
        D_S_W2{i}{p}=D_f_a{i}*kron(I(:,p) , Sa{i}')+A*D_S_W2{i-1}{p};
        D_S_B2{i}{p}=D_f_a{i}*I(:,p)+A*D_S_B2{i-1}{p};
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

gradient_Perf=[D_J_param1,D_J_param2];


end



function z=reward(s,FF)
s=s(1:3,1);
dest=FF*[8;8;-3];
F=(s-dest)'*(s-dest);
z=10*exp(-F/((FF^2)*36));
end