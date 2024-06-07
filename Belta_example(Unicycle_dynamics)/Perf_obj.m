function J1 = Perf_obj(s0, param, q, M, P, T)

n=size(s0,1);

index=0;
Weights{1}=reshape(param(1,index+1:index+M*(n+1)), [n+1, M])';  
index=index+M*(n+1);
Biases{1}= param(1,index+1:index+M)';
index=index+M;
Weights{2}=reshape(param(1,index+1:index+M*P), [M, P])';
index=index+M*P;
Biases{2}=param(1,index+1:index+P)';


discount=0.9;

Sab{1}=Weights{1}*[s0;0]+Biases{1};
Sa{1}=tanh(Sab{1});
a{1}=Weights{2}*Sa{1}+Biases{2};
s{1}=Dynamics(s0,a{1},q);
J1=discount^(1)*reward(s{1});
for i=2:T
    Sab{i}=Weights{1}*[s{i-1};i-1]+Biases{1};
    Sa{i}=tanh(Sab{i});
    a{i}=Weights{2}*Sa{i}+Biases{2};
    s{i}=Dynamics(s{i-1},a{i},q);
    J1=J1+discount^(i)*reward(s{i});
end
end




function z=reward(s)
s=s(1:2,1);
dest=[8;8];
F=(s-dest)'*(s-dest);
z=10*exp(-F/36);
end