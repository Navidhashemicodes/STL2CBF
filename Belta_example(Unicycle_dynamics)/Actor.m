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