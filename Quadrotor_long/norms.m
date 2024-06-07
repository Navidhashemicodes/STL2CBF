function  y = norms(x,j)
n=size(x,2);
y=zeros(1,n);
for i=1:n
    y(i)=norm(x(:,i),j);
end
end