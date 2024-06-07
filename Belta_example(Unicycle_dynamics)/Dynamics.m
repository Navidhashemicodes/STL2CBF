function snext=Dynamics(s,a,q)
a(1,1)=0.5*a(1,1);
a(2,1)=0.5*a(2,1);

snext=zeros(3,1);
snext(1,1)=q*s(1,1) + (  2 * sigmoid(a(1,1))  /   tanh(a(2,1))  )  *  (  sin(s(3,1)+0.5*tanh(a(2,1)))  -  sin(s(3,1))  );
snext(2,1)=q*s(2,1) + (  2 * sigmoid(a(1,1))  /   tanh(a(2,1))  )  *  (  cos(s(3,1))  -  cos(s(3,1)+0.5*tanh(a(2,1)))  );
snext(3,1)=q*s(3,1) + 0.5*tanh(a(2,1));

end