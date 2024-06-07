clear
clc
close all
load('rho=0.3.mat')
iter=22500;



SS1=zeros(iter-1,len);
SS2=zeros(iter-1,1);
SS3=zeros(iter-1,12);
for j=1:iter-1
    SS1(j,:)=param{j}; 
    SS2(j,:)=lam_param{j};
    SS3(j,:)=beta_param{j}(1,[1:10,21,22]);
end


figure
for j=1:len
    plot(SS1(:,j),'Linewidth',2)
    hold on
end


figure
for j=1:1
    plot(SS2(:,j),'Linewidth',2)
    hold on
end

figure
for j=1:10+2
    plot(SS3(:,j),'Linewidth',2)
    hold on
end



% clear
% clc
% figure
% load('candidate2_0_00rho_cutafter25000.mat')
% subplot(2,2,1)
% Param=param{1};
% plotTraj([0.5 0.6],[0.5 0.6],0.01,M,P,Param,T);
% 
% NNN=10000;
% for j=1:NNN
%     thet=2*pi*j/NNN;
%     r=[cos(thet);sin(thet)];
%     x1(:,j)=[5;5]+sqrt(1.5)*r;
%     x2(:,j)=[2;8]+sqrt(1.5)*r;
%     x3(:,j)=[8;2]+sqrt(1.5)*r;
% end
% 
% plot(x1(1,:),x1(2,:));
% hold on
% plot(x2(1,:),x2(2,:));
% hold on
% plot(x3(1,:),x3(2,:));
% hold on
% plot(8,8,'*');
% axis equal
% xlim([-1,10])
% ylim([-1,10])
% 
% 
% 
% subplot(2,2,2)
% dd=true;
% ii=iter;
% ll=0;
% while dd
%     ii=ii-1;
%     if J2{ii}>0.999*rho
%         ll=ll+1;
%         J1pos(ll)=J1s(ii);indices(ll)=ii;
%     end
%     if ii==2
%         dd=false;
%     end
% end
% index=find(J1pos==max(J1pos));index=index(1);
% Param=param{indices(index)};
% 
% plotTraj([0.5 0.6],[0.5 0.6],0.01,M,P,Param,T);
% 
% NNN=10000;
% for j=1:NNN
%     thet=2*pi*j/NNN;
%     r=[cos(thet);sin(thet)];
%     x1(:,j)=[5;5]+sqrt(1.5)*r;
%     x2(:,j)=[2;8]+sqrt(1.5)*r;
%     x3(:,j)=[8;2]+sqrt(1.5)*r;
% end
% 
% plot(x1(1,:),x1(2,:));
% hold on
% plot(x2(1,:),x2(2,:));
% hold on
% plot(x3(1,:),x3(2,:));
% hold on
% plot(8,8,'*');
% 
% axis equal
% xlim([-1,10])
% ylim([-1,10])
% 
% 
% load('candidate2_0_3rho.mat')
% subplot(2,2,3)
% dd=true;
% ii=iter;
% ll=0;
% while dd
%     ii=ii-1;
%     if J2{ii}>0.999*rho
%         ll=ll+1;
%         J1pos(ll)=J1s(ii);indices(ll)=ii;
%     end
%     if ii==2
%         dd=false;
%     end
% end
% index=find(J1pos==max(J1pos));index=index(1);
% Param=param{indices(index)};
% 
% plotTraj([0.5 0.6],[0.5 0.6],0.01,M,P,Param,T);
% 
% NNN=10000;
% for j=1:NNN
%     thet=2*pi*j/NNN;
%     r=[cos(thet);sin(thet)];
%     x1(:,j)=[5;5]+sqrt(1.5)*r;
%     x2(:,j)=[2;8]+sqrt(1.5)*r;
%     x3(:,j)=[8;2]+sqrt(1.5)*r;
% end
% 
% plot(x1(1,:),x1(2,:));
% hold on
% plot(x2(1,:),x2(2,:));
% hold on
% plot(x3(1,:),x3(2,:));
% hold on
% plot(8,8,'*');
% 
% axis equal
% xlim([-1,10])
% ylim([-1,10])
% 
% 
% load('candidate2_0_75rho.mat')
% subplot(2,2,4)
% dd=true;
% ii=iter;
% ll=0;
% while dd
%     ii=ii-1;
%     if J2{ii}>0.999*rho
%         ll=ll+1;
%         J1pos(ll)=J1s(ii);indices(ll)=ii;
%     end
%     if ii==2
%         dd=false;
%     end
% end
% index=find(J1pos==max(J1pos));index=index(1);
% Param=param{indices(index)};
% 
% plotTraj([0.5 0.6],[0.5 0.6],0.01,M,P,Param,T);
% 
% NNN=10000;
% for j=1:NNN
%     thet=2*pi*j/NNN;
%     r=[cos(thet);sin(thet)];
%     x1(:,j)=[5;5]+sqrt(1.5)*r;
%     x2(:,j)=[2;8]+sqrt(1.5)*r;
%     x3(:,j)=[8;2]+sqrt(1.5)*r;
% end
% 
% plot(x1(1,:),x1(2,:));
% hold on
% plot(x2(1,:),x2(2,:));
% hold on
% plot(x3(1,:),x3(2,:));
% hold on
% plot(8,8,'*');
% axis equal
% xlim([-1,10])
% ylim([-1,10])

% clear
% clc
% load('candidate2_0_00rho_cutafter25000.mat')
% subplot(1,3,1)
% dd=true;
% ii=iter;
% ll=0;
% while dd
%     ii=ii-1;
%     if J2{ii}>0.999*rho
%         ll=ll+1;
%         J1pos(ll)=J1s(ii);indices(ll)=ii;
%     end
%     if ii==2
%         dd=false;
%     end
% end
% index=find(J1pos==max(J1pos));index=index(1);
% Param=param{indices(index)};
% Lam=lam_param{indices(index)}+10;
% Beta=beta_param{indices(index)};
% NN=100000;
% for iii=1:NN
%     s0=[0.5+0.1*rand(2,1);pi/2];
%     [~, TheB]=STL_obj(s0,Param, Lam,Beta,M,P,T);
%     STL1(iii)=min(TheB);
% end
% plot(STL1(1:5000),'green')
% 
% 
% Param=param{indices(index)};
% 
% clearvars -except STL1
% load('candidate2_0_3rho.mat')
% subplot(1,3,2)
% dd=true;
% ii=iter;
% ll=0;
% while dd
%     ii=ii-1;
%     if J2{ii}>0.999*rho
%         ll=ll+1;
%         J1pos(ll)=J1s(ii);indices(ll)=ii;
%     end
%     if ii==2
%         dd=false;
%     end
% end
% index=find(J1pos==max(J1pos));index=index(1);
% Param=param{indices(index)};
% Lam=lam_param{indices(index)}+10;
% Beta=beta_param{indices(index)};
% NN=100000;
% for iii=1:NN
%     s0=[0.5+0.1*rand(2,1);pi/2];
%     [~, TheB]=STL_obj(s0,Param, Lam,Beta,M,P,T);
%     STL2(iii)=min(TheB);
% end
% plot(STL2(1:5000),'green')
% clearvars -except STL2 STL1
% clc
% load('candidate2_0_75rho.mat')
% subplot(1,3,3)
% dd=true;
% ii=iter;
% ll=0;
% while dd
%     ii=ii-1;
%     if J2{ii}>0.999*rho
%         ll=ll+1;
%         J1pos(ll)=J1s(ii);indices(ll)=ii;
%     end
%     if ii==2
%         dd=false;
%     end
% end
% index=find(J1pos==max(J1pos));index=index(1);
% Param=param{indices(index)};
% Lam=lam_param{indices(index)}+10;
% Beta=beta_param{indices(index)};
% NN=100000;
% for iii=1:NN
%     s0=[0.5+0.1*rand(2,1);pi/2];
%     [~, TheB]=STL_obj(s0,Param, Lam,Beta,M,P,T);
%     STL3(iii)=min(TheB);
% end
% 
% plot(STL3(1:5000),'green')
% 
% STLM1=mean(STL1);
% STLC1=cov(STL1);
% disp(STLM1);
% disp(STLC1);
% 
% STLM2=mean(STL2);
% STLC2=cov(STL2);
% disp(STLM2);
% disp(STLC2);
% 
% STLM3=mean(STL3);
% STLC3=cov(STL3);
% disp(STLM3);
% disp(STLC3);
% 
% figure
% t=-2:0.001:2;
% y1=exp(-(t-STLM1).^2/STLC1);
% y2=exp(-(t-STLM2).^2/STLC2);
% y3=exp(-(t-STLM3).^2/STLC3);
% 
% plot(t,y1)
% hold on
% plot(t,y2)
% hold on
% plot(t,y3)

