function [JP, J_STL_h, J_STL_s] = expectations_Quadrotor(rrs,pphis,tthetas, deltas, A_mean, B, numbers, M, Param, Lam, Bet, T)

rng(0)
FF=0.1*0.5*0.25;


rob_hard=zeros(1,numbers);
rob_soft=zeros(1,numbers);
rewards=zeros(1,numbers);
for ii=1:numbers
    rr=rrs(1)+rand*(rrs(2)-rrs(1));
    pphi=pphis(1)+rand*(pphis(2)-pphis(1));
    ttheta=tthetas(1)+rand*(tthetas(2)-tthetas(1));
    delta_model= deltas(1)+rand*(deltas(2)-deltas(1));
    initxx = 2+ 0.5*rr*cos(ttheta)*cos(pphi);
    inityy = 2+ 0.5*rr*cos(ttheta)*sin(pphi);
    initzz =    0.5*rr*sin(ttheta);
    s0 = FF* [initxx;inityy;initzz;zeros(3,1)];
    A=A_mean+delta_model*eye(6);
    
    rob_hard(ii) = robustness_Quadrotor(s0,Param, M, T, A, B);

    [rob_soft(ii), ~] = STL_obj(s0, Param, Lam, Bet, M, T, A, B);
    rewards(ii) = Perf_obj(s0, Param, M, T, A, B);
end

JP.minvalue=min(rewards);JP.maxvalue=max(rewards);JP.meanvalue=mean(rewards);JP.covar=cov(rewards);
J_STL_h.minvalue=min(rob_hard);J_STL_h.maxvalue=max(rob_hard);J_STL_h.meanvalue=mean(rob_hard);J_STL_h.covar=cov(rob_hard);
J_STL_s.minvalue=min(rob_soft);J_STL_s.maxvalue=max(rob_soft);J_STL_s.meanvalue=mean(rob_soft);J_STL_s.covar=cov(rob_soft);

end