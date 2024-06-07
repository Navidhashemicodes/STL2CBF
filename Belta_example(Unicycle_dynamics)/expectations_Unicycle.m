function [JP, J_STL_h, J_STL_s] = expectations_Unicycle(initxx,inityy,initalpha, model_UC, numbers, M, P, Param, Lam, Bet, T)

rng(0)


rob_hard=zeros(1,numbers);
rob_soft=zeros(1,numbers);
rewards=zeros(1,numbers);
for ii=1:numbers
    q= model_UC(1)+rand*(model_UC(2)-model_UC(1));
    s0 = [initxx(1);inityy(1);initalpha(1)]+rand(3,1).*[initxx(2)-initxx(1);  inityy(2)-inityy(1);  initalpha(2)-initalpha(1)];
    
    rob_hard(ii) = robustness_Belta(s0,Param, q, M, P, T);
    [rob_soft(ii), ~] = STL_obj(s0, Param, Lam, Bet, q, M, P, T);
    rewards(ii) = Perf_obj(s0, Param, q, M, P, T);
end

JP.minvalue=min(rewards);JP.maxvalue=max(rewards);JP.meanvalue=mean(rewards);JP.covar=cov(rewards);
J_STL_h.minvalue=min(rob_hard);J_STL_h.maxvalue=max(rob_hard);J_STL_h.meanvalue=mean(rob_hard);J_STL_h.covar=cov(rob_hard);
J_STL_s.minvalue=min(rob_soft);J_STL_s.maxvalue=max(rob_soft);J_STL_s.meanvalue=mean(rob_soft);J_STL_s.covar=cov(rob_soft);


figure(1)
hist(rewards,1000,'normalized')

figure(2)
hist(rob_soft,1000,'normalized')

figure(3)
hist(rob_hard,1000,'normalized')

end