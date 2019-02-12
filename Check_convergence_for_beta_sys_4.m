% Check that the number of simulations is enough for convergence
Expected_costs=convergence(j).Exp_costs;
ix = randperm(length(Expected_costs));
Expected_costs_shuffled = Expected_costs(ix);

parfor i_sim=1:input.N_sim
    meanEC(i_sim)=mean(Expected_costs_shuffled(1:i_sim));
    meaEC_plus(i_sim)=mean(Expected_costs_shuffled(1:i_sim))+std(Expected_costs_shuffled(1:i_sim))/sqrt(i_sim);
    meaEC_minus(i_sim)=mean(Expected_costs_shuffled(1:i_sim))-std(Expected_costs_shuffled(1:i_sim))/sqrt(i_sim);
end

plot(1:input.N_sim,meanEC,'.-r')
hold on
plot(1:input.N_sim,meaEC_plus,'.-k')
plot(1:input.N_sim,meaEC_minus,'.-k')

legend('\mu','\mu+\sigma/N^{0.5}','\mu-\sigma/N^{0.5}')
xlabel('Number of simulations N_{sim}')
ylabel('E[C_{tot}](N_{sim}')
    
    
    
    
    
    