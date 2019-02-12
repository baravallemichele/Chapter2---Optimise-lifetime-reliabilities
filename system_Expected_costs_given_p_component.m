function [ Expected_costs ] = system_Expected_costs_given_p_component( p_component,beta_system,input )
%Calcualte the system expected csots given p_component and beta_system
%   input
%       p : design parameter
%Initialize variables
Expected_costs=zeros(1,input.N_sim);
    parfor i_sim=1:input.N_sim
             disp(i_sim/input.N_sim)
             % for series or 1-component systes the failure of the component leads to the consequences H_sys
             if  or(input.sys_type(i_sim)==0, input.num_el(i_sim)==1)       
                      Contruction_costs= (1+input.num_el(i_sim).*input.mCI(i_sim).*p_component(i_sim));
                      Obsolescence_costs=(1+input.num_el(i_sim).*input.mCI(i_sim).*p_component(i_sim)+input.mA(i_sim)).*input.mOmega./log(1+input.mGamma_reftime);
                      Failure_costs=(1+input.num_el(i_sim).*input.mCI(i_sim).*p_component(i_sim)+input.mHsys(i_sim)).*normcdf(-beta_system)./log(1+input.mGamma_reftime);
                      Expected_costs(i_sim)=Contruction_costs+Obsolescence_costs+Failure_costs;
             else
                      %parallel system:
                      %         failure of one component but survival of system    --> (C0/n+CI*p+Cdir)Pfcomp*Pfcomp given system survival
                      %         failure of one component and failure of the system --> (C0+n*CI*p+Hsys)*Pfsys given component failure
                      Prob_comp_fail_sys_surv = input.Prob_comp_failure_system_survival(i_sim);
                      Contruction_costs= (1+input.num_el(i_sim).*input.mCI(i_sim).*p_component(i_sim));
                      Obsolescence_costs=(1+input.num_el(i_sim).*input.mCI(i_sim).*p_component(i_sim)+input.mA(i_sim)).*input.mOmega./log(1+input.mGamma_reftime);
                      System_Failure_costs=(1+input.num_el(i_sim).*input.mCI(i_sim).*p_component(i_sim)+input.mHsys(i_sim)).*normcdf(-beta_system)./log(1+input.mGamma_reftime);
                      Component_Failure_System_Survival_costs=(1./input.num_el(i_sim)+input.mCI(i_sim).*p_component(i_sim)+input.mCdir(i_sim)).*Prob_comp_fail_sys_surv./log(1+input.mGamma_reftime);
                      Expected_costs(i_sim)=Contruction_costs+Obsolescence_costs+System_Failure_costs+Component_Failure_System_Survival_costs;
             end
    end
end
    


