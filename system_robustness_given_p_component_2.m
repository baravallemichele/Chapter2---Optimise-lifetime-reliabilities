function [ robustness_index ] = system_robustness_given_p_component_2( p_component,beta_component,input )
%find the Robustness index given p_component
%   input
%       p : design parameter
%Initialize variables 
robustness_index=zeros(1,input.N_sim);

    parfor i_sim=1:input.N_sim
        disp(i_sim/input.N_sim)
        beta_system = system_beta_given_p_component( p_component(i_sim),input,i_sim );
        % for series or 1-component systes the failure of the component leads to the consequences H_sys
        if  or(input.sys_type(i_sim)==0, input.num_el(i_sim)==1)
            robustness_index(i_sim)=-999;
        else
            %parallel system:
            %         failure of one component but survival of system --> (C0/n+CI*p+Cdir)Pfcomp*Pfcomp given system survival
            %         failure of one component and failure of the system --> (C0+n*CI*p+Hsys)*Pfsys given component failure
            % Calculation of Pfsys|compfailure
            Prob_comp_fail_sys_surv = input.Prob_comp_failure_system_survival(i_sim);
            Indirect_risk=input.mHsys(i_sim).*normcdf(-beta_system);
            Direct_risk=input.mCdir(i_sim).*Prob_comp_fail_sys_surv;
            robustness_index(i_sim)=Direct_risk/(Direct_risk+Indirect_risk);
        end
    end
end
    


