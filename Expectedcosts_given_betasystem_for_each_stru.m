function [ Exp_cost ] = Expectedcosts_given_betasystem_for_each_stru( beta_system,input,rv,i_sim )
%for a given i_sim it gives the Exp_cost associated to a beta_system
%to be used for optimizing each single sistem

p_start=input.start_p_comp(i_sim); 
%find p_component giving beta_system(pc)=beta_system
myfun=@(pp_component) (system_beta_given_p_component( pp_component,input,i_sim ) - beta_system)^2;


[p_component,fval]=fminsearch(myfun,p_start,input.options_fminsearch);
while fval>(0.05^2) %Check if it is the minimum
   p_start=p_start-0.25;
   myfun=@(pp_component) (system_beta_given_p_component( pp_component,input,i_sim ) - beta_system)^2;
   [p_component,fval]=fminsearch(myfun,p_start,input.options_fminsearch);
end
%calculate Prob_comp_fail_sys_surv

                    if  or(input.sys_type(i_sim)==0, input.num_el(i_sim)==1) %Series or 1-element systems
                        Prob_comp_failure_system_survival=0;
                    else %Parallel systems
                    % Matrix-based System Reliability Analysis (MSR)
                    % In: System reliability and sensitivity under statistical
                    % dependence by matrix-based system reliability method.
                    % Structural Safety. Volume 31, Issue 2, March 2009, Pages 148-156    
                        %run FORM to get the alpha factors
                    [ beta_sys,alpha_system ] = system_beta_given_p_component( p_component,input,i_sim );
                    [ beta_comp,alpha_component ] = component_beta_given_p_component( p_component,input,i_sim );
                    %Correlation coefficients
                    rho_comp_comp=1-alpha_component(2)^2;%correlation between failure of two components
                    rho_comp_sys=alpha_system(1:6)*alpha_component(1:6)';%correlation between failure of one component and failure of system
                    
                    % MSR script modified for avoiding to create
                    % a input file *.m for each iteration
                    num_el=input.num_el(i_sim);
                    dum2=[];
                    dum2=[(0:num_el) 0 -(num_el+1)];
                    sys_def={dum2,'link'};
                    beta=[beta_comp.*ones(1,num_el) beta_sys]';
                    R=zeros(num_el+1,num_el+1);
                    for r=1:(num_el+1)
                        for c=r:(num_el+1)
                            if c==r
                            R(r,c)=1;
                            else if c<(num_el+1)
                                R(r,c)=rho_comp_comp;
                                else
                                R(r,c)=rho_comp_sys;
                                end
                            end
                        end
                    end
                    n_CSRV=1; integ_method='direct';
                    %%%%%%%%%% lines of the script developed by Bora
                    %%%%%%%%%% Gencturk from FERUMsystems_MSR package at http://projects.ce.berkeley.edu/ferum/Download/download.html
                           if (length(sys_def{1})==1 && (abs(sys_def{1})<=1)) || isempty(sys_def{1})==1;
                                    sys_type='component';
                                else
                                    n_comp=max(abs(sys_def{1}));
                                    [C]=event_matrix(n_comp);
                                    [c_sys,sys_type,n_cut_sets]=sys_event(sys_def,C);
                                end
                                R=R+R'; R=R-diag(diag(R))+eye(length(R));
                                n=length(R);
                                [r,R_DS,res_check,iter_results]=gen_DS_solver(R,n_CSRV);
                                [Pf_direct]=failure_prob(beta,r,c_sys,sys_type,sys_def,integ_method);
                    %%%%%%%%%%%%%%
                    Prob_comp_failure_system_survival=Pf_direct;
                    end


                    if  or(input.sys_type(i_sim)==0, input.num_el(i_sim)==1)
                        Contruction_costs= (1+input.num_el(i_sim).*input.mCI(i_sim).*p_component);
                        Obsolescence_costs=(1+input.num_el(i_sim).*input.mCI(i_sim).*p_component+input.mA(i_sim)).*input.mOmega./log(1+input.mGamma_reftime);
                        Failure_costs=(1+input.num_el(i_sim).*input.mCI(i_sim).*p_component+input.mHsys(i_sim)).*normcdf(-beta_system)./log(1+input.mGamma_reftime);
                        Exp_cost=Contruction_costs+Obsolescence_costs+Failure_costs;
                    else
                        Contruction_costs= (1+input.num_el(i_sim).*input.mCI(i_sim).*p_component);
                        Obsolescence_costs=(1+input.num_el(i_sim).*input.mCI(i_sim).*p_component+input.mA(i_sim)).*input.mOmega./log(1+input.mGamma_reftime);
                        System_Failure_costs=(1+input.num_el(i_sim).*input.mCI(i_sim).*p_component+input.mHsys(i_sim)).*normcdf(-beta_system)./log(1+input.mGamma_reftime);
                        Component_Failure_System_Survival_costs=(1./input.num_el(i_sim)+input.mCI(i_sim).*p_component+input.mCdir(i_sim)).*Prob_comp_failure_system_survival./log(1+input.mGamma_reftime);
                        Exp_cost=Contruction_costs+Obsolescence_costs+System_Failure_costs+Component_Failure_System_Survival_costs;
                    end

end
    


