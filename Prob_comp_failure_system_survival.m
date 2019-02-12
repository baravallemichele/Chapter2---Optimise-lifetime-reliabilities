function [ Prob_comp_failure_system_survival ] = Prob_comp_failure_system_survival( p_component,input,rv )
%find the Pf of the events of system surviving with failed components
%   input
%       p_component : design parameter of the component
% Initialize varaibles
    Prob_comp_failure_system_survival=zeros(input.N_sim,1);
% Calcualte Prob_comp_failure_system_survival
    parfor i_sim=1:input.N_sim
        if  or(input.sys_type(i_sim)==0, input.num_el(i_sim)==1) %Series or 1-element systems
             Prob_comp_failure_system_survival(i_sim)=0;
        else %Parallel systems
            % Matrix-based System Reliability Analysis (MSR)
            % In: System reliability and sensitivity under statistical
            % dependence by matrix-based system reliability method.
            % Structural Safety. Volume 31, Issue 2, March 2009, Pages 148-156
                % Run FORM to get the alpha factors for calcualting
                % correlation beteen failure in different ref_time periods
                    [ beta_sys,alpha_system ] = system_beta_given_p_component( p_component(i_sim),input,i_sim );
                    [ beta_comp,alpha_component ] = component_beta_given_p_component( p_component(i_sim),input,i_sim );
                %Correlation coefficients
                    rho_comp_comp=1-alpha_component(2)^2;%correlation between failure of two components (Only the resistance x(2) is independent from component to component)
                    rho_comp_sys=alpha_system(1:6)*alpha_component(1:6)';%correlation between failure of one component and failure of system
                % MSR script modified for avoiding to create
                % a input file *.m for each iteration
                    num_el=input.num_el(i_sim);
                    dum2=[];
                    dum2=[(0:num_el) 0 -(num_el+1)];
                    sys_def={dum2,'link'};
                    beta=[beta_comp.*ones(1,num_el) beta_sys]';
                    R=zeros(num_el+1,num_el+1); %Correlation matrix
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
                Prob_comp_failure_system_survival(i_sim)=Pf_direct;
            % end MSR method
        end %End if
    end %End parfor
end


