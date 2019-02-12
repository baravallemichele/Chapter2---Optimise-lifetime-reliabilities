function [ p_component ] = design_system( beta_system_req,input )
%the function finds the p_component that guives exactly beta_system=beta_system_req

%initialize varaibles
    p_component=zeros(1,input.N_sim);
%Options of the search algorithms    
    options_fminsearch = optimset('Display','Notify','TolFun',1e-4);
    options_fminunc = optimoptions('fminunc','Algorithm','quasi-newton','Display','Notify','StepTolerance',1e-4);
    %simulations
    parfor i_sim=1:input.N_sim
           p_start=input.start_p_comp(i_sim); %Start p_component value
           %pp_component = dummy variable
           myfun=@(pp_component) (system_beta_given_p_component( pp_component,input,i_sim ) - beta_system_req).^2;
           %[p_out]=fminsearch(myfun,p_start,options_fminsearch);
           p_component(i_sim)=fminunc(myfun,p_start,options_fminunc);
           %Check
           p_startt=p_start;
           while system_beta_given_p_component( p_component(i_sim),input,i_sim )>beta_system_req+.05%check
               p_startt=p_startt-0.5;
               disp('The design of the system was not successful \beta_{sys}~=\beta_{sys_req}')
               disp(system_beta_given_p_component( p_component(i_sim),input,i_sim ))
               disp(beta_system_req)
               p_component(i_sim)=fminsearch(myfun,p_startt,options_fminsearch);
           end
           while system_beta_given_p_component( p_component(i_sim),input,i_sim )<beta_system_req-.05%check
               p_startt=p_startt+0.5;
               disp('The design of the system was not successful \beta_{sys}~=\beta_{sys_req}')
               disp(system_beta_given_p_component( p_component(i_sim),input,i_sim ))
               disp(beta_system_req)
               p_component(i_sim)=fminsearch(myfun,p_startt,options_fminsearch);
           end
    end
end