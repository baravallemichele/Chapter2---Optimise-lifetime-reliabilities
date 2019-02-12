function [ p_component ] = design_component( beta_component_req,input )
%the function finds the p that gives exactly beta_component=beta_component_req
    %Gradient based optimization options
    options_fminunc = optimoptions('fminunc','Algorithm','quasi-newton','Display','Notify','StepTolerance',1e-4);
    parfor i_sim=1:input.N_sim
       % Function: squared differences
       myfun=@(pp_component) (component_beta_given_p_component( pp_component,input,i_sim ) - beta_component_req).^2; %pp_component = dummy variable
       %estimated design as starting point for searching algorithm
       xrd=1*exp(-0.4*beta_component_req*input.COV_XR(i_sim));
       rd=1*exp(-0.7*beta_component_req*input.COV_R(i_sim));
       gsd=1; gpd=1;
       xqd=1*exp(+0.4*beta_component_req*input.COV_XQ(i_sim));
       a=pi/sqrt(6)/input.stddev_Q(i_sim); u=input.mean_Q(i_sim)-0.577/a;
       qd=u-(1/a)*log(-log(normcdf(0.7*beta_component_req)));
       p_start=1/(xrd*rd).*((1-input.aQ(i_sim))*(input.aG(i_sim)*gsd+(1-input.aG(i_sim))*gpd)+input.aQ(i_sim)*qd*xqd);%start point for search algorithm
       % Find p_component givign beta_component=beta_component_req
       [p_component(i_sim),delta_beta_at_opt]=fminunc(myfun,p_start,options_fminunc);
       % Check results
       p_startt=p_start+0.2;
       while delta_beta_at_opt>(0.05^2) 
           disp('The design of the component was not successful \beta_{comp}~=\beta_{comp_req}')
           [p_component(i_sim),delta_beta_at_opt]=fminsearch(myfun,p_startt,input.options_fminsearch);
           p_startt=p_startt+0.2;
       end
    end
end