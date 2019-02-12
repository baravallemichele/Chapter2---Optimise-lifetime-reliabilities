function [ beta_comp,alpha_component ] = component_beta_given_p_component( p_component,input,i_sim )
%find the reliability of the component given a design of the component: p_component
%   input
%       p_component : design parameter of the component
%		limit state function: p_component*XR*R-(1-aQ)*[aG*GS+(1-aG)*GP]-aQ*XQ*Q;
    %Define global varaibles
    global probdata gfundata
    %Load FERUM options
    LoadFerumOptions;
    %Load limit state functions
    limit_state_functions;
    for ip=1:length(p_component)
           %Random vars    (        distr.     mean                 std                     start point
            probdata(i_sim).marg = [2,         1,                   input.COV_XR(i_sim),    1-input.COV_XR(i_sim),          nan, nan, nan, nan, 0;...
                                    2,         1,                   input.COV_R(i_sim),     1-1.5*input.COV_R(i_sim),       nan, nan, nan, nan, 0;...
                                    1,         1,                   input.COV_GS(i_sim),    1,                              nan, nan, nan, nan, 0;...
                                    1,         1,                   input.COV_GP(i_sim),    1,                              nan, nan, nan, nan, 0;...
                                    2,         1,                   input.COV_XQ(i_sim),    1+input.COV_XQ(i_sim),          nan, nan, nan, nan, 0;...
                                    15,        input.mean_Q(i_sim), input.stddev_Q(i_sim),  input.mean_Q(i_sim)+2*input.stddev_Q(i_sim),   nan, nan, nan, nan, 0;...
                                    1,         1,                   999,                    1,                              nan, nan, nan, nan, 0;...
                                    1,         1,                   999,                    1,                              nan, nan, nan, nan, 0];
            %Random variables: 	x1:XR - model uncertainties 
            %					x2:R  - resistance 
            %					x3:GS - self-weight
            %					x4:GP - permanent load
            %					x5:XQ - time-independent part of the load
            %					x6:Q  - Load maxima in the reference period selected
            % Correlation matrix (square matrix with dimension equal to number of r.v.'s)
            probdata(i_sim).correlation = eye(size(probdata(i_sim).marg,1));
            % Determine the parameters,the mean and standard deviation associated with the distribution of each random variable
            probdata(i_sim).parameter = distribution_parameter(probdata(i_sim).marg);   
            % define the parameter values
            gfundata(1).thetag = [p_component(ip) input.aG(i_sim) input.aQ(i_sim)];
            % FORM
            [formresults] = form(1,probdata(i_sim),analysisopt,gfundata,femodel,randomfield);
            beta_comp(ip) = formresults.beta1;
            alpha_component(ip,:)=formresults.alpha;
    end
end

