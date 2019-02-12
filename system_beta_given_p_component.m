function [ beta_sys,alpha_system ] = system_beta_given_p_component( p_component,input,i_sim )
%find the reliability of the SYSTEM given a design p of a COMPONENT
%   input
%       p_component : design parameter of the component
% Define global varaibles
    global probdata gfundata
% Initialize variables
    beta_sys=zeros(1,length(p_component));
% Load limit state function and FERUM options
    LoadFerumOptions;
    limit_state_functions;
    for ip=1:length(p_component)
        if  p_component(ip)>0
            if or(input.sys_type(i_sim)==0,input.num_el(i_sim)==1)%series systems or one-component systems
                % Calculate the component reliability
                    %Random vars    (       distr.     mean                 std                     start point
                    probdata(i_sim).marg =[ 2,         1,                   input.COV_XR(i_sim),    1*exp(-0.6*4*input.COV_XR(i_sim)),  nan, nan, nan, nan, 0;...
                                            2,         1,                   input.COV_R(i_sim),     1*exp(-0.6*4*input.COV_R(i_sim)),   nan, nan, nan, nan, 0;...
                                            1,         1,                   input.COV_GS(i_sim),    1,                                  nan, nan, nan, nan, 0;...
                                            1,         1,                   input.COV_GP(i_sim),    1,                                  nan, nan, nan, nan, 0;...
                                            2,         1,                   input.COV_XQ(i_sim),    1*exp(0.3*4*input.COV_XQ(i_sim)),   nan, nan, nan, nan, 0;...
                                            15,        input.mean_Q(i_sim), input.stddev_Q(i_sim),  1.5*(input.mean_Q(i_sim)+2*input.stddev_Q(i_sim)),   nan, nan, nan, nan, 0;...
                                            1,         1,                   999,                    1,                                  nan, nan, nan, nan, 0;...
                                            1,         1,                   999,                    1,                                  nan, nan, nan, nan, 0];

                    % Correlation matrix (square matrix with dimension equal to number of r.v.'s)
                    probdata(i_sim).correlation = eye(size(probdata(i_sim).marg,1));
                    % Determine the parameters,the mean and standard deviation associated with the distribution of each random variable
                    probdata(i_sim).parameter = distribution_parameter(probdata(i_sim).marg);
                    % define the parameter values
                    gfundata(1).thetag = [p_component(ip) input.aG(i_sim) input.aQ(i_sim)];
                    % FORM
                    [formresults] = form(1,probdata(i_sim),analysisopt,gfundata,femodel,randomfield);
                    results.form = formresults;
                    % output structure
                    beta1_component=results.form.beta1;
                    beta_sys(ip)=beta1_component;
                    alpha_system(ip,:)=formresults.alpha;
                % End component reliability
                
                % Calculate series system reliability based on component reliability
                    if input.num_el(i_sim)>1 %Series systems with 2 or more components
                        alphaR_component=results.form.alpha(2);
                        %Use exact formula for series with equal components and equal
                        %correlation from Roscoe et al., Structural Safety 2015.
                        rho=(1-alphaR_component^2); %correlation between components
                        b=beta1_component; %component reliability index
                        myfun=@(t) (1-(1-normcdf(-(b-t.*sqrt(rho))./(sqrt(1-rho)))).^input.num_el(i_sim)).*normpdf(t);
                        Pf_sys=integral(myfun,-Inf,Inf,'RelTol',1e-15,'AbsTol',1e-15);
                        beta_sys(ip)=-norminv(Pf_sys);
                    end % series systems with 2 or more components
                
            else %ductile parallel systems with 2 or more components
                %Calculate system reliability
                n=input.num_el(i_sim); %number of elements
                %Limit state function:
                %g=p*xr*(x2+x7+x8+..+x15)-(1-aQ)*(aG*n*gs+(1-aG)*n*gp)-aQ*xq*(n*q)
                %Random vars    (       distr.     mean                     std                       start point
                probdata(i_sim).marg =[ 2,         1,                       input.COV_XR(i_sim),      1*exp(-0.6*4*input.COV_XR(i_sim)),  nan, nan, nan, nan, 0;...     % XR
                                        2,         1,                       input.COV_R(i_sim),       1*exp(-0.6*4*n*input.COV_R(i_sim)),   nan, nan, nan, nan, 0;...   % R element 1
                                        1,         n,                       n*input.COV_GS(i_sim),    n,                            nan, nan, nan, nan, 0;...           % GS
                                        1,         n,                       n*input.COV_GP(i_sim),    n,                            nan, nan, nan, nan, 0;...           % GP
                                        2,         1,                       input.COV_XQ(i_sim),      1*exp(0.3*4*input.COV_XQ(i_sim)),  nan, nan, nan, nan, 0;...      % XQ
                                        15,        n*input.mean_Q(i_sim),   n*input.stddev_Q(i_sim),  n*1.5*(input.mean_Q(i_sim)+2*input.stddev_Q(i_sim)),   nan, nan, nan, nan, 0;...  % Q
                                        2,         1,                       input.COV_R(i_sim),       1*exp(-0.6*4*input.COV_R(i_sim)),   nan, nan, nan, nan, 0;...     % R element 2
                                        2,         1,                       input.COV_R(i_sim),       1*exp(-0.6*4*input.COV_R(i_sim)),   nan, nan, nan, nan, 0;...     % R element 3
                                        2,         1,                       input.COV_R(i_sim),       1*exp(-0.6*4*input.COV_R(i_sim)),   nan, nan, nan, nan, 0;...     % ...
                                        2,         1,                       input.COV_R(i_sim),       1*exp(-0.6*4*input.COV_R(i_sim)),   nan, nan, nan, nan, 0;...
                                        2,         1,                       input.COV_R(i_sim),       1*exp(-0.6*4*input.COV_R(i_sim)),   nan, nan, nan, nan, 0;...
                                        2,         1,                       input.COV_R(i_sim),       1*exp(-0.6*4*input.COV_R(i_sim)),   nan, nan, nan, nan, 0;...
                                        2,         1,                       input.COV_R(i_sim),       1*exp(-0.6*4*input.COV_R(i_sim)),   nan, nan, nan, nan, 0;...
                                        2,         1,                       input.COV_R(i_sim),       1*exp(-0.6*4*input.COV_R(i_sim)),   nan, nan, nan, nan, 0;...
                                        2,         1,                       input.COV_R(i_sim),       1*exp(-0.6*4*input.COV_R(i_sim)),   nan, nan, nan, nan, 0];       % R element 10

                % Correlation matrix (square matrix with dimension equal to number of r.v.'s)
                probdata(i_sim).correlation = eye(size(probdata(i_sim).marg,1));
                % Determine the parameters,the mean and standard deviation associated with the distribution of each random variable
                probdata(i_sim).parameter = distribution_parameter(probdata(i_sim).marg);
                % vector with 1s and 0s. (e.g. [1 1 0 0 0 0 0 0 0 0]=system with 2 components)
                sys=zeros(1,10); sys(1:n)=1;
                % define the parameter values
                gfundata(3).thetag = [p_component(ip) input.aG(i_sim) input.aQ(i_sim) sys];
                % FORM
                [formresults] = form(3,probdata(i_sim),analysisopt,gfundata,femodel,randomfield);
                % output structure
                beta_sys(ip)=formresults.beta1;
                alpha_system(ip,:)=formresults.alpha;
            end
        else %case with p_component<0
            beta_sys(ip)=-999;
        end %end c
    end %and cycle on ip
end%end function


















