%% Further Prodata Ferum input
global probdata gfundata 
for i=1:input.N_sim
    % Type of joint distribution: Nataf
    probdata(i).transf_type = 3;
    % Numerical computation of the modified Nataf correlation matrix
    probdata(i).Ro_method   = 1;
    % Flag for computation of sensitivities w.r.t. means, standard deviations,
    % parameters and correlation coefficients: 1 = yes 0 = no
    probdata(i).flag_sens   = 0;
end
%% Further LSF Ferum input
for i=1:3
    gfundata(i).evaluator  = 'basic';
    gfundata(i).type       = 'expression';   % Do no change this field!

    % Flag for computation of sensitivities w.r.t. thetag parameters of the limit-state function
    % 1: all sensitivities assessed, 0: no sensitivities assessment
    gfundata(i).flag_sens  = 0;
end
%% Further Ferum options

% Ferum analysis options
analysisopt.analysistype = 10;               %FORM analysis
analysisopt.echo_flag = 0;                   % 1: FERUM interactive mode, 0: FERUM silent mode
analysisopt.multi_proc = 1;                  % 1: block_size g-calls sent simultaneously
                                             %    - gfunbasic.m is used and a vectorized version of gfundata.expression is available.
                                             %      The number of g-calls sent simultaneously (block_size) depends on the memory
                                             %      available on the computer running FERUM.
                                             %    - gfunxxx.m user-specific g-function is used and able to handle block_size computations
                                             %      sent simultaneously, on a cluster of PCs or any other multiprocessor computer platform.
                                             % 0: g-calls sent sequentially
%analysisopt.block_size           = 500000;   % Number of g-calls to be sent simultaneously

% FORM analysis options
analysisopt.ig_max               = 1000;      % Maximum number of iterations allowed in the search algorithm
analysisopt.il_max               = 10;        % Maximum number of line iterations allowed in the search algorithm
analysisopt.e1                   = 1e-2;     % Tolerance on how close design point is to limit-state surface
analysisopt.e2                   = 1e-2 ;    % Tolerance on how accurately the gradient points towards the origin
analysisopt.step_code            = 0;        % 0: step size by Armijo rule, otherwise: given value is the step size
analysisopt.Recorded_u           = 0;        % 0: u-vector not recorded at all iterations, 1: u-vector recorded at all iterations
analysisopt.Recorded_x           = 0;        % 0: x-vector not recorded at all iterations, 1: x-vector recorded at all iterations
analysisopt.grad_flag            = 'ddm';    % 'ddm': direct differentiation, 'ffd': forward finite difference
analysisopt.ffdpara              = 1000;     % Parameter for computation of FFD estimates of gradients - Perturbation = stdv/analysisopt.ffdpara;
                                             % Recommended values: 1000 for basic limit-state functions, 50 for FE-based limit-state functions
analysisopt.ffdpara_thetag       = 1000;     % Parameter for computation of FFD estimates of dbeta_dthetag
                                             % perturbation = thetag/analysisopt.ffdpara_thetag if thetag ~= 0 or 1/analysisopt.ffdpara_thetag if thetag == 0;
                                             % Recommended values: 1000 for basic limit-state functions, 100 for FE-based limit-state functions
                                             
femodel = [];
randomfield = []; 