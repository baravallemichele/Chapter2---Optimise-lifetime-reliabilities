%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scripts behind calculations presented in the article:
% A risk-based approach for calibration of design codes
% Journal
% https://doi.org/
%     Michele Baravalle a,*
%     Jochen Köhler a
% a)Dept. of Structural Engineering, Norwegian University of Science & Technology, Richard Birkelands vei 1A, 7491 Trondheim, Norway
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2017 Michele Baravalle, NTNU
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or (at
% your option) any later version.

% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Details: 
% This scripts optimizes lifetime target relaibilitity at component and 
% system level with a risk-based approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
clf
clear global
%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Costs 
    %FIxed construction costs are equal to 1MU, all other costs are
    %expressed in multiples of C0
    %The eman values of the costs are considered unknown and represented by
    %a unifrm distribution in the following intervals
    input.lim_mCI=[1e-4 1e-3];  %limits of E[CI] slope of the construction costs
    input.lim_mHsys=[4 9];      %limits of E[Hsys] costs due to system failure
    input.lim_mCdir=[0 2];      %limits of E[Cdir] direct consequences associated to component failure in parallel system (given system survival)
    input.lim_mA=[0 1];         %limtis of E[A] costs of demolition due to obsolescence
% Other input    
    input.design_life=50;       %expected design life in years
    input.ref_time=input.design_life; %reference time for reliability calculations with the time-integrated approach
    input.mOmega=1/input.design_life*input.ref_time;  %E[Omega] expected design life in the ref_time units
    input.mGamma_1yr=0.03;      %E[yearly interest rate] yearly mean real societal interest rate
    input.mGamma_reftime=(1+input.mGamma_1yr)^input.ref_time-1; %E[interest refate referring to design life]
% Define the domain where the total Expected costs are to be estiamted
    dom_beta_min=1.5; %minimum value of relaibility index
    dom_beta_max=4.25; %maximum values of relaibility index
    num_beta_points=10; %number of points where estiamting the expected costs in [dom_beta_min,dom_beta_max]
% System configurations
    input.num_el=    [1   2   3   4   5  6  7  8  9  10]; %OBS: do not insert more than 10 elements in the systems!
    input.w_num_el=  [30  30  20  10  5  5  0  0  0  0]; %Probability mass function for num_el
    %             =  [100 0   0   0   0  0  0  0  0  0]; %for only systems with one component
    input.w_sys_type=[.5  1-.5];  %Probability mass function for Series and Parallel systems 
    %               =[1   0];     %Series only 
    %               =[0   1];     %Parallel only 
    %               =[.75 .25];   %Series 75- Parallel 25 
% Limit state fucntion parameters and random variables
    input.lim_aG=[.6 1];          %Limits of the aG factor - Uniformly distributed
    input.lim_aQ=[0.1 0.9];       %Limits of the aQ factor - Uniformly distributed
    input.lim_COV_XR=[.05 .20];   %Limits of the COV_XR - Uniformly distributed
    input.lim_COV_R=[.05 0.25];   %Limits of the COV_R - Uniformly distributed
    input.lim_COV_GS=[0.04 0.10]; %Limits of the COV_GS - Uniformly distributed
    input.lim_COV_GP=[0.05 0.20]; %Limits of the COV_GP - Uniformly distributed
    input.lim_COV_XQ=[0.10 0.30]; %Limits of the COV_XQ - Uniformly distributed
    input.lim_COV_Q_1yr=[0.10 0.90]; %Limits of the COV_Q - Uniformly distributed OBS:Yearly maxima, mean value of yearly maxima = 1!!!
% Simulation parameters
    input.N_sim=2500;             %Number of Monte Carlo simulations (M in the article)   
% Optimization parameters
    input.options_fminsearch = optimset('Display','Notify','TolFun',1e-5,'TolX',1e-3,'MaxFunEvals',1000');
 
%% SCRIPTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Latin hypercube sampling of the random variables
    X = lhsdesign(input.N_sim,13) ;   
% Sampling of "structural systems"
    input.w_num_el=input.w_num_el./sum(input.w_num_el);
    num_el=[]; %vector with number of elements in the system
    for i=1:length(input.num_el)
        num_el=[num_el input.num_el(i).*ones(1,input.N_sim.*input.w_num_el(i))]; 
    end
    input.num_el=num_el;
% Vector with system type (0=series, 1=parallel)
    dum=X(:,1); dum(dum<=input.w_sys_type(1))=0; dum(dum>input.w_sys_type(1))=1;
    input.sys_type=dum;
% Sampling of LSF and random varaibles' parameters
    input.mCI=input.lim_mCI(1)+X(:,2).*(input.lim_mCI(2)-input.lim_mCI(1));
    input.mHsys=input.lim_mHsys(1)+X(:,3).*(input.lim_mHsys(2)-input.lim_mHsys(1));
    input.mCdir=input.lim_mCdir(1)+X(:,4).*(input.lim_mCdir(2)-input.lim_mCdir(1));
    input.mA=input.lim_mA(1)+X(:,5).*(input.lim_mA(2)-input.lim_mA(1));
    input.aG=input.lim_aG(1)+X(:,6).*(input.lim_aG(2)-input.lim_aG(1));
    input.aQ=input.lim_aQ(1)+X(:,7).*(input.lim_aQ(2)-input.lim_aQ(1));
    input.COV_XR=input.lim_COV_XR(1)+X(:,8).*(input.lim_COV_XR(2)-input.lim_COV_XR(1));
    input.COV_R=input.lim_COV_R(1)+X(:,9).*(input.lim_COV_R(2)-input.lim_COV_R(1));
    input.COV_GS=input.lim_COV_GS(1)+X(:,10).*(input.lim_COV_GS(2)-input.lim_COV_GS(1));
    input.COV_GP=input.lim_COV_GP(1)+X(:,11).*(input.lim_COV_GP(2)-input.lim_COV_GP(1));
    input.COV_XQ=input.lim_COV_XQ(1)+X(:,12).*(input.lim_COV_XQ(2)-input.lim_COV_XQ(1));
    input.COV_Q_1yr=input.lim_COV_Q_1yr(1)+X(:,13).*(input.lim_COV_Q_1yr(2)-input.lim_COV_Q_1yr(1));
% Calculate mean, std dev and cov of maxima of Q over the selected reference period
    input.stddev_Q=input.COV_Q_1yr.*1; %mean of eyarly maxima is 1!
    input.mean_Q=(1+sqrt(6)/pi*input.COV_Q_1yr*log(input.ref_time));
    input.COV_Q=input.stddev_Q./input.mean_Q;
    rando=100;
% Name of case for saving files
    name_case=strcat('Ser',num2str(input.w_sys_type(1)*100),'%_Par',num2str(input.w_sys_type(2)*100),'%_Hsys=',num2str(input.lim_mHsys(1)),'-',num2str(input.lim_mHsys(2)),'_CI=',num2str(input.lim_mCI(1)),'-',num2str(input.lim_mCI(2)),'_Cdir=',num2str(input.lim_mCdir(1)),'-',num2str(input.lim_mCdir(2)));
    
%% Section 1 - Calibration for level 1 code
% Calculate total Expected costs in 10 points in [dom_beta_min,dom_beta_max]
    figure('Name',name_case)
    subplot(1,2,2)
    j=0; beta_component=zeros(1,num_beta_points); Expected_costs_of_system=zeros(1,num_beta_points);
    for beta_component_req=linspace(dom_beta_min,dom_beta_max,num_beta_points)
       j=j+1;
       beta_component(j)=beta_component_req; %Beta at component level referring to input.ref_time 
       %Find the component design giving the system reliability requirement
            p_component=design_component(beta_component_req,input);
       %Find the probability of system survival with component failure for parallel systems     
            input.Prob_comp_failure_system_survival=Prob_comp_failure_system_survival( p_component,input,rando );
       %Calculate expected costs     
            Expected_costs_of_system(j) = system_Expected_costs_given_p_component_2( p_component,beta_component_req,input );
       start_p_comp(j,:)=p_component;%save p_components to use as starting point for design of the systems later
    end
    plot(beta_component,Expected_costs_of_system,'og')
    hold on
% Find spline interpolant  
    x=dom_beta_min:0.001:dom_beta_max;
    interp_ExpCost_values = spline(beta_component,Expected_costs_of_system,x);
    plot(x,interp_ExpCost_values,'--b');
% Find minimum of the spline interpolant    
    ExpCost_opt=min(interp_ExpCost_values); beta_opt_level1=x(interp_ExpCost_values==ExpCost_opt)
    text(beta_opt_level1-0.5,ExpCost_opt,strcat('\beta_{c,opt}^{(T_l_i_f_e )} = ',num2str(0.01*round(100*beta_opt_level1))),'Color','blue')
    xlabel('\beta_{c}^{(T_l_i_f_e )}'); ylabel('E_{\Delta,\Theta}[C_t_o_t]');grid on; title('Level 1 Code')
%% Section 2 - Calibration for level 3 code
%for each given beta design the system wi
    j=0; beta_system=zeros(1,num_beta_points); Expected_costs_of_system=zeros(1,num_beta_points);
subplot(1,2,1)
    for beta_system_req=linspace(dom_beta_min,dom_beta_max,num_beta_points)
       j=j+1;
       beta_system(j)=beta_system_req; %Beta at system level referring to input.ref_time 
       input.start_p_comp=start_p_comp(j,:); %Starting point of the searching algorithm
       %Find the component design giving the system reliability requirement
            p_component=design_system( beta_system_req,input );
       %Find the probability of system survival with component failure for parallel systems     
            input.Prob_comp_failure_system_survival=Prob_comp_failure_system_survival(p_component,input ,rando); 
       %input.Prob_comp_failure_system_survival=zeros(1,input.N_sim); %Assume zero probability of component failure and system survival.
            Exp_costs= system_Expected_costs_given_p_component( p_component,beta_system_req,input );
            Expected_costs_of_system(j) = mean(Exp_costs);
       %Save variables for checking convergence
            convergence(j).p_component=p_component;
            convergence(j).beta_system=beta_system_req;
            convergence(j).Exp_costs=Exp_costs;
    end
    plot(beta_system,Expected_costs_of_system,'og')
    hold on
% Find spline interpolant      
    x=dom_beta_min:0.001:dom_beta_max;
    interp_ExpCost_values = spline(beta_system,Expected_costs_of_system,x);
    plot(x,interp_ExpCost_values,'--b');
% Find minimum of the spline interpolant     
    ExpCost_opt=min(interp_ExpCost_values); beta_opt_level3=x(interp_ExpCost_values==ExpCost_opt)
    text(beta_opt_level3-0.5,ExpCost_opt,strcat('\beta_{sys,opt}^{(T_l_i_f_e )} = ',num2str(0.01*round(100*beta_opt_level3))),'Color','blue')
    xlabel('\beta_{sys}^{(T_l_i_f_e )}'); ylabel('E_{\Delta,\Theta}[C_t_o_t]');grid on; title('Level 3 Code')
%% Section 3 - Check robustness index and relaibility of structures designed with  beta_opt
% Check robustness index of structures designed with  beta_opt_level1
    p_component_with_beta_opt_level1=design_component(beta_opt_level1,input);
    input.Prob_comp_failure_system_survival=Prob_comp_failure_system_survival( p_component_with_beta_opt_level1,input,rando );
    robustness_index_with_beta_opt_level1 = system_robustness_given_p_component_2( p_component_with_beta_opt_level1,beta_opt_level1,input );
    parfor i_sim=1:input.N_sim
        beta_sys_with_beta_opt_level1(i_sim) = system_beta_given_p_component( p_component_with_beta_opt_level1(i_sim),input,i_sim )
    end
% Check robustness index of structures designed with  beta_opt_level3
    p_component_with_beta_opt_level3=design_system(beta_opt_level3,input);
    input.Prob_comp_failure_system_survival=Prob_comp_failure_system_survival( p_component_with_beta_opt_level3,input,rando );
    robustness_index_with_beta_opt_level3 = system_robustness_given_p_component( p_component_with_beta_opt_level3,beta_opt_level3,input );
    parfor i_sim=1:input.N_sim
        beta_comp_with_beta_opt_level3(i_sim) = component_beta_given_p_component( p_component_with_beta_opt_level3(i_sim),input,i_sim )
    end
% plots Robustness
    figure('Name',name_case)
    subplot(2,1,2)
    histogram(robustness_index_with_beta_opt_level1~=-999)
    grid on; xlabel('I_r_o_b');ylabel('f_{I_r_o_b}'); title('Design with \beta_{c,opt}^{(T_l_i_f_e )}')
    subplot(2,1,1)
    histogram(robustness_index_with_beta_opt_level3~=-999)
    grid on; xlabel('I_r_o_b');ylabel('f_{I_r_o_b}'); title('Design with \beta_{sys,opt}^{(T_l_i_f_e )}')
% plots beta 
    figure('Name',name_case)
    subplot(2,1,2)
    histogram(beta_sys_with_beta_opt_level1)
    grid on; xlabel('\beta_{sys}^{(T_l_i_f_e )}');ylabel('f_{B_s_y_s}'); title('Design with \beta_{c,opt}^{(T_l_i_f_e )}')
    subplot(2,1,1)
    histogram(beta_comp_with_beta_opt_level3)
    grid on; xlabel('\beta_c^{(T_l_i_f_e )}');ylabel('f_{B_c}'); title('Design with \beta_{sys,opt}^{(T_l_i_f_e )}')
%% Section 4 - Optimize each single structure
    figure('Name',name_case)
    i=1;
% Initialize variables    
    beta_system_opt=zeros(1,input.N_sim);           %Optimal system beta
    Exp_cost_given_beta_sys=zeros(1,input.N_sim);   %Exp costs given ebta sys    
    beta_component_opt=zeros(1,input.N_sim);        %Optimal beta component
    Exp_cost_given_beta_comp=zeros(1,input.N_sim);  %Exp costs given beta component
    beta_start=beta_system(round(0.5*num_beta_points)); %Starting beta for search algorithm
    input.start_p_comp=start_p_comp(find(beta_system==beta_start),:); %start_p for search algorith 
    parfor i_sim=1:input.N_sim
        myfun=@(beta_sys) Expectedcosts_given_betasystem_for_each_stru(beta_sys,input,rando,i_sim);
        [beta_system_opt(i_sim),Exp_cost_given_beta_sys(i_sim)]=fminsearch(myfun,beta_start,input.options_fminsearch);
        myfun=@(beta_c) Expectedcosts_given_betacomponent_for_each_stru(beta_c,input,rando,i_sim);
        [beta_component_opt(i_sim),Exp_cost_given_beta_comp(i_sim)]=fminsearch(myfun,beta_start,input.options_fminsearch);
    end
% Plot
subplot(2,2,1)
histfit(beta_system_opt)%,'Normalization','pdf')
xlabel('\beta_{sys,opt}')
ylabel('f_{B_{sys,opt}} (\beta_{sys,opt})')
subplot(2,2,2)
histfit(beta_component_opt)%,'Normalization','pdf')
xlabel('\beta_{c,opt}')
ylabel('f_{B_{c,opt}} (\beta_{c,opt})')
subplot(2,2,3)
plot(beta_system_opt,Exp_cost_given_beta_sys,'o')
xlabel('\beta_{sys,opt}')
ylabel('E[C(\beta_{sys,opt})]')
subplot(2,2,4)
plot(beta_component_opt,Exp_cost_given_beta_comp,'o')
xlabel('\beta_{c,opt}')
ylabel('E[C(\beta_{c,opt})]') 
savefig(strcat(name_case,'- Optimiz each stru.fig') )

%% Section 9 - Check convergence at beta_sys=4
figure('Name',name_case)
     j=round(0.5*num_beta_points); %convergence at beta_sys=3.5
     Check_convergence_for_beta_sys_4;
