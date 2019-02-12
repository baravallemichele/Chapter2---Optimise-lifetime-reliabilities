global gfundata
%% One variable load only
% limit state function: p_component*XR*R-(1-aQ)*[aG*GS+(1-aG)*GP]-aQ*XQ*Q;
%Random variables: 	x1:XR - model uncertainties 
%					x2:R  - resistance 
%					x3:GS - self-weight
%					x4:GP - permanent load
%					x5:XQ - time-independent part of the load
%					x6:Q  - Load maxima in the reference period selected
%Parameters:        1 - p_component
%                   2 - aQ
%                   3 - aG
gfundata(1).parameter = 'yes'; 
    gfundata(1).expression = 'gfundata(1).thetag(1)*x(1)*x(2)-(1-gfundata(1).thetag(3))*( gfundata(1).thetag(2)*x(3)+(1-gfundata(1).thetag(2))*x(4) )-gfundata(1).thetag(3)*x(5)*x(6)';
    % Give explicit gradient expressions with respect to the involved quantities (in the order x(1), x(2), ...) if DDM is used:
        gfundata(1).dgdq = { 'gfundata(1).thetag(1)*x(2)' ;
                             'gfundata(1).thetag(1)*x(1)' ;
                             '-(1-gfundata(1).thetag(3))*gfundata(1).thetag(2)';
                             '-(1-gfundata(1).thetag(3))*(1-gfundata(1).thetag(2))';
                             '-gfundata(1).thetag(3)*x(6)';
                             '-gfundata(1).thetag(3)*x(5)';
                             '0';
                             '0'};
    % Give explicit gradient expressions with respect to the limit-state function parameters 
    %(in the order thetag(1), thetag(2), ...) if DDM is used:
        gfundata(1).dgthetag = {'x(1)*x(2)';
                                '-(1-gfundata(1).thetag(3))*(x(3)-x(4))';
                                'gfundata(1).thetag(2)*x(3) +(1-gfundata(1).thetag(2))*x(4)-x(5)*x(6)'};
%% Rackwitz simple lsf 
% limit state function: p_component*R-Q;
%Random variables: 	x1:R -  resistance 
%					x2:Q  - Load maxima in the reference period selected
%Parameters         1 - p_component
gfundata(2).parameter = 'yes'; 
    gfundata(2).expression = 'gfundata(2).thetag(1)*x(1)-x(2)';
    % Give explicit gradient expressions with respect to the involved quantities (in the order x(1), x(2), ...) if DDM is used:
        gfundata(2).dgdq = { 'gfundata(2).thetag(1)' ;
                             '-1'};
    % Give explicit gradient expressions with respect to the limit-state function parameters 
    %(in the order thetag(1), thetag(2), ...) if DDM is used:
        gfundata(2).dgthetag = {'x(2)'};
        
%% Parallel system with ductile elements
% limit state function: g=p_component*XR*(s1*R1+s2*R2+s3*R3+..+s10*R10)-(1-aQ)*(aG*n*GS+(1-aG)*n*GP)-aQ*XQ*(n*Q)
%Random variables: 	x1:XR - model uncertainties 
%					x2:R1 - resistance element 1
%					x3:GS - self-weight
%					x4:GP - permanent load
%					x5:XQ - time-independent part of the load
%					x6:Q  - Load maxima in the reference period selected
%                   x7:R2 - resistance element 2
%                   x8:R3
%                   ...
%                   x15:R10 resistance element 10
%Parameters:        1 - p_component
%                   2 - aQ
%                   3 - aG
%                   4 - s1
%                   5 - s2
%                   ...
%                   13 - s10 
gfundata(3).parameter = 'yes'; 
    gfundata(3).expression = 'gfundata(3).thetag(1)*x(1)*(x(2)*gfundata(3).thetag(4)+x(7)*gfundata(3).thetag(5)+x(8)*gfundata(3).thetag(6)+x(9)*gfundata(3).thetag(7)+x(10)*gfundata(3).thetag(8)+x(11)*gfundata(3).thetag(9)+x(12)*gfundata(3).thetag(10)+x(13)*gfundata(3).thetag(11)+x(14)*gfundata(3).thetag(12)+x(15)*gfundata(3).thetag(13))-(1-gfundata(3).thetag(3))*( gfundata(3).thetag(2)*x(3)+(1-gfundata(3).thetag(2))*x(4) )-gfundata(3).thetag(3)*x(5)*x(6)';
    % Give explicit gradient expressions with respect to the involved quantities (in the order x(1), x(2), ...) if DDM is used:
        gfundata(3).dgdq = { 'gfundata(3).thetag(1)*(x(2)*gfundata(3).thetag(4)+x(7)*gfundata(3).thetag(5)+x(8)*gfundata(3).thetag(6)+x(9)*gfundata(3).thetag(7)+x(10)*gfundata(3).thetag(8)+x(11)*gfundata(3).thetag(9)+x(12)*gfundata(3).thetag(10)+x(13)*gfundata(3).thetag(11)+x(14)*gfundata(3).thetag(12)+x(15)*gfundata(3).thetag(13))' ;
                             'gfundata(3).thetag(1)*x(1)*gfundata(3).thetag(4)' ;
                             '-(1-gfundata(3).thetag(3))*gfundata(3).thetag(2)';
                             '-(1-gfundata(3).thetag(3))*(1-gfundata(3).thetag(2))';
                             '-gfundata(3).thetag(3)*x(6)';
                             '-gfundata(3).thetag(3)*x(5)';
                             'gfundata(3).thetag(1)*x(1)*gfundata(3).thetag(5)' ;
                             'gfundata(3).thetag(1)*x(1)*gfundata(3).thetag(6)' ;
                             'gfundata(3).thetag(1)*x(1)*gfundata(3).thetag(7)' ;
                             'gfundata(3).thetag(1)*x(1)*gfundata(3).thetag(8)' ;
                             'gfundata(3).thetag(1)*x(1)*gfundata(3).thetag(9)' ;
                             'gfundata(3).thetag(1)*x(1)*gfundata(3).thetag(10)' ;
                             'gfundata(3).thetag(1)*x(1)*gfundata(3).thetag(11)' ;
                             'gfundata(3).thetag(1)*x(1)*gfundata(3).thetag(12)' ;
                             'gfundata(3).thetag(1)*x(1)*gfundata(3).thetag(13)' };
    % Give explicit gradient expressions with respect to the limit-state function parameters 
    %(in the order thetag(1), thetag(2), ...) if DDM is used:
        gfundata(3).dgthetag = {'x(1)*(x(2)*gfundata(3).thetag(4)+x(7)*gfundata(3).thetag(5)+x(8)*gfundata(3).thetag(6)+x(9)*gfundata(3).thetag(7)+x(10)*gfundata(3).thetag(8)+x(11)*gfundata(3).thetag(9)+x(12)*gfundata(3).thetag(10)+x(13)*gfundata(3).thetag(11)+x(14)*gfundata(3).thetag(12)+x(15)*gfundata(3).thetag(13))';
                                '-(1-gfundata(3).thetag(3))*(x(3)-x(4))';
                                'gfundata(3).thetag(2)*x(3) +(1-gfundata(3).thetag(2))*x(4)-x(5)*x(6)';
                                'gfundata(3).thetag(1)*x(1)*x(2)';
                                'gfundata(3).thetag(1)*x(1)*x(7)';
                                'gfundata(3).thetag(1)*x(1)*x(8)';
                                'gfundata(3).thetag(1)*x(1)*x(9)';
                                'gfundata(3).thetag(1)*x(1)*x(10)';
                                'gfundata(3).thetag(1)*x(1)*x(11)';
                                'gfundata(3).thetag(1)*x(1)*x(12)';
                                'gfundata(3).thetag(1)*x(1)*x(13)';
                                'gfundata(3).thetag(1)*x(1)*x(14)';
                                'gfundata(3).thetag(1)*x(1)*x(15)' };
                            
                            
                            