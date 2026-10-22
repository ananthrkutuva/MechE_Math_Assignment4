%Runs numerical integration arbitrary RK method
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
%              have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%BT_struct: a struct that contains the Butcher tableau
%   BT_struct.A: matrix of a_{ij} values
%   BT_struct.B: vector of b_i values
%   BT_struct.C: vector of c_i values
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in,tspan,X0,h_ref,BT_struct)
    %same as the fixed step integrators from assignment 3, just calls
    %explicit_RK_step with whatever butcher tableau gets passed in
    %smallest number of steps where the actual step size is still <= h_ref
    num_steps = ceil((tspan(2)-tspan(1))/h_ref);
    h_avg = (tspan(2)-tspan(1))/num_steps;

    t_list = linspace(tspan(1),tspan(2),num_steps+1)';

    %each row is X at one time step (X is a column so it gets transposed)
    X_list = zeros(num_steps+1,length(X0));
    X_list(1,:) = X0';

    num_evals = 0;
    XA = X0;
    for n = 1:num_steps
        [XB,step_evals] = explicit_RK_step(rate_func_in,t_list(n),XA,h_avg,BT_struct);
        num_evals = num_evals + step_evals;

        X_list(n+1,:) = XB';
        XA = XB;
    end
end
