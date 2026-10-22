%This function computes the value of X at the next time step
%for any arbitrary RK method
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
%              have the form: dXdt = rate_func_in(t,X) (t is before X)
%t: the value of time at the current step
%XA: the value of X(t)
%h: the time increment for a single step i.e. delta_t = t_{n+1}- t_{n}
%BT_struct: a struct that contains the Butcher tableau
%   BT_struct.A: matrix of a_{ij} values
%   BT_struct.B: vector of b_i values
%   BT_struct.C: vector of c_i values
%OUTPUTS:
%XB: the approximate value for X(t+h) (the next step)
% formula depends on the integration method used
% num_evals: A count of the number of times that you called
%            rate_func_in when computing the next step
function [XB, num_evals] = explicit_RK_step(rate_func_in,t,XA,h,BT_struct)
    % initialize the m x s matrix of zeros (K)
    m = length(XA);
    s = length(BT_struct.B);
    K = zeros(m, s);

    %fill in K one column at a time, column i is k_i
    %A is lower triangular for an explicit method, so row i of A only has
    %weights for k_1 to k_{i-1}. the columns of K that haven't been
    %filled in yet are still zero so they don't add anything to the sum
    for i = 1:s
        %sum of a_ij*k_j for j = 1 to i-1 (eqn 5), as a dot product
        sum_val1 = K*(BT_struct.A(i,:)');

        %slope at the ith intermediate guess
        K(:,i) = rate_func_in(t+BT_struct.C(i)*h, XA+h*sum_val1);
    end

    %weighted average of all the k's using B (eqn 6)
    sum_val2 = K*(BT_struct.B');
    XB = XA + h*sum_val2;

    %rate_func_in gets called once per stage (once per column of K)
    num_evals = s;
end
