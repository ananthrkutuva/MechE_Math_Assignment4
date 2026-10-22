%Butcher tableaus for the three explicit RK methods we're comparing
%(tableaus are from the wikipedia list of runge-kutta methods)
%picked one 2nd, one 3rd and one 4th order method so the global truncation
%error should scale like h^2, h^3 and h^4
%INPUTS:
%   None
%OUTPUTS:
%BT_list: cell array of Butcher tableau structs (one per method)
%   BT_list{i}.A: matrix of a_{ij} values
%   BT_list{i}.B: vector of b_i values
%   BT_list{i}.C: vector of c_i values
%method_names: name of each method (same order as BT_list)
function [BT_list,method_names] = define_RK_methods()
    %Heun's method (2nd order, 2 stages)
    %forward euler guess for the end of the step, then average the slope
    %at the start and the slope at that guess
    Heun = struct();
    Heun.C = [0, 1];
    Heun.B = [1/2, 1/2];
    Heun.A = [0, 0;...
              1, 0];

    %Kutta's third order method (3rd order, 3 stages)
    Kutta3 = struct();
    Kutta3.C = [0, 1/2, 1];
    Kutta3.B = [1/6, 2/3, 1/6];
    Kutta3.A = [  0, 0, 0;...
                1/2, 0, 0;...
                 -1, 2, 0];

    %classic RK4 (4th order, 4 stages)
    %slope at the start, two at the midpoint, one at the end
    %weighted 1-2-2-1 (same weights as simpson's rule)
    RK4 = struct();
    RK4.C = [0, 1/2, 1/2, 1];
    RK4.B = [1/6, 1/3, 1/3, 1/6];
    RK4.A = [  0,   0, 0, 0;...
             1/2,   0, 0, 0;...
               0, 1/2, 0, 0;...
               0,   0, 1, 0];

    BT_list = {Heun, Kutta3, RK4};
    method_names = {'Heun''s Method','Kutta''s Third Order','Classic RK4'};
end
