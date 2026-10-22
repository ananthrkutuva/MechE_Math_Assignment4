close all;
% global truncation error experiment for the three RK methods
% integrates the planetary motion IVP over three orbits for 31 step sizes
% (even in log space) and compares the last point to the exact solution
% from compute_planetary_motion at the end time
% fits error = k*h^p for each method to get p

% figures for the report, saved one folder up so they don't end up in the repo
fig_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Assignment 4 Report Figures');
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% same orbit as run_RK_integration_example
orbit_params = struct();
orbit_params.m_sun = 1;
orbit_params.m_planet = 1;
orbit_params.G = 40;
V0 = [8;0;0;1.5];

rate_func = @(t,V) gravity_rate_func(t,V,orbit_params);

% period from kepler's third law (see run_RK_integration_example)
energy0 = 0.5*norm(V0(3:4))^2 - orbit_params.G*orbit_params.m_sun/norm(V0(1:2));
a_orbit = -orbit_params.G*orbit_params.m_sun/(2*energy0);
orbit_period = 2*pi*sqrt(a_orbit^3/(orbit_params.G*orbit_params.m_sun));
tspan = [0,3*orbit_period];

% exact state at the end time (scalar time in, column vector out)
V_true_end = compute_planetary_motion(tspan(2),V0,orbit_params);

% past h = 0.1 the steps are too big to get around the sun (close
% approach) and the error is just the size of the orbit, so stop there
h_ref_list = logspace(-4,-1,31);

[BT_list,method_names] = define_RK_methods();

% only fit where it's really scaling like h^p. bigger steps and the higher
% order terms still matter, and RK4 hits round off noise around 1e-11
% (that's about how accurate the secant method solution is anyway)
h_fit_max = 1e-2;
error_fit_min = 1e-10;

h_avg_list = zeros(3,length(h_ref_list));
calls_list = zeros(3,length(h_ref_list));
error_list = zeros(3,length(h_ref_list));
p_h = zeros(3,1);
k_h = zeros(3,1);

for m = 1:3
    for i = 1:length(h_ref_list)
        [~,X_list,h_avg,num_evals] = explicit_RK_fixed_step_integration(rate_func,tspan,V0,h_ref_list(i),BT_list{m});

        % last row of X_list is the final value, turned back into a column
        error_list(m,i) = norm(X_list(end,:)' - V_true_end);
        h_avg_list(m,i) = h_avg;
        calls_list(m,i) = num_evals;
    end

    fit_spots = h_avg_list(m,:) <= h_fit_max & error_list(m,:) > error_fit_min;
    [p_h(m),k_h(m)] = generate_error_fit(h_avg_list(m,fit_spots),error_list(m,fit_spots));
end

% blue/red/green, same colors as the example plot
colors = [0.00, 0.25, 0.85; 0.90, 0.15, 0.05; 0.2235, 1.0000, 0.0784];
% fit lines are all black, different dashes so you can tell which is which
line_list = {'-','--','-.'};

fig1 = figure('Color','w','Position',[100,100,1000,700]);
for m = 1:3
    loglog(h_avg_list(m,:),error_list(m,:),'o','Color',colors(m,:),'MarkerFaceColor',colors(m,:), ...
        'MarkerEdgeColor','k','MarkerSize',7,'DisplayName',[method_names{m},' data']);
    hold on;
    loglog(h_avg_list(m,:),k_h(m)*h_avg_list(m,:).^p_h(m),line_list{m},'Color','k','LineWidth',2, ...
        'DisplayName',sprintf('%s fit, $p$ = %.2f',method_names{m},p_h(m)));
end
grid on; box on;
xlabel('Step size, $h$ (-)','Interpreter','latex');
ylabel('Global truncation error, $\epsilon$ (-)','Interpreter','latex');
title('Global Truncation Error vs. Step Size for Three Runge-Kutta Methods','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','FontSize',16);
% data goes from bottom left to top right. RK4 is the steepest line so
% the bottom right corner (big h, tiny error) is the emptiest spot, the top
% left one ends up covering some of the heun points
legend('Interpreter','latex','Location','southeast','FontSize',14);
exportgraphics(fig1,fullfile(fig_dir,'RK_global_truncation_step_size_fit.png'),'Resolution',300);

% table for the report (p is dimensionless)
fprintf('\nTable: global truncation error scaling with step size, p from error = k*h^p\n');
fprintf('(t = [0,%.4f] (3 orbits), fit over h <= %.0e and error > %.0e, all values dimensionless (-))\n', ...
    tspan(2),h_fit_max,error_fit_min);
fprintf('%-22s %10s %16s\n','Method','Stages','Estimate p (-)');
for m = 1:3
    fprintf('%-22s %10d %16.4f\n',method_names{m},length(BT_list{m}.B),p_h(m));
end

% also saved as a csv so it's easy to drop into the word table
global_h_table = table(method_names',p_h,'VariableNames',{'Method','p'});
writetable(global_h_table,fullfile(fig_dir,'RK_global_truncation_step_size_p_table.csv'));
