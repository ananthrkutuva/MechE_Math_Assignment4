close all;
% runs each of the three RK methods on the planetary motion IVP with the
% same step size and plots them against the exact solution from
% compute_planetary_motion (orion's secant method version)
% two subplots, x(t) on top and y(t) on the bottom, for three orbits

% figures for the report, saved one folder up so they don't end up in the repo
fig_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Assignment 4 Report Figures');
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% same system as the usage example in compute_planetary_motion
% starts at the far point of the orbit (moving slower than a circular orbit
% would) so the orbit comes out as an ellipse
orbit_params = struct();
orbit_params.m_sun = 1;
orbit_params.m_planet = 1;
orbit_params.G = 40;
V0 = [8;0;0;1.5];

% wrap the rate function so it's only a function of t and V
rate_func = @(t,V) gravity_rate_func(t,V,orbit_params);

% period of the orbit from kepler's third law, T = 2*pi*sqrt(a^3/(G*m_sun))
% the semi major axis a comes from the energy per unit mass of the orbit,
% v^2/2 - G*m_sun/r = -G*m_sun/(2a), which stays constant the whole time
energy0 = 0.5*norm(V0(3:4))^2 - orbit_params.G*orbit_params.m_sun/norm(V0(1:2));
a_orbit = -orbit_params.G*orbit_params.m_sun/(2*energy0);
orbit_period = 2*pi*sqrt(a_orbit^3/(orbit_params.G*orbit_params.m_sun));
fprintf('orbital period T = %.4f\n',orbit_period);

% three full orbits
tspan = [0,3*orbit_period];

% one step size for all three methods. small enough that heun's method
% (the worst one) still lands on top of the exact solution after 3 orbits
h_ref = 0.02;

[BT_list,method_names] = define_RK_methods();

% exact solution on a fine grid so the black line is smooth
t_fine = linspace(tspan(1),tspan(2),3000);
V_exact = compute_planetary_motion(t_fine,V0,orbit_params);

% blue/red/green, same colors as the truncation plot
colors = [0.00, 0.25, 0.85; 0.90, 0.15, 0.05; 0.2235, 1.0000, 0.0784];
marker_list = {'o','s','d'};
size_list = [9,8,8];

% there are 1700+ steps, so only every 25th step gets a marker (otherwise
% it's just a solid blob). each method starts at a different step so the
% markers don't land on top of each other
marker_spacing = 25;
marker_offset = [1,9,17];

% labels for the two subplots (column 1 of V is x, column 2 is y)
comp_names = {'x','y'};

fig = figure('Color','w','Position',[100,100,1100,900]);
for c = 1:2
    subplot(2,1,c);
    hold on; box on; grid on;

    true_plot = plot(t_fine,V_exact(:,c),'k-','LineWidth',2,'DisplayName','Exact solution');

    for m = 1:3
        [t_list,X_list,h_avg,num_evals] = explicit_RK_fixed_step_integration(rate_func,tspan,V0,h_ref,BT_list{m});
        idx = marker_offset(m):marker_spacing:length(t_list);
        num_plots(m) = plot(t_list(idx),X_list(idx,c),marker_list{m},'Color',colors(m,:), ...
            'MarkerFaceColor',colors(m,:),'MarkerEdgeColor','k','MarkerSize',size_list(m), ...
            'DisplayName',sprintf('%s, $h$ = %.4f',method_names{m},h_avg));

        % how far off each method is by the end, just to check
        if c == 1
            V_true = compute_planetary_motion(t_list,V0,orbit_params);
            fprintf('%s h = %.4f: max position error = %.3e, %d rate function calls\n', ...
                method_names{m},h_avg,max(vecnorm(X_list(:,1:2)-V_true(:,1:2),2,2)),num_evals);
        end
    end

    % the orbit spends most of its time at the far end (x near 8), so the
    % top of the x plot is full. leave headroom above the data and put the
    % legend in one row across the top so it doesn't cover anything
    xlim(tspan);
    if c == 1
        ylim([-4,13]);
    else
        ylim([-6,8.5]);
    end
    xlabel('Time, $t$ (-)','Interpreter','latex');
    ylabel(sprintf('Position, $%s(t)$ (-)',comp_names{c}),'Interpreter','latex');
    title(sprintf('$%s$-Component of the Planet''s Position',comp_names{c}),'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex','FontSize',16);
    legend([true_plot,num_plots],'Interpreter','latex','Location','north', ...
        'Orientation','horizontal','NumColumns',2,'FontSize',13);
end

sgtitle('Runge-Kutta Approximations of Planetary Motion over Three Orbits', ...
    'Interpreter','latex','FontSize',18);
exportgraphics(fig,fullfile(fig_dir,'RK_example_integration.png'),'Resolution',300);
