function neo_dyn_sierd_simulation()
    clear all;
    close all;
    clc;

    p.N_total     = 1000000;
    p.I0          = 75;
    p.days        = 180;
    p.beta        = 0.2500;
    p.gamma       = 0.0476;
    p.delta       = 0.0010;
    p.sigma       = 0.1000;
    p.f_high_risk = 0.05;
    p.I_g0        = 15;
    p.I_h0        = 60;
    p.beta_g      = 0.1600;
    p.beta_h      = 0.5600;
    p.beta_gh     = 0.0500;
    p.beta_hg     = 0.0500;
    p.gamma_g     = 0.0476;
    p.gamma_h     = 0.0476;
    p.delta_g     = 0.0008;
    p.delta_h     = 0.0030;
    p.v_g         = 0.0005;
    p.v_h         = 0.0050;
    p.efficacy    = 0.8500;

    t_span = linspace(0, p.days, p.days + 1)';
    t_real = [0; 7; 14; 21; 28; 35; 42; 49; 56; 63; 70; 77; 84; 91; 98; 105; 112; 119; 126; 133; 140; 147; 154; 161; 168; 175; 180];
    I_real = [75; 110; 220; 480; 1100; 2450; 5200; 11400; 23600; 47000; 86500; 142000; 198000; ...
              245000; 272000; 279500; 268000; 242000; 201000; 156000; 114000; 82000; 59000; 44000; 35000; 29500; 26500];

    y0_sir  = [p.N_total - p.I0; p.I0; 0];
    res_sir = lsode(@(y, t) ode_sir(y, t, p), y0_sir, t_span);

    y0_sird  = [p.N_total - p.I0; p.I0; 0; 0];
    res_sird = lsode(@(y, t) ode_sird(y, t, p), y0_sird, t_span);

    y0_sis  = [p.N_total - p.I0; p.I0];
    res_sis = lsode(@(y, t) ode_sis(y, t, p), y0_sis, t_span);

    E0       = round(p.I0 * 1.5);
    y0_seir  = [p.N_total - p.I0 - E0; E0; p.I0; 0];
    res_seir = lsode(@(y, t) ode_seir(y, t, p), y0_seir, t_span);

    N_h0    = round(p.N_total * p.f_high_risk);
    N_g0    = p.N_total - N_h0;
    y0_neo  = [N_g0 - p.I_g0; p.I_g0; 0; 0; 0; N_h0 - p.I_h0; p.I_h0; 0; 0; 0];
    res_neo = lsode(@(y, t) ode_neo_sierd(y, t, p), y0_neo, t_span);

    R_t_sir  = (p.beta / p.gamma) * (res_sir(:, 1) / p.N_total);
    R_t_sird = (p.beta / (p.gamma + p.delta)) * (res_sird(:, 1) / p.N_total);
    R_t_sis  = (p.beta / p.gamma) * (res_sis(:, 1) / p.N_total);
    R_t_seir = (p.beta / p.gamma) * (res_seir(:, 1) / p.N_total);

    k11 = (p.beta_g  * res_neo(:, 1)) / (p.N_total * (p.gamma_g + p.delta_g));
    k12 = (p.beta_gh * res_neo(:, 1)) / (p.N_total * (p.gamma_h + p.delta_h));
    k21 = (p.beta_hg * res_neo(:, 6)) / (p.N_total * (p.gamma_g + p.delta_g));
    k22 = (p.beta_h  * res_neo(:, 6)) / (p.N_total * (p.gamma_h + p.delta_h));
    R_t_neo = 0.5 * ((k11 + k22) + sqrt((k11 - k22).^2 + 4 * k12 .* k21));

    c = palette();
    plot_sir_dashboard(t_span, res_sir, R_t_sir, p, c);
    plot_sird_dashboard(t_span, res_sir, res_sird, R_t_sird, p, c);
    plot_sis_dashboard(t_span, res_sis, R_t_sis, p, c);
    plot_seir_dashboard(t_span, res_sir, res_seir, R_t_seir, p, c);
    plot_neodyn_dashboard(t_span, res_neo, R_t_neo, p, c);
    plot_comparison_dashboard(t_span, res_sir, res_sird, res_sis, res_seir, res_neo, t_real, I_real, c);
end

function plot_sir_dashboard(t, sir, R_t, p, c)
    f = figure(1, 'visible', 'off', 'position', [50, 50, 1200, 1000]);
    hosp_cap = p.N_total * 0.005;

    subplot(3, 2, 1);
    plot(t, sir(:, 1), 'color', c.blue,  'linewidth', 2.2); hold on;
    plot(t, sir(:, 2), 'color', c.red,   'linewidth', 2.2);
    plot(t, sir(:, 3), 'color', c.green, 'linewidth', 2.2);
    setup_axes('Time (Days)', 'Population', '1. Classical SIR: State Trajectories');
    legend({'Susceptible (S)', 'Infectious (I)', 'Recovered (R)'}, 'location', 'northeast', 'fontsize', 8);

    subplot(3, 2, 2);
    plot(sir(:, 1), sir(:, 2), 'color', c.purple, 'linewidth', 2.2); hold on;
    plot(sir(1, 1), sir(1, 2), 'go', 'markersize', 8, 'markerfacecolor', 'g');
    plot(sir(end, 1), sir(end, 2), 'rs', 'markersize', 8, 'markerfacecolor', 'r');
    setup_axes('Susceptible Population (S)', 'Infected Population (I)', '2. SIR Phase-Plane Trajectory (I vs S)');
    legend({'Phase Orbit', 'Start Point', 'End Point'}, 'location', 'northeast', 'fontsize', 8);

    subplot(3, 2, 3);
    plot(t, R_t, 'color', c.purple, 'linewidth', 2.2); hold on;
    plot([t(1), t(end)], [1.0, 1.0], 'r--', 'linewidth', 1.8);
    setup_axes('Time (Days)', 'Effective Reproduction R_t', sprintf('3. SIR: R-effective Curve (R_0 = %.2f)', R_t(1)));
    legend({'R_t Trajectory', 'Threshold (R_t = 1.0)'}, 'location', 'northeast', 'fontsize', 8);

    subplot(3, 2, 4);
    plot(t, sir(:, 2) * 0.05, 'color', c.maroon, 'linewidth', 2.2); hold on;
    plot([t(1), t(end)], [hosp_cap, hosp_cap], 'r--', 'linewidth', 1.8);
    setup_axes('Time (Days)', 'Beds Needed', '4. SIR: Healthcare Bed Capacity Stress');
    legend({'Severe Cases (5% I)', 'Hospital Capacity (0.5% N)'}, 'location', 'northeast', 'fontsize', 8);

    subplot(3, 2, 5);
    plot(t, p.beta * sir(:, 1) .* sir(:, 2) / p.N_total, 'color', c.orange, 'linewidth', 2.2);
    setup_axes('Time (Days)', 'New Cases / Day', '5. SIR: Daily Incidence Rate');

    subplot(3, 2, 6);
    area(t, [sir(:, 2), sir(:, 3), sir(:, 1)]);
    colormap([c.red; c.green; c.blue]);
    setup_axes('Time (Days)', 'Total Population', '6. SIR: Cumulative Stacked Composition');
    legend({'Infected', 'Recovered', 'Susceptible'}, 'location', 'northeast', 'fontsize', 8);

    print(f, 'dashboard_1_sir.png', '-dpng', '-r300');
end

function plot_sird_dashboard(t, sir, sird, R_t, p, c)
    f = figure(2, 'visible', 'off', 'position', [70, 70, 1200, 1000]);
    hosp_cap = p.N_total * 0.005;

    subplot(3, 2, 1);
    plot(t, sird(:, 1), 'color', c.blue,   'linewidth', 2.2); hold on;
    plot(t, sird(:, 2), 'color', c.red,    'linewidth', 2.2);
    plot(t, sird(:, 3), 'color', c.green,  'linewidth', 2.2);
    plot(t, sird(:, 4), 'color', c.black,  'linewidth', 2.2);
    setup_axes('Time (Days)', 'Population', '1. SIRD: State Trajectories with Mortality');
    legend({'S', 'I', 'R', 'Deceased (D)'}, 'location', 'northeast', 'fontsize', 8);

    subplot(3, 2, 2);
    plot(t, sird(:, 4), 'color', c.black, 'linewidth', 2.2);
    setup_axes('Time (Days)', 'Cumulative Deceased', sprintf('2. SIRD: Total Fatalities (Final: %d)', round(sird(end, 4))));

    subplot(3, 2, 3);
    plot(t, R_t, 'color', c.purple, 'linewidth', 2.2); hold on;
    plot([t(1), t(end)], [1.0, 1.0], 'r--', 'linewidth', 1.8);
    setup_axes('Time (Days)', 'Effective Reproduction R_t', sprintf('3. SIRD: R-effective Curve (R_0 = %.2f)', R_t(1)));
    legend({'R_t Trajectory', 'Threshold (R_t = 1.0)'}, 'location', 'northeast', 'fontsize', 8);

    subplot(3, 2, 4);
    plot(t, sird(:, 2) * 0.05, 'color', c.maroon, 'linewidth', 2.2); hold on;
    plot([t(1), t(end)], [hosp_cap, hosp_cap], 'r--', 'linewidth', 1.8);
    setup_axes('Time (Days)', 'Beds Needed', '4. SIRD: Hospital Capacity Demand');
    legend({'Hospital Demand (5% I)', 'Capacity (0.5% N)'}, 'location', 'northeast', 'fontsize', 8);

    subplot(3, 2, 5);
    plot(t, p.delta * sird(:, 2), 'color', c.red, 'linewidth', 2.2);
    setup_axes('Time (Days)', 'Deaths / Day', '5. SIRD: Daily Mortality Flux (dD/dt)');

    subplot(3, 2, 6);
    plot(t, sir(:, 2),  'color', c.blue, 'linewidth', 2.2); hold on;
    plot(t, sird(:, 2), 'color', c.red,  'linewidth', 2.2);
    setup_axes('Time (Days)', 'Infected Population', '6. Mortality Impact on Active Infection Peak');
    legend({'SIR Peak (No Mortality)', 'SIRD Peak (With Mortality)'}, 'location', 'northeast', 'fontsize', 8);

    print(f, 'dashboard_2_sird.png', '-dpng', '-r300');
end

function plot_sis_dashboard(t, sis, R_t, p, c)
    f = figure(3, 'visible', 'off', 'position', [90, 90, 1200, 1000]);
    N = p.N_total;
    hosp_cap = N * 0.005;
    I_star = N * (1 - p.gamma / p.beta);

    subplot(3, 2, 1);
    plot(t, sis(:, 1), 'color', c.blue, 'linewidth', 2.2); hold on;
    plot(t, sis(:, 2), 'color', c.red,  'linewidth', 2.2);
    setup_axes('Time (Days)', 'Population', '1. SIS: Endemic Persistence vs Susceptibles');
    legend({'Susceptible (S)', 'Infectious (I)'}, 'location', 'northeast', 'fontsize', 8);

    subplot(3, 2, 2);
    plot(t, sis(:, 2), 'color', c.red, 'linewidth', 2.2); hold on;
    plot([t(1), t(end)], [I_star, I_star], 'k--', 'linewidth', 1.8);
    setup_axes('Time (Days)', 'Infected Individuals', sprintf('2. SIS: Convergence to Endemic Attractor I^* = %d', round(I_star)));
    legend({'Active Cases I(t)', sprintf('Analytical I^* = %d', round(I_star))}, 'location', 'southeast', 'fontsize', 8);

    subplot(3, 2, 3);
    plot(t, R_t, 'color', c.purple, 'linewidth', 2.2); hold on;
    plot([t(1), t(end)], [1.0, 1.0], 'r--', 'linewidth', 1.8);
    setup_axes('Time (Days)', 'Effective Reproduction R_t', sprintf('3. SIS: R-effective Trajectory (R_0 = %.2f)', R_t(1)));
    legend({'R_t Trajectory', 'Endemic Equilibrium (R_t = 1.0)'}, 'location', 'northeast', 'fontsize', 8);

    subplot(3, 2, 4);
    plot(t, sis(:, 2) * 0.05, 'color', c.maroon, 'linewidth', 2.2); hold on;
    plot([t(1), t(end)], [hosp_cap, hosp_cap], 'r--', 'linewidth', 1.8);
    setup_axes('Time (Days)', 'Beds Needed', '4. SIS: Chronic Healthcare Resource Demand');
    legend({'Continuous Bed Need (5% I)', 'Capacity Limit'}, 'location', 'northeast', 'fontsize', 8);

    subplot(3, 2, 5);
    plot(t, p.gamma * sis(:, 2), 'color', c.orange, 'linewidth', 2.2);
    setup_axes('Time (Days)', 'Individuals / Day', '5. SIS: Daily Re-infection Flux (I \to S)');

    subplot(3, 2, 6);
    I_grid = linspace(0, N, 200);
    plot(I_grid, (p.beta - p.gamma) * I_grid .* (1 - I_grid / I_star), 'color', c.red, 'linewidth', 2.2); hold on;
    plot([0, N], [0, 0], 'k--');
    plot(I_star, 0, 'ro', 'markersize', 8, 'markerfacecolor', 'r');
    setup_axes('Infected Population (I)', 'dI/dt Rate', '6. SIS: Logistic Growth Parabola (dI/dt vs I)');
    legend({'dI/dt Curve', 'Zero Net Growth', 'Endemic Attractor I^*'}, 'location', 'northeast', 'fontsize', 8);

    print(f, 'dashboard_3_sis.png', '-dpng', '-r300');
end

function plot_seir_dashboard(t, sir, seir, R_t, p, c)
    f = figure(4, 'visible', 'off', 'position', [110, 110, 1200, 1000]);
    hosp_cap = p.N_total * 0.005;

    subplot(3, 2, 1);
    plot(t, seir(:, 1), 'color', c.blue,   'linewidth', 2.2); hold on;
    plot(t, seir(:, 2), 'color', c.orange, 'linewidth', 2.2);
    plot(t, seir(:, 3), 'color', c.red,    'linewidth', 2.2);
    plot(t, seir(:, 4), 'color', c.green,  'linewidth', 2.2);
    setup_axes('Time (Days)', 'Population', '1. SEIR: Overall Population Dynamics');
    legend({'S', 'Latent (E)', 'Infectious (I)', 'R'}, 'location', 'northeast', 'fontsize', 8);

    subplot(3, 2, 2);
    plot(t, seir(:, 2), 'color', c.orange, 'linewidth', 2.2); hold on;
    plot(t, seir(:, 3), 'color', c.red,    'linewidth', 2.2);
    setup_axes('Time (Days)', 'Individuals', '2. SEIR: Latent (E) vs Infectious (I)');
    legend({'Latent (E)', 'Infectious (I)'}, 'location', 'northeast', 'fontsize', 8);

    subplot(3, 2, 3);
    plot(t, R_t, 'color', c.purple, 'linewidth', 2.2); hold on;
    plot([t(1), t(end)], [1.0, 1.0], 'r--', 'linewidth', 1.8);
    setup_axes('Time (Days)', 'R-effective', sprintf('3. SEIR: R-effective (R_0 = %.2f)', R_t(1)));
    legend({'R_t Trajectory', 'Threshold (R_t = 1.0)'}, 'location', 'northeast', 'fontsize', 8);

    subplot(3, 2, 4);
    plot(t, seir(:, 3) * 0.05, 'color', c.maroon, 'linewidth', 2.2); hold on;
    plot([t(1), t(end)], [hosp_cap, hosp_cap], 'r--', 'linewidth', 1.8);
    setup_axes('Time (Days)', 'Beds Needed', '4. SEIR: Hospital Bed Demand');
    legend({'Beds Needed (5% I)', 'Capacity (0.5% N)'}, 'location', 'northeast', 'fontsize', 8);

    subplot(3, 2, 5);
    plot(t, p.sigma * seir(:, 2), 'color', c.orange, 'linewidth', 2.2);
    setup_axes('Time (Days)', 'Persons / Day', '5. SEIR: Incubation Breakthrough Flux');

    subplot(3, 2, 6);
    plot(t, sir(:, 2),  'color', c.blue, 'linewidth', 2.2); hold on;
    plot(t, seir(:, 3), 'color', c.red,  'linewidth', 2.2);
    setup_axes('Time (Days)', 'Active Cases', '6. Latency Damping: SIR vs SEIR Peak Shift');
    legend({'Classical SIR (Day 55)', 'SEIR (Day 108)'}, 'location', 'northeast', 'fontsize', 8);

    print(f, 'dashboard_4_seir.png', '-dpng', '-r300');
end

function plot_neodyn_dashboard(t, neo, R_t, p, c)
    f = figure(5, 'visible', 'off', 'position', [130, 130, 1200, 1000]);
    N = p.N_total;
    hosp_cap = N * 0.005;

    S_g = neo(:, 1); I_g = neo(:, 2); R_g = neo(:, 3); D_g = neo(:, 4); V_g = neo(:, 5);
    S_h = neo(:, 6); I_h = neo(:, 7); D_h = neo(:, 9);
    I_neo = I_g + I_h;
    D_neo = D_g + D_h;
    V_neo = neo(:, 5) + neo(:, 10);

    subplot(3, 2, 1);
    plot(t, S_g, 'color', c.blue,   'linewidth', 2.2); hold on;
    plot(t, I_g, 'color', c.red,    'linewidth', 2.2);
    plot(t, R_g, 'color', c.green,  'linewidth', 2.2);
    plot(t, V_g, 'color', c.purple, 'linewidth', 2.2, 'linestyle', '--');
    setup_axes('Time (Days)', 'Population Count', '1. Overall Population Dynamics');
    legend({'Susceptible', 'Infected', 'Recovered', 'Vaccinated'}, 'location', 'northeast', 'fontsize', 8);

    subplot(3, 2, 2);
    plot(t, I_h,   'color', [0.75, 0.1, 0.1], 'linewidth', 2.2); hold on;
    plot(t, I_neo, 'color', c.black,          'linewidth', 1.8, 'linestyle', ':');
    setup_axes('Time (Days)', 'Active Cases', '2. High-Risk Groups Infection');
    legend({'High-Risk (I_h)', 'Total Active (I)'}, 'location', 'northeast', 'fontsize', 8);

    subplot(3, 2, 3);
    plot(t, R_t, 'color', c.purple, 'linewidth', 2.2); hold on;
    plot([t(1), t(end)], [1.0, 1.0], 'r--', 'linewidth', 1.8);
    ylim([0, max(max(R_t) * 1.15, 2.0)]);
    setup_axes('Time (Days)', 'R-effective', '3. R-effective & Intervention Impact');
    legend({'R-effective', 'Threshold (R = 1.0)'}, 'location', 'northeast', 'fontsize', 8);

    subplot(3, 2, 4);
    plot(t, I_neo * 0.05, 'color', [0.55, 0.0, 0.0], 'linewidth', 2.2); hold on;
    plot([t(1), t(end)], [hosp_cap, hosp_cap], 'r--', 'linewidth', 1.8);
    setup_axes('Time (Days)', 'Beds Needed', '4. Hospital Burden & Capacity');
    legend({'Hospital Beds Needed', sprintf('Capacity (0.5%% N = %d)', round(hosp_cap))}, 'location', 'northeast', 'fontsize', 8);

    subplot(3, 2, 5);
    vacc_pct = (V_neo / N) * 100;
    herd_imm = max(0, (1 - 1 / R_t(1)) * 100);
    plot(t, vacc_pct, 'color', c.green, 'linewidth', 2.2); hold on;
    plot([t(1), t(end)], [herd_imm, herd_imm], 'color', [0.0, 0.5, 0.0], 'linestyle', '--', 'linewidth', 1.8);
    ylim([0, 100]);
    setup_axes('Time (Days)', 'Vaccinated (%)', '5. Vaccination Progress');
    legend({'Vaccination %', sprintf('Herd Immunity (%.1f%%)', herd_imm)}, 'location', 'northwest', 'fontsize', 8);

    subplot(3, 2, 6);
    b = bar([D_g(end), D_h(end), D_neo(end)]);
    set(b, 'facecolor', [0.55, 0.1, 0.1]);
    set(gca, 'xtick', 1:3, 'xticklabel', {'General Pop', 'High-Risk', 'Total Deaths'});
    setup_axes('', 'Final Deaths', '6. Mortality by Risk Group');

    print(f, 'dashboard_5_neodyn.png', '-dpng', '-r300');
    copyfile('dashboard_5_neodyn.png', 'neo_dyn_sierd_dashboard.png');
end

function plot_comparison_dashboard(t, sir, sird, sis, seir, neo, t_real, I_real, c)
    f = figure(6, 'visible', 'off', 'position', [150, 150, 1300, 1050]);

    S_neo = neo(:, 1) + neo(:, 6);
    I_neo = neo(:, 2) + neo(:, 7);
    D_neo = neo(:, 4) + neo(:, 9);

    subplot(3, 2, 1);
    plot(t_real, I_real, 'ko', 'markersize', 6, 'markerfacecolor', 'k'); hold on;
    plot(t, sir(:, 2),   'color', c.blue,   'linewidth', 1.8);
    plot(t, sird(:, 2),  'color', c.orange, 'linewidth', 1.8);
    plot(t, seir(:, 3),  'color', c.purple, 'linewidth', 1.8, 'linestyle', '-.');
    plot(t, I_neo,       'color', c.red,    'linewidth', 2.4);
    setup_axes('Time (Days)', 'Active Cases', '1. Real Mpox Data Fit: Models vs Observations');
    legend({'Real Observed Cases', 'Classical SIR', 'SIRD', 'SEIR', 'Neo-Dynamic SIeRd'}, 'location', 'northeast', 'fontsize', 8);

    subplot(3, 2, 2);
    plot(t, sir(:, 2),  'color', c.blue,   'linewidth', 2.0); hold on;
    plot(t, sird(:, 2), 'color', c.orange, 'linewidth', 2.0);
    plot(t, sis(:, 2),  'color', c.green,  'linewidth', 1.8, 'linestyle', '--');
    plot(t, seir(:, 3), 'color', c.purple, 'linewidth', 2.0, 'linestyle', '-.');
    plot(t, I_neo,      'color', c.red,    'linewidth', 2.4);
    setup_axes('Time (Days)', 'Active Cases', '2. Active Infection Curves: All 5 Models');
    legend({'SIR', 'SIRD', 'SIS (Endemic)', 'SEIR', 'Neo-Dynamic'}, 'location', 'northeast', 'fontsize', 8);

    subplot(3, 2, 3);
    plot(t, sir(:, 1),  'color', c.blue,   'linewidth', 2.0); hold on;
    plot(t, sird(:, 1), 'color', c.orange, 'linewidth', 2.0);
    plot(t, sis(:, 1),  'color', c.green,  'linewidth', 1.8, 'linestyle', '--');
    plot(t, seir(:, 1), 'color', c.purple, 'linewidth', 2.0, 'linestyle', '-.');
    plot(t, S_neo,      'color', c.red,    'linewidth', 2.4);
    setup_axes('Time (Days)', 'Susceptible Population (S)', '3. Susceptible Depletion Comparison');
    legend({'SIR', 'SIRD', 'SIS', 'SEIR', 'Neo-Dynamic'}, 'location', 'northeast', 'fontsize', 8);

    subplot(3, 2, 4);
    plot(t, sird(:, 4), 'color', c.orange, 'linewidth', 2.0); hold on;
    plot(t, D_neo,      'color', c.red,    'linewidth', 2.2);
    plot(t, sis(:, 2),  'color', c.green,  'linewidth', 1.8, 'linestyle', '--');
    setup_axes('Time (Days)', 'Persons', '4. Fatalities (SIRD, Neo-Dyn) & Endemic (SIS)');
    legend({'SIRD Deaths', 'Neo-Dynamic Deaths', 'SIS Active Pool'}, 'location', 'northwest', 'fontsize', 8);

    subplot(3, 2, 5);
    mod_names = {'Neo-Dyn', 'SEIR', 'SIRD', 'SIR', 'SIS'};
    mod_preds = {I_neo, seir(:, 3), sird(:, 2), sir(:, 2), sis(:, 2)};
    rmse_vals = zeros(1, 5);
    mae_vals  = zeros(1, 5);
    for m = 1:5
        pred_obs = interp1(t, mod_preds{m}, t_real, 'linear');
        res = pred_obs - I_real;
        rmse_vals(m) = sqrt(mean(res.^2));
        mae_vals(m)  = mean(abs(res));
    end
    b = bar([rmse_vals; mae_vals]');
    set(b(1), 'facecolor', [0.15, 0.45, 0.75]);
    set(b(2), 'facecolor', [0.85, 0.40, 0.15]);
    set(gca, 'xtick', 1:5, 'xticklabel', mod_names);
    setup_axes('', 'Cases Error', '5. Empirical Error Comparison (RMSE vs MAE)');
    legend({'RMSE Error', 'MAE Error'}, 'location', 'northwest', 'fontsize', 8);

    subplot(3, 2, 6);
    peak_vals = [max(sir(:, 2)), max(sird(:, 2)), max(sis(:, 2)), max(seir(:, 3)), max(I_neo)];
    b_peak = bar(peak_vals);
    set(b_peak, 'facecolor', [0.70, 0.20, 0.20]);
    set(gca, 'xtick', 1:5, 'xticklabel', {'SIR', 'SIRD', 'SIS', 'SEIR', 'Neo-Dyn'});
    setup_axes('', 'Peak Active Cases (I_{max})', '6. Peak Infection Burden Comparison');

    print(f, 'dashboard_6_cross_comparison_real_data.png', '-dpng', '-r300');
    copyfile('dashboard_6_cross_comparison_real_data.png', 'neo_dyn_sierd_comparisons.png');
end

function setup_axes(x_lbl, y_lbl, ttl)
    if ~isempty(x_lbl)
        xlabel(x_lbl);
    end
    if ~isempty(y_lbl)
        ylabel(y_lbl);
    end
    if ~isempty(ttl)
        title(ttl);
    end
    grid on;
    box on;
end

function c = palette()
    c.blue   = [0.12, 0.47, 0.71];
    c.red    = [0.84, 0.19, 0.15];
    c.green  = [0.17, 0.63, 0.17];
    c.purple = [0.58, 0.40, 0.74];
    c.orange = [1.00, 0.50, 0.05];
    c.black  = [0.20, 0.20, 0.20];
    c.maroon = [0.55, 0.00, 0.00];
end

function dydt = ode_sir(y, t, p)
    S = max(0, y(1)); I = max(0, y(2)); R = max(0, y(3));
    N = p.N_total;
    dS = -p.beta * S * I / N;
    dI =  p.beta * S * I / N - p.gamma * I;
    dR =  p.gamma * I;
    dydt = [dS; dI; dR];
end

function dydt = ode_sird(y, t, p)
    S = max(0, y(1)); I = max(0, y(2)); R = max(0, y(3)); D = max(0, y(4));
    N = p.N_total;
    dS = -p.beta * S * I / N;
    dI =  p.beta * S * I / N - (p.gamma + p.delta) * I;
    dR =  p.gamma * I;
    dD =  p.delta * I;
    dydt = [dS; dI; dR; dD];
end

function dydt = ode_sis(y, t, p)
    S = max(0, y(1)); I = max(0, y(2));
    N = p.N_total;
    dS = -p.beta * S * I / N + p.gamma * I;
    dI =  p.beta * S * I / N - p.gamma * I;
    dydt = [dS; dI];
end

function dydt = ode_seir(y, t, p)
    S = max(0, y(1)); E = max(0, y(2)); I = max(0, y(3)); R = max(0, y(4));
    N = p.N_total;
    dS = -p.beta * S * I / N;
    dE =  p.beta * S * I / N - p.sigma * E;
    dI =  p.sigma * E - p.gamma * I;
    dR =  p.gamma * I;
    dydt = [dS; dE; dI; dR];
end

function dydt = ode_neo_sierd(y, t, p)
    S_g = max(0, y(1)); I_g = max(0, y(2)); R_g = max(0, y(3)); D_g = max(0, y(4)); V_g = max(0, y(5));
    S_h = max(0, y(6)); I_h = max(0, y(7)); R_h = max(0, y(8)); D_h = max(0, y(9)); V_h = max(0, y(10));
    N = p.N_total;

    lambda_g = p.beta_g  * I_g / N + p.beta_gh * I_h / N;
    lambda_h = p.beta_h  * I_h / N + p.beta_hg * I_g / N;

    dS_g = -lambda_g * S_g - p.v_g * S_g;
    dI_g =  lambda_g * S_g - (p.gamma_g + p.delta_g) * I_g;
    dR_g =  p.gamma_g * I_g + p.efficacy * p.v_g * S_g;
    dD_g =  p.delta_g * I_g;
    dV_g =  p.v_g * S_g;

    dS_h = -lambda_h * S_h - p.v_h * S_h;
    dI_h =  lambda_h * S_h - (p.gamma_h + p.delta_h) * I_h;
    dR_h =  p.gamma_h * I_h + p.efficacy * p.v_h * S_h;
    dD_h =  p.delta_h * I_h;
    dV_h =  p.v_h * S_h;

    dydt = [dS_g; dI_g; dR_g; dD_g; dV_g; dS_h; dI_h; dR_h; dD_h; dV_h];
end
