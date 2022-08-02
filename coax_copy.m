%{
Purdue Space Program - Liquids
Research and Development - coaxial bi-swirl sizing code
Created by: Jacob Bell
In this early version, a  mix of both method's presented by Beyvel and Bazarov are used. 
Essentially just means that forms of
certain equations are expressed in varying ways (see document on the
drive for direct equation comparison).
%}
clc;
clear; 

% Global Constants
gravity = 32.3; % [ft/s^2] Acceleration due to gravity

% Inner Element constants and assumptions
m_dot_LOX = 1; % [lbm/s] Oxidizer mass flow  
delta_p_LOX = 10800; % [lbf/ft^2] Oxidizer pressure drop 
LOX_dens = 71.168; % [lb/ft^3] LOX density 
inner_num_inlets = 3; % [N/A] number of tangential inlets 
spray_angle = 30; % [degrees] Desired spray half angle 
kin_visc_LOX= 2.362 * 10^-6; % [ft^2/s] kinematic viscosity of LOX
K_guess = 3;
inner_wall_thck = .005104; % [ft] wall thickness of the inner element, ~ 1/16"
coeff_nozzle_open = 3; % coefficient of nozzle opening, from Beyvel pg. 263 
    % should be (2-5), Bazarov says 3x. This value is the ratio of "swirl
    % arm"/nozzle radius

% Initial Calculations

lcv = 1;
while lcv < 100
    inner_filling_eff = fzero(@(inner_filling_eff) K_guess - (((1-inner_filling_eff) * sqrt(2)) ...
        / (inner_filling_eff * sqrt(inner_filling_eff))), .4); % [N/A] Beyvel and Orzechowski, Eq. 5-65
    
    inner_disc_coeff = inner_filling_eff * sqrt(inner_filling_eff / (2 - inner_filling_eff)); % [N/A] Beyvel and Orzechowski, Eq. 5-66
    
    fcn = @(inner_S) ((sqrt(1 - inner_disc_coeff^2 * K_guess^2)) ...
    - (inner_S * sqrt(inner_S^2 - inner_disc_coeff^2 * K_guess^2)) ...
    - (inner_disc_coeff^2 * K_guess^2 * log((1 ...
    + sqrt(1 - inner_disc_coeff^2 * K_guess^2)) ...
    / (inner_S + sqrt(inner_S^2 - inner_disc_coeff^2 ...
    * K_guess^2)))) - inner_disc_coeff); % Beyvel Eq. 5-58

    inner_S = fzero(fcn, .8);

    inner_swirl_diam = sqrt(4 * m_dot_LOX / (pi * inner_disc_coeff * ...
        sqrt(2 * LOX_dens * delta_p_LOX * gravity))); % Beyvel, Eq. 5-82
    
    inner_R = (coeff_nozzle_open * inner_swirl_diam) / 2; % "swirl arm", center axis to inlet
    
    inner_inlet_diam = sqrt((2 * inner_R * inner_swirl_diam) / (inner_num_inlets * K_guess)); % Beyvel, Eq. 5-84

    Re_LOX = (4 * m_dot_LOX)/ (pi * LOX_dens * ...
        kin_visc_LOX * sqrt(inner_num_inlets) * inner_inlet_diam); % Beyvel, Eq. 5-86
    
    inner_frict_coeff = exp((25.8/log(Re_LOX)^2.58) - 2); % Friction Coefficient Beyvel, Eq. 5-87
    
    K_frict = (inner_R * (inner_swirl_diam / 2)) / ((inner_num_inlets * ...
        (inner_inlet_diam / 2)^2) + (inner_frict_coeff / 2) * ...
        inner_R * (inner_R - (inner_swirl_diam / 2))); % K value taking friction into account Beyvel, Eq. 5-81

    visc_spray_angle = atand((2 * inner_disc_coeff * K_frict) / sqrt((1 + inner_S)^2 - ...
    (4 * inner_disc_coeff^2 * K_guess^2))); % Beyvel, Eq. 5-80

    eq_filling_eff = fzero(@(eq_filling_eff) K_frict - (((1-eq_filling_eff) * sqrt(2)) ...
        / (eq_filling_eff * sqrt(eq_filling_eff))), .4); % [N/A] Beyvel and Orzechowski, Eq. 5-65

    eq_disc_coeff = eq_filling_eff * sqrt(eq_filling_eff / ...
        (2 - eq_filling_eff)); % temporary discharge coeff from K_frict value

    inner_visc_disc_coeff = eq_disc_coeff / (sqrt(1 + ...
        eq_disc_coeff^2 * (K_guess^2 / coeff_nozzle_open^2))); % Bazarov Eq. 99

    inner_visc_swirl_diam = sqrt(4 * m_dot_LOX / (pi * inner_visc_disc_coeff * ...
        sqrt(2 * LOX_dens * delta_p_LOX * gravity))); % Beyvel, Eq. 5-82

    K_visc = (2 * inner_R * inner_visc_swirl_diam) / (inner_num_inlets * ...
        inner_inlet_diam^2); % Beyvel 5-83

    lcv = lcv + 1;
    if ((K_guess - K_visc) / K_visc < .04)
        K_final = K_visc;
        inner_diam_final = inner_swirl_diam;
        external_nozzle_diam = inner_diam_final + (2 * inner_wall_thck); % rough guess, wall = 1/16"
        break
    else
        K_guess = K_visc;
    end
end

% Inner swirler final parameter calcs
inner_nozzle_length = 2 * inner_swirl_diam; % Suggested by Bazarov pg. 76
inner_chamber_length = 3 * inner_R; % From Bazarov, should be (2-3) * R
inner_inlet_length = (inner_inlet_diam / 2) * 3; % Suggested by Bazarov to be (3-4) * inlet radius
inner_chamber_diam = (2 * inner_R) + inner_inlet_diam;

% TIME FOR PART TWO!!!!!!

% Outer Element constants and assumptions
m_dot_FUEL = .313; % [lbm/s] Oxidizer mass flow  
delta_p_FUEL = 10800; % [lbf/ft^2] Oxidizer pressure drop 
FUEL_dens = 50.57; % [lb/ft^3] LOX density 
outer_num_inlets = 4; % [N/A] number of tangential inlets 
kin_visc_FUEL= 2.362 * 10^-6; % [ft^2/s] kinematic viscosity of LOX
coeff_nozzle_open = 3; % coefficient of nozzle opening, from Beyvel pg. 263 
    % should be (2-5), Bazarov says 3x (for closed). This value is the ratio of "swirl
    % arm"/nozzle radius
permitted_vortex_rad = (external_nozzle_diam / 2) + .00098425; % value comes from bazarov, Pg. 76

% INITIAL VALUES
outer_swirl_diam = 2 * permitted_vortex_rad; % first estimation of outer diameter
effective_flow_area = m_dot_FUEL / sqrt(2 * FUEL_dens * delta_p_FUEL * gravity); % [] Beyvel, Eq. 5-63
outlet_area = (pi / 4) * outer_swirl_diam^2; % [] 
outer_disc_coeff = effective_flow_area / outlet_area; % [N/A]

lcv = 1;
while lcv < 100
    fcn = @(outer_filling_eff) (outer_filling_eff * sqrt(outer_filling_eff ...
        / (2 - outer_filling_eff)) - outer_disc_coeff); % Beyvel Eq. 5-66

    outer_filling_eff = fzero(fcn, .7);

    outer_K = (((1 - outer_filling_eff) * sqrt(2)) ...
        / (outer_filling_eff * sqrt(outer_filling_eff))); % [N/A] Beyvel, Eq. 5-65

    rel_vortex_rad = .1827 * log(7 * outer_K); % from empirical plot

    outer_swirl_diam_new = 2 * (permitted_vortex_rad / rel_vortex_rad); 

    effective_flow_area = m_dot_FUEL / sqrt(2 * FUEL_dens * delta_p_FUEL * gravity); % [] Beyvel, Eq. 5-63

    outlet_area = (pi / 4) * outer_swirl_diam_new^2; % [] 

    outer_disc_coeff = effective_flow_area / outlet_area; % [N/A]

     if abs((outer_swirl_diam_new - outer_swirl_diam) / outer_swirl_diam) < .005
         break
     else
         outer_swirl_diam = outer_swirl_diam_new;
     end
     lcv = lcv + 1;
end

% Outer element final parameters
outer_R = (outer_swirl_diam * coeff_nozzle_open) / 2;
outer_inlet_diam = sqrt((2 * outer_R * outer_swirl_diam) / (outer_num_inlets * outer_K)); % Beyvel, Eq. 5-84
outer_nozzle_length = 2 * outer_swirl_diam; % Suggested by Bazarov pg. 76
outer_chamber_length = outer_R; % can really be anything reasonable
outer_inlet_length = (outer_inlet_diam / 2) * 3; % Suggested by Bazarov to be (3-4) * inlet radius
outer_chamber_diam = (2 * outer_R) + outer_inlet_diam;

% Inner swirler output table
inner_params = ["Inner Swirl Diam."; "Inner Inlet Diam."; "Inner Chamber Lg."; "Inner Inlet Lg."; "Inner Chamber diam."];
inner_value = [inner_diam_final * 12; inner_inlet_diam * 12; inner_chamber_length * 12; inner_inlet_length * 12; inner_chamber_diam * 12];
units = ["[in]"; "[in]"; "[in]"; "[in]"; "[in]"];
inner_output = table(inner_params, inner_value, units);

% Outer swirler output table
outer_params = ["Outer Swirl Diam."; "Outer Inlet Diam."; "Outer Chamber Lg."; "Outer Inlet Lg."; "Outer Chamber diam."];
outer_value = [outer_swirl_diam_new * 12; outer_inlet_diam * 12; outer_chamber_length * 12; outer_inlet_length * 12; outer_chamber_diam * 12];
outer_output = table(outer_params, outer_value, units);

















