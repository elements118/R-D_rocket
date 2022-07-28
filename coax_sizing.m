clc;
clear all; 

% Global Constants
gravity = 32.3; % [ft/s^2] Acceleration due to gravity

% Inner Element constants and assumptions
m_dot_LOX = 1; % [lbm/s] Oxidizer mass flow  
delta_p_LOX = 14400; % [lbf/ft^2] Oxidizer pressure drop 
LOX_dens = 71.168; % [lb/ft^3] LOX density 
inner_num_inlets = 4; % [N/A] number of tangential inlets 
spray_angle = 30; % [degrees] Desired spray half angle 
kin_visc_LOX= 2.362 * 10^-6; % [ft^2/s] kinematic viscosity of LOX
K_guess = 1.8;

% Initial Calculations
inner_filling_eff = fzero(@(inner_filling_eff) K_guess - (((1-inner_filling_eff) * sqrt(2)) ...
    / (inner_filling_eff * sqrt(inner_filling_eff))), .4); % [N/A] Beyvel and Orzechowski, Eq. 5-65

inner_disc_coeff = inner_filling_eff * sqrt(inner_filling_eff / (2 - inner_filling_eff)); % [N/A] Beyvel and Orzechowski, Eq. 5-66

inner_swirl_diam = sqrt(4 * m_dot_LOX / (pi * inner_disc_coeff * sqrt(2 * LOX_dens * delta_p_LOX * gravity))); % Beyvel, Eq. 5-82

inner_nozzle_length = 2 * inner_swirl_diam; % Suggested by Bazarov

R = 2 * inner_swirl_diam; % [N/A] Find this ratio, from Beyvel pg. 263 should be (2-5) times inner swirl diam

inner_chamber_length = 3 * R; % From Bazarov, should be (2-3) * R

inner_inlet_diam = sqrt((2 * R * inner_swirl_diam) / (inner_num_inlets * K_guess)); % Beyvel, Eq. 5-84

inner_inlet_length = (inner_inlet_diam / 2) * 3; % Suggested by Bazarov to be (3-4) * inlet radius

fcn = @(inner_S) ((sqrt(1 - inner_disc_coeff^2 * K_guess^2)) ...
    - (inner_S * sqrt(inner_S^2 - inner_disc_coeff^2 * K_guess^2)) ...
    - (inner_disc_coeff^2 * K_guess^2 * log((1 ...
    + sqrt(1 - inner_disc_coeff^2 * K_guess^2)) ...
    / (inner_S + sqrt(inner_S^2 - inner_disc_coeff^2 ...
    * K_guess^2)))) - inner_disc_coeff); % Beyvel Eq. 5-58

inner_S = fzero(fcn, 1);

Re_LOX = (4 * m_dot_LOX)/ (pi * LOX_dens * ...
    kin_visc_LOX * sqrt(inner_num_inlets) * inner_inlet_diam); % Beyvel, Eq. 5-86

inner_frict_coeff = exp((25.8/log(Re_LOX)^2.58) - 2); % Beyvel, Eq. 5-87

K_visc = (R * (inner_swirl_diam / 2)) / ((inner_num_inlets * ...
    (inner_inlet_diam / 2)^2) + (inner_frict_coeff / 2) * R * (R - (inner_swirl_diam / 2))); % Beyvel, Eq. 5-81

visc_spray_angle = atand((2 * inner_disc_coeff * K_visc) / sqrt((1 + inner_S)^2 - ...
    (4 * inner_disc_coeff^2 * K_guess^2))); % Beyvel, Eq. 5-80

% Inner Element constants and assumptions
inner_wall_thck = .005104; % [ft] wall thickness of the inner element
m_dot_FUEL = 2.5; % [lbm/s] mass flow of fuel
delta_p_FUEL = 14400; % [lbf/ft^2] pressure drop of fuel 
FUEL_dens = 50.566; % [lbm/ft^3] fuel density
outer_num_inlets = 4; % [N/A] number of inlets of the outer element
kin_visc_FUEL = 2.368*10^-5; % [ft^2/s] kinematic viscosity of the fuel
vortex_gap = .01181 / 12; % [ft] permitted gas vortex radius, Bazarov page 76

lcv = 1;
tf = []; % creating empty matrix
for tf_guess = .001:.01:.20
    outer_swirl_diam = (inner_swirl_diam + ...
        (2 * inner_wall_thck) + (2 * vortex_gap) + (2 * tf_guess)); % [] Outer swirl outlet diameter

    effective_flow_area = m_dot_FUEL / sqrt(2 * FUEL_dens * delta_p_FUEL * gravity); % [] Beyvel, Eq. 5-63

    outlet_area = (pi / 4) * outer_swirl_diam^2; % [] 

    outer_disc_coeff = effective_flow_area / outlet_area; % [N/A]

    fcn = @(outer_filling_eff) (outer_filling_eff * sqrt(outer_filling_eff / (2 - outer_filling_eff)) - outer_disc_coeff);

    outer_filling_eff = fzero(fcn, .7);

    outer_K = (((1 - outer_filling_eff) * sqrt(2)) ...
        / (outer_filling_eff * sqrt(outer_filling_eff))); % [N/A] Beyvel, Eq. 5-65

    fcn = @(outer_S) (sqrt(1 - outer_disc_coeff^2 * outer_K^2)) ...
        - (outer_S * sqrt(outer_S^2 - outer_disc_coeff^2 * outer_K^2)) ...
        - (outer_disc_coeff^2 * outer_K^2 * log((1 ...
        + sqrt(1 - outer_disc_coeff^2 * outer_K^2)) ...
        / (outer_S + sqrt(outer_S^2 - outer_disc_coeff^2 ...
        * outer_K^2)))) - outer_disc_coeff; % Beyvel Eq. 5-58

     outer_S = fzero(fcn, .8);

     outer_spray_angle = atan((2 * outer_disc_coeff * outer_K) / (sqrt((1 + outer_S)^2 - ...
        (4 * outer_disc_coeff^2 * outer_K^2)))); % [] Beyvel Eq. 5-80

     outer_tf = .5 * (outer_swirl_diam - (outer_S * outer_swirl_diam)); % [] outer film thickness

     tf(lcv) = outer_tf;

     lcv = lcv + 1;
end 




























