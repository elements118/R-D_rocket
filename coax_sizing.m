clc;
clear all; 

% Inner Element constants and assumptions
m_dot_LOX = ; % [] Oxidizer mass flow  
delta_p_LOX = ; % [ Oxidizer pressure drop 
LOX_dens = ; % [] LOX density 
inner_num_inlets = ; % [N/A] number of tangential inlets 
R = ; % [N/A] Find this ratio
spray_angle = ; % [degrees] Desired spray half angle 
kin_visc_LOX= ; % [] kinematic viscosity of LOX

% Initial Calculations
lcv = 0;
angle_diff = [];
for K_guess = 0:.05:10
    inner_filling_eff = fzero(@(inner_filling_eff) K_guess - (((1-inner_filling_eff) * sqrt(2)) ...
        / (inner_filling_eff * sqrt(inner_filling_eff))), .3); % [N/A] Beyvel and Orzechowski, Eq. 5-65

    inner_disc_coeff = inner_filling_eff * sqrt(inner_filling_eff / (2 - inner_filling_eff)); % [N/A] Beyvel and Orzechowski, Eq. 5-66

    inner_swirl_diam = sqrt(4 * m_dot_LOX / (pi * inner_disc_coeff * sqrt(2 * LOX_dens * delta_p_LOX))); % Beyvel, Eq. 5-82

    inner_inlet_diam = sqrt((2 * R * inner_swirl_diam) / (inner_num_inlets * K_guess)); % Beyvel, Eq. 5-84

    inner_S = fzero(@(inner_S) ((sqrt(1-(inner_disc_coeff^2 * K_guess^2))) - ...
        (inner_S * sqrt(inner_S^2-(inner_disc_coeff^2 * K_guess^2))) - ...
        (inner_disc_coeff^2 * K_guess^2 * log((1 + sqrt(1 - inner_disc_coeff^2 * K_guess^2))/(S + ...
        sqrt(inner_S^2 - inner_disc_coeff^2 * K_guess^2)))) - inner_disc_coeff), .8); % Justify this

    Re_LOX = 4 * m_dot_LOX / (pi * LOX_dens * ...
        kin_visc_LOX * sqrt(inner_num_inlets) * inner_inlet_diam); % Beyvel, Eq. 5-86

    inner_frict_coeff = exp((25.8/log(Re_LOX)^2.58) - 2); % Beyvel, Eq. 5-87

    K_visc = (R * (inner_swirl_diam / 2)) / ((inner_num_inlets * ...
        (inner_inlet_diam / 2)^2) + (inner_frict_coeff / 2) * R * (R - (inner_swirl_diam / 2))); % Beyvel, Eq. 5-81

    visc_spray_angle = atan((2 * inner_disc_coeff * K_visc) / sqrt((1 + inner_S)^2 - ...
        (4 * inner_disc_coeff^2 * K_guess^2))); % Beyvel, Eq. 5-75

    % Calculate difference between desired and actual spray angle 
    angle_diff(lcv) = abs((visc_spray_angle - spray_angle) / spray_angle);
    lcv = lcv + 1;
end

% Inner Element constants and assumptions
inner_wall_thck = ; % [] wall thickness of the inner element
m_dot_FUEL = ; % [] mass flow of fuel
delta_p_FUEL ; % [] pressure drop of fuel 
FUEL_dens = ; % [] fuel density
outer_num_inlets = ; % [] number of inlets of the outer element
kin_visc_FUEL = ; % [] kinematic viscosity of the fuel
vortex_gap = .01181; % [] permitted gas vortex radius, Bazarov page 76

lcv = 0;
tf = []; % creating empty matrix
for tf_guess = 0:.01:.20
    outer_swirl_diam = (inner_swirl_diam + ...
        (2 * inner_wall_thick) + (2 * vortex_gap) + (2 * tf_guess)); % [] Outer swirl outlet diameter
    effective_flow_area = m_dot_FUEL / sqrt(2 * FUEl_dens * delta_p_FUEL); % [] Beyvel, Eq. 5-63

    outlet_area = (pi / 4) * outer_swirl_diam ^ 2; % [] 

    outer_disc_coeff = effective_flow_area / outlet_area; % [N/A]

    outer_filling_eff = fzero(@(outer_filling_eff) ...
        (outer_filling_eff * sqrt(outer_filling_eff / ...
        (2 - outer_filling_eff)) - outer_disc_coeff), .3); % [N/A] Beyvel, Eq. 5-66
    
    outer_K = (((1 - outer_filling_eff) * sqrt(2)) ...
        / (outer_filling_eff * sqrt(outer_filling_eff))); % [N/A] Beyvel, Eq. 5-65

     outer_S = fzero(@(outer_S) ((sqrt(1-(outer_disc_coeff^2 * outer_K^2))) - ...
        (outer_S * sqrt(outer_S^2-(outer_disc_coeff^2 * outer_K^2))) - ...
        (outer_disc_coeff^2 * outer_K^2 * log((1 + sqrt(1 - outer_disc_coeff^2 * outer_K^2))/(S + ...
        sqrt(outer_S^2 - outer_disc_coeff^2 * outer_K^2)))) - outer_disc_coeff), .8); % [] Justify this

     outer_spray_angle = atan((2 * outer_disc_coeff * outer_K) / (sqrt((1 + inner_S)^2 - ...
        (4 * inner_disc_coeff^2 * K_guess^2)))); % [] Beyvel Eq. 5-80

     outer_tf = (.5 * outer_swirl_diam) - (outer_S * outer_swirl_diam); % [] outer film thickness

     tf_diff(lcv) = abs((outer_tf - tf_guess) / tf_guess);
     lcv = lcv + 1;
end 




























