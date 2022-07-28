clc;
clear all; 

%{ 
In this version, the Beyvel and Orzechowski method is followed for the 
inner element, and the method outlined in Bazarov's "design and dynamics
of jet and swirl injectors" is followed for the outer element
%}

% Global Constants
gravity = 32.3; % [ft/s^2] Acceleration due to gravity

% Inner Element constants and assumptions
m_dot_LOX = 1; % [lbm/s] Oxidizer mass flow  
delta_p_LOX = 14400; % [lbf/ft^2] Oxidizer pressure drop 
LOX_dens = 71.168; % [lb/ft^3] LOX density 
inner_num_inlets = 3; % [N/A] number of tangential inlets 
spray_angle = 60; % [degrees] Desired spray half angle 
kin_visc_LOX= 2.362 * 10^-6; % [ft^2/s] kinematic viscosity of LOX
K_guess = 3.5;
inner_wall_thck = .005104; % [ft] wall thickness of the inner element, ~ 1/16"
coeff_nozzle_open = 3; % coefficient of nozzle opening, from Beyvel pg. 263 
    % should be (2-5), Bazarov says 3x

% Initial Calculations

lcv = 1;
matrix = [];

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
    
    R = (coeff_nozzle_open * inner_swirl_diam) / 2;

    inner_nozzle_length = 2 * inner_swirl_diam; % Suggested by Bazarov pg. 76

    inner_chamber_length = 3 * R; % From Bazarov, should be (2-3) * R
    
    inner_inlet_diam = sqrt((2 * R * inner_swirl_diam) / (inner_num_inlets * K_guess)); % Beyvel, Eq. 5-84
    
    inner_inlet_length = (inner_inlet_diam / 2) * 3; % Suggested by Bazarov to be (3-4) * inlet radius

    Re_LOX = (4 * m_dot_LOX)/ (pi * LOX_dens * ...
        kin_visc_LOX * sqrt(inner_num_inlets) * inner_inlet_diam); % Beyvel, Eq. 5-86
    
    inner_frict_coeff = exp((25.8/log(Re_LOX)^2.58) - 2); % Beyvel, Eq. 5-87
    
    K_visc = (R * (inner_swirl_diam / 2)) / ((inner_num_inlets * ...
        (inner_inlet_diam / 2)^2) + (inner_frict_coeff / 2) * R * (R - (inner_swirl_diam / 2))); % Beyvel, Eq. 5-81
    
    inner_visc_filling_eff = fzero(@(inner_visc_filling_eff) K_visc - (((1-inner_visc_filling_eff) * sqrt(2)) ...
        / (inner_visc_filling_eff * sqrt(inner_visc_filling_eff))), .4); % [N/A] Beyvel and Orzechowski, Eq. 5-65

    inner_visc_disc_coeff = inner_visc_filling_eff * sqrt(inner_visc_filling_eff / ...
        (2 - inner_visc_filling_eff));

    inner_visc_swirl_diam = sqrt(4 * m_dot_LOX / (pi * inner_visc_disc_coeff * ...
        sqrt(2 * LOX_dens * delta_p_LOX * gravity))); % Beyvel, Eq. 5-82

    K_visc = (2 * R * inner_visc_swirl_diam) / (inner_num_inlets * ...
        inner_inlet_diam^2); % Beyvel 5-83

    visc_spray_angle = atand((2 * inner_visc_disc_coeff * K_visc) / sqrt((1 + inner_S)^2 - ...
    (4 * inner_visc_disc_coeff^2 * K_guess^2))); % Beyvel, Eq. 5-80

    external_nozzle_diam = inner_swirl_diam + (2 * inner_wall_thck); % rough guess, wall = 1/16"

    lcv = lcv + 1;
    matrix(lcv,1) = K_guess;
    matrix(lcv,2) = K_visc;
    matrix(lcv,3) = inner_S;
    matrix(lcv,4) = visc_spray_angle;
    matrix(lcv,5) = inner_disc_coeff;

    if (K_guess - K_visc) / K_visc < .04
        break
    else
        K_guess = K_visc;
    end
end
























