%{
Purdue Space Program - Liquids
Research and Development - coaxial bi-swirl sizing code
Created by: Jacob Bell
In this early version, a  mix of both method's presented by Beyvel and Bazarov are used. 
Essentially just means that forms of
certain equations are expressed in varying ways (see document on the
drive for direct equation comparison). I'd like to emphasize the point that
many things are assumed initially, including the radius of swirling
"coefficient of nozzle opening" and chamber length.
%}
clc;
clear; 

% Global Constants
gravity = 32.3; % [ft/s^2] Acceleration due to gravity

%% Inner Element constants and assumptions
m_dot_LOX = 1; % [lbm/s] Oxidizer mass flow  
delta_p_LOX = 10800; % [lbf/ft^2] Oxidizer pressure drop 
LOX_dens = 71.168; % [lb/ft^3] LOX density 
in_inlets = 3; % [N/A] number of tangential inlets 
spray_angle = 30; % [degrees] Desired spray half angle 
kin_visc_LOX= 2.362 * 10^-6; % [ft^2/s] kinematic viscosity of LOX
K_guess = 3;
inner_wall_thck = .005104; % [ft] wall thickness of the inner element, ~ 1/16"
coeff_nozzle_open = 2; % coefficient of nozzle opening, from Beyvel pg. 263 
    % should be (2-5), Bazarov says 3x. This value is the ratio of "swirl
    % arm"/nozzle radius

%% Inner Loop
lcv = 1;
while lcv < 100
    in_fill_eff = fzero(@(in_fill_eff) K_guess - (((1-in_fill_eff) * sqrt(2)) ...
        / (in_fill_eff * sqrt(in_fill_eff))), .4); % [N/A] Beyvel and Orzechowski, Eq. 5-65
    
    in_disc_coeff = in_fill_eff * sqrt(in_fill_eff / (2 - in_fill_eff)); % [N/A] Beyvel and Orzechowski, Eq. 5-66
    
    fcn = @(in_s) ((sqrt(1 - in_disc_coeff^2 * K_guess^2)) ...
    - (in_s * sqrt(in_s^2 - in_disc_coeff^2 * K_guess^2)) ...
    - (in_disc_coeff^2 * K_guess^2 * log((1 ...
    + sqrt(1 - in_disc_coeff^2 * K_guess^2)) ...
    / (in_s + sqrt(in_s^2 - in_disc_coeff^2 ...
    * K_guess^2)))) - in_disc_coeff); % Beyvel Eq. 5-58

    in_s = fzero(fcn, .8);

    in_swirl_diam = sqrt(4 * m_dot_LOX / (pi * in_disc_coeff * ...
        sqrt(2 * LOX_dens * delta_p_LOX * gravity))); % Beyvel, Eq. 5-82
    
    in_R = (coeff_nozzle_open * in_swirl_diam) / 2; % "swirl arm", center axis to inlet
    
    in_inlet_diam = sqrt((2 * in_R * in_swirl_diam) / (in_inlets * K_guess)); % Beyvel, Eq. 5-84

    Re_LOX = (4 * m_dot_LOX)/ (pi * LOX_dens * ...
        kin_visc_LOX * sqrt(in_inlets) * in_inlet_diam); % Beyvel, Eq. 5-86
    
    in_frict_coeff = exp((25.8/log(Re_LOX)^2.58) - 2); % inlet hydraulic loss Coefficient Beyvel, Eq. 5-87
    
    K_frict = (in_R * (in_swirl_diam / 2)) / ((in_inlets * ...
        (in_inlet_diam / 2)^2) + (in_frict_coeff / 2) * ...
        in_R * (in_R - (in_swirl_diam / 2))); % K value taking hydraulic loss into account Beyvel, Eq. 5-81

    eq_filling_eff = fzero(@(eq_filling_eff) K_frict - (((1-eq_filling_eff) * sqrt(2)) ...
        / (eq_filling_eff * sqrt(eq_filling_eff))), .4); % [N/A] Beyvel and Orzechowski, Eq. 5-65

    eq_disc_coeff = eq_filling_eff * sqrt(eq_filling_eff / ...
        (2 - eq_filling_eff)); % temporary discharge coeff from K_frict value

    in_visc_disc_coeff = eq_disc_coeff / (sqrt(1 + ...
        eq_disc_coeff^2 * (K_guess^2 / coeff_nozzle_open^2))); % Bazarov Eq. 99. Takes into account angular momentum losses

    in_visc_swirl_diam = sqrt(4 * m_dot_LOX / (pi * in_visc_disc_coeff * ...
        sqrt(2 * LOX_dens * delta_p_LOX * gravity))); % Beyvel, Eq. 5-82

    in_R = (coeff_nozzle_open * in_visc_swirl_diam) / 2;

    K_visc = (2 * in_R * in_visc_swirl_diam) / (in_inlets * ...
        in_inlet_diam^2); % Beyvel 5-83

    visc_spray_angle = atand((2 * in_disc_coeff * K_visc) / sqrt((1 + in_s)^2 - ...
    (4 * in_visc_disc_coeff^2 * K_guess^2))); % Beyvel, Eq. 5-80\

    lcv = lcv + 1;
    if ((K_guess - K_visc) / K_visc < .03)
        K_final = K_visc;
        in_diam_final = in_visc_swirl_diam;
        ex_nozzle_diam = in_diam_final + (2 * inner_wall_thck); % rough guess, wall = 1/16"
        break
    else
        K_guess = K_visc;
    end
end

%% Inner swirler final parameter calcs
in_nozzle_length = 2 * in_swirl_diam; % Suggested by Bazarov pg. 76
in_cham_length = 3 * in_R; % From Bazarov, should be (2-3) * R
in_inlet_length = (in_inlet_diam / 2) * 3; % Suggested by Bazarov to be (3-4) * inlet radius
in_cham_diam = (2 * in_R) + in_inlet_diam;

% TIME FOR PART TWO!!!!!!

%% Outer Element constants and assumptions
m_dot_FUEL = .313; % [lbm/s] Fuel mass flow  
delta_p_FUEL = 10800; % [lbf/ft^2] Fuel pressure drop 
FUEL_dens = 27.0; % [lb/ft^3] Fuel density 
out_num_inlets = 4; % [N/A] number of tangential inlets 
kin_visc_FUEL= 2.362 * 10^-6; % [ft^2/s] kinematic viscosity of Fuel
coeff_nozzle_open = 3.5; % coefficient of nozzle opening, from Beyvel pg. 263 
    % should be (2-5), Bazarov says 3x (for closed). This value is the ratio of "swirl
    % arm"/nozzle radius
permitted_vortex_rad = (ex_nozzle_diam / 2) + .00098425; % value comes from bazarov, Pg. 76 is .3mm

%% INITIAL OUTER VALUES using minimum values for outlet area 
film_thick_guess = .00098425; % same as gap, just a guess
out_swirl_diam = 2 * permitted_vortex_rad + 2 * film_thick_guess; % first estimation of outer diameter
effective_flow_area = m_dot_FUEL / sqrt(2 * FUEL_dens * delta_p_FUEL * gravity); % [] Beyvel, Eq. 5-63
outlet_area = (pi / 4) * out_swirl_diam^2; % [] 
out_disc_coeff = effective_flow_area / outlet_area; % [N/A]

lcv = 1;

%% Outer Loop
while lcv < 100
    fcn = @(out_fill_eff) (out_fill_eff * sqrt(out_fill_eff ...
        / (2 - out_fill_eff)) - out_disc_coeff); % Beyvel Eq. 5-66

    out_fill_eff = fzero(fcn, .7);

    outer_K = (((1 - out_fill_eff) * sqrt(2)) ...
        / (out_fill_eff * sqrt(out_fill_eff))); % [N/A] Beyvel, Eq. 5-65

    out_disc_coeff = out_fill_eff * sqrt(out_fill_eff / (2 - out_fill_eff));

    fcn = @(outer_S) ((sqrt(abs(1 - out_disc_coeff^2 * outer_K^2))) ...
    - (outer_S * sqrt(abs(outer_S^2 - out_disc_coeff^2 * outer_K^2))) ...
    - (out_disc_coeff^2 * outer_K^2 * log((1 ...
    + sqrt(abs(1 - out_disc_coeff^2 * outer_K^2))) ...
    / (outer_S + sqrt(abs(in_s^2 - out_disc_coeff^2 ...
    * outer_K^2))))) - out_disc_coeff); % Beyvel Eq. 5-58

    outer_S = fzero(fcn, .7);

    film_thick = .5 * (out_swirl_diam - (outer_S * out_swirl_diam));

    out_swirl_diam_new = (2 * permitted_vortex_rad) + (2 * film_thick);

     if abs((out_swirl_diam_new - out_swirl_diam) / out_swirl_diam) < .005
         break
     else
         out_swirl_diam = out_swirl_diam_new;
     end
     lcv = lcv + 1;
end

%% Outer element final parameters
outer_R = (out_swirl_diam * coeff_nozzle_open) / 2;
out_inlet_diam = sqrt((2 * outer_R * out_swirl_diam) / (out_num_inlets * outer_K)); % Beyvel, Eq. 5-84
outer_nozzle_length = 2 * out_swirl_diam; % Suggested by Bazarov pg. 76
out_cham_length = outer_R; % can really be anything reasonable
out_inlet_length = (out_inlet_diam / 2) * 3; % Suggested by Bazarov to be (3-4) * inlet radius
outer_chamber_diam = (2 * outer_R) + out_inlet_diam;
outer_check = (2 * outer_R * out_swirl_diam) / (out_num_inlets * out_inlet_diam^2);

%% Inner swirler output table
inner_params = ["Inner Swirl Diam."; "Inner Inlet Diam."; "Inner Chamber Lg."; "Inner Inlet Lg."; "Inner Chamber diam."; "Discharge Coef."];
inner_value = [in_diam_final * 12; in_inlet_diam * 12; in_cham_length * 12; in_inlet_length * 12; in_cham_diam * 12; in_visc_disc_coeff];
units = ["[in]"; "[in]"; "[in]"; "[in]"; "[in]"; "[N/A]"];
inner_output = table(inner_params, inner_value, units);

%% Outer swirler output table
outer_params = ["Outer Swirl Diam."; "Outer Inlet Diam."; "Outer Chamber Lg."; "Outer Inlet Lg."; "Outer Chamber diam."];
outer_value = [out_swirl_diam_new * 12; out_inlet_diam * 12; out_cham_length * 12; out_inlet_length * 12; outer_chamber_diam * 12];
outer_units = ["[in]"; "[in]"; "[in]"; "[in]"; "[in]"];
outer_output = table(outer_params, outer_value, outer_units);

















