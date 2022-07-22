clc; 
clear all; 

K_guess = 2.4;
inner_filling_eff = fzero(@(inner_filling_eff) (K_guess - ...
    (((1-inner_filling_eff) * sqrt(2)) / (inner_filling_eff * sqrt(inner_filling_eff)))), .4); % [N/A] Beyvel and Orzechowski, Eq. 5-65

inner_disc_coeff = inner_filling_eff * sqrt(inner_filling_eff / (2 - inner_filling_eff));

fcn = @(inner_S) ((sqrt(1 - inner_disc_coeff^2 * K_guess^2)) ...
    - (inner_S * sqrt(inner_S^2 - inner_disc_coeff^2 * K_guess^2)) ...
    - (inner_disc_coeff^2 * K_guess^2 * log((1 ...
    + sqrt(1 - inner_disc_coeff^2 * K_guess^2)) ...
    / (inner_S + sqrt(inner_S^2 - inner_disc_coeff^2 ...
    * K_guess^2)))) - inner_disc_coeff);
inner_S = fzero(fcn, 1);



