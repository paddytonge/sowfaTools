% Constants
k = 0.4;                 % von Karman constant
rho = 1.225;
thetav = 300;
g = 9.81;                % gravitational acceleration (m/s^2)
z = 10;                  % Height at which wind speed is measured (m)
z_h = 80;
z0 = 0.001;                % Surface roughness length (m)
U_obs = u1(1);           % Observed wind speed at height z (m/s)
gamma1 = 16;             % Coefficient for unstable conditions
beta1 = 5;               % Coefficient for stable conditions

%Heat flux if a cooling rate applied
% cooling_rate = -1.38888e-4;  % Cooling rate (K/s)
% cp = 1005;               % Specific heat capacity of air (J/(kg·K))
% rho = 1.225;             % Density of air at sea level (kg/m^3)
% flux = cooling_rate * ((cp * rho)/g);  % Sensible heat flux (W/m^2)

%Heat flux applied
flux = 0.04;

%Initial u_star guess
tau_u = uw1(1);          % Assuming first index corresponds to the surface
tau_v = vw1(1);          % Assuming first index corresponds to the surface
u_star_guess = sqrt((tau_u^2 + tau_v^2)/rho); % Initial guess for u_star (m/s)
tolerance = 1e-4;        % Convergence tolerance

% Main execution
u_star = solve_u_star(k, z, z0, U_obs, gamma1, beta1, u_star_guess, flux, g);
fprintf('Calculated friction velocity u_* = %f m/s\n', u_star);

zABL = 970;
cell = round(zABL*(length(heights)/1280));
Deltat = t1(cell)-t1(1);
DeltaU = u1(cell)-u1(1);
DeltaV = v1(cell)-v1(1);
Deltaz = heights(cell)-0;
    
L_final = (-thetav*((u_star)^3))/(k*g*flux);
fprintf('Calculated Obukhov Length L = %f m/s\n', L_final);

zeta_final = z_h/L_final;
fprintf('Calculated Stability Parameter zeta = %f m/s\n', zeta_final);
    
Ri = ((g/thetav)*((Deltat)*(Deltaz)))/((DeltaU)^2+(DeltaV)^2);
fprintf('Calculated Bulk Reynolds Number Ri_B = %f m/s\n', Ri);

% Function to compute Psi_M based on stability parameter zeta
function psi = Psi_M(zeta, gamma1, beta1)
    if zeta < 0  % Unstable
        x = (1 - (gamma1 * zeta))^(1/4);
        psi = (2 * log((1 + x) / 2)) + (log((1 + x^2) / 2)) - (2 * atan(x)) + (pi / 2);
    elseif zeta > 0  % Stable
        psi = -beta1 * zeta;
    else  % Neutral
        psi = 0;
    end
end

% Iterative solver for u_star
function u_star = solve_u_star(k, z, z0, U_obs, gamma1, beta1, u_star_guess, flux, g)
    u_star = u_star_guess;
    max_iter = 1000;
    iter = 1;
    
    for iter = 1:max_iter
        L = -u_star^3 / (k * g * flux); % Calculate Monin-Obukhov length
        zeta = z / L;
        psi = Psi_M(zeta, gamma1, beta1);
        u_star = ((U_obs*k)/(log(z / z0)))+psi;  
    end
end
