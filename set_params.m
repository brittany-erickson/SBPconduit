function p = set_params(scale,flare)
%%
%% DEFINE UNSCALED:
    p.L = 1500;             % (km) Default 3000
    p.g = 10;              % m/s^2
    p.rho_tilde = 3000;     % kg/(m^3)
    p.T = 1100;             % K
    p.R = 462;              % J/kg^{-1} K
    
    p.nt = 0.02;            % unitless, default 0.01; % try 0.6 and L = 2000 and shoudl be gas everywhere and all negative e-vales. try larger valkue like 0.02 or 0.05
    p.eta = 0*10^4;           % 10^3 (Pa-s) viscosity
    p.Pc = 0;%1000;            % Pa between 10^3 to 10^6
    p.Patm = 1e5;           % atmospheric pressure (Pa)
    p.Kliq = 1e10;       % liquid compressibility (Pa)
    p.sw = 4e-6;            % 1/(sqrt(Pa)) 
    p.f0 = 0.001; %.001;    % kg/(s^2 m^2). f corresponds to drag from shear stress on conduit walls. 
    
    p.phiscale = 0.01;      % small values mean transition sharper.  If 0 does a step function
    p.tau = 1; %.01;      % s
    
    p.v_exit = 100;         % optional velocity BC (m/s)
  
    p.nsigma1 = 0.05;       % background gas mass fraction at base (not pressure dependent)
    p.nsigma2 = 0.0025;       % background gas mass fraction at top
    p.zDn =500;             % (m) gas frac will start to increase at zDn m below the surface (z = p.L-zD). this is the depth at which the nonequil gas frac is 1% greater than sigma1. 
    
    if p.nsigma2>p.nsigma1
        p.zstar_n = (-p.zDn)/atanh(((101/100)*p.nsigma1 - p.nsigma2)/(p.nsigma2-p.nsigma1));
    else
        p.zstar_n  = p.zDn;
    end
    
    if flare
       p.r1  = 5;           % (m)  conduit radius from z = 0 (at depth) to approxmately z = p.L-zD
        p.zD =300;         % (m) conduit will start to flare at zD m below the surface (z = p.L-zD). this is the depth at which the radius is 1% greater than r1. 
        p.r2 = 10;          % (m) conduit radius at Earth surface (z = L)
        p.zstar = (-p.zD)/atanh(((101/100)*p.r1 - p.r2)/(p.r2-p.r1));

    else
        
        p.r0 = 10;          % (m)  conduit radius at z = 0
        p.rL = 10;          % (m)  conduit radius at z = L
    end
    
%params for magmastatic code
p.GasExsolution=(p.nt/p.sw)^2;
p.GasLookup=0;
p.rescale=0;

if scale

    p.L = p.L/1000;                     % (km)
    p.g = p.g*1e3;                      % 10^{-3)m/s^2
    p.rho_tilde = p.rho_tilde/1e6;      % 10^6 kg/(m^3)
    
    p.eta = p.eta/10^3;                 % (10^3 Pa-s) viscosity
    p.Pc = p.Pc/10^6;                   % 10^6 Pa
    p.Kliq = p.Kliq/10^6;               % 10^6 Pa
    p.sw = p.sw*1e3;                    % 1/(10^3 sqrt(Pa)) 
    p.f0 = p.f0*1e-3;                   % 10^3 kg/(s^2 m^2). f corresponds to drag from shear stress on conduit walls. 
    p.tau = p.tau/1000;                 % 10^3 s
  
    p.z_c = p.z_c/1000;                 % km - depth at which coeff. are evaluated for constant coefficient case.

end
