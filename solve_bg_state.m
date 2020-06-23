function bg = solve_bg_state(dz, tol, LEFT, RIGHT, Nmax, bc_type, p, flare)

% Nonlinear shooting method that returns fields that solve steady state
% background equations.
% i.e. solve dvbar_dz = f1(t,vbar,Pbar), dPbar_dz = f2(t,vbar,Pbar) on z in [0, L], with boundary conditions
% vbar(L) =v_exit or ceq_bar (choked flow), Pbar(0) =  p.rhobar*p.g*p.L + p.Pc
% iterate a max of Nmax times, over vbar(0) = v0^k, i.e. velocity at base
% of conduit, until err = |vbar(L)-v_exit| < tol.
% same is true for magmastatic, except we iterate over pressure to satisfy
% atmospheric pressure top BC
% [LEFT, RIGHT] is interval over which we bisect to find vbar(0).
% tol is error tolerance for for choked flow b.c.

% scale = false;
% p = set_params(scale);

m = (LEFT + RIGHT)/2; %initial guess is midpoint.

% set initial guesses for vbar and Pbar  at base of conduit:
if strcmp(bc_type,'magmastatic')
    vbar0 = 0;
elseif strcmp(bc_type,'velocity')
    vbar0 = 0;
else
    vbar0 = m;    
end

%get background fields for basal pressure 
F0 = get_bg_fields(0,0,p.rho_tilde*p.g*p.L + p.Pc + p.Patm, p, flare);
%guess pressure based on gas present (or not)
Pbar0 = F0.rhobar*p.g*p.L + p.Pc + p.Patm;

err = 1e12;     % initial error for err = abs(v_exit - vbar(L))
diff = 100;     % difference between current guess vbar0 and updated guess m

opt = odeset('reltol',1e-12,'abstol',1e-12); % tolerances for ODE solver
iter = 1; 

if strcmp(bc_type,'choked_flow')

    while err > tol && diff > tol && iter <= Nmax 

                Y0 = [vbar0; Pbar0];                                                % set initial values for ODE solver
                [z, Y] = ode45(@shooting_function,[0:dz:p.L], Y0,opt,p, flare);     % solve ODE
                vbar = Y(:,1); Pbar = Y(:,2);                                       % unpack solution
                F = get_bg_fields(z,vbar,Pbar,p, flare);                            % get fields evaluated at vbar, Pbar
            

                if isreal(F.ceqbar) && isreal(vbar) && F.ceqbar(end) >= vbar(end)  && abs(z(end) - p.L) < 1e-9  %then initial vel too low  
                    err  = F.ceqbar(end) - vbar(end);
                    LEFT = m;
                else
                    err = 1e12;
                    RIGHT = m;  % initial velocity too high
          
                end


                m = (LEFT + RIGHT)/2;
  
                if isreal(vbar0) && isreal(m) && abs(z(end)-p.L)<1e-9
                    diff = abs(vbar0-m);
                else
                    diff = 100;
                end
                
                disp('still working')
                
                vbar0 = m;
                iter = iter + 1;

    end
     if err < tol
         disp(['converged in ' num2str(iter) ' iterations'])
     else
         disp(['did not converge in ' num2str(Nmax) ' iterations'])
     end
    
else
    
     while err > tol && diff > tol && iter <= Nmax

                Y0 = [vbar0;Pbar0];                                             % set initial values for ODE solver
                [z, Y] = ode45(@shooting_function,[0:dz:p.L], Y0,opt,p,flare);  % solve ODE
                vbar = Y(:,1); Pbar = Y(:,2);                                   % unpack solution
                F = get_bg_fields(z,vbar,Pbar,p,flare);                         % get fields evaluated at vbar, Pbar

        if strcmp(bc_type,'magmastatic')
            %in case of magmastatic iterate over basal pressure
            
                if Pbar(end) < p.Patm
                    err = abs(p.Patm - Pbar(end)); 
                    LEFT = m;                        
                else
                    RIGHT = m;                                                  % initial velocity too high
                    err = abs(p.Patm - Pbar(end));%1e12;           
                end
                
                m = (LEFT + RIGHT)/2;
        else
            %otherwise, iterate over basal velocity
                if isreal(vbar) && vbar(end) < p.v_exit && abs(z(end)-p.L)<1e-9   % then initial vel too low
                        err = abs(p.v_exit - vbar(end));
                        LEFT = m;      
                else
                    RIGHT = m;                                                  % initial velocity too high
                    err = 1e12;           
                end
                
                m = (LEFT + RIGHT)/2;

                if isreal(vbar0) && isreal(m)
                    diff = abs(vbar0-m);
                else
                    diff = 100;
                end                
        end
                    
                disp('still working')
                
                if strcmp(bc_type,'magmastatic')
                    Pbar0 = m;
                else
                    vbar0 = m;
                end
                
                iter = iter + 1;

     end
     if err < tol
         disp(['converged in ' num2str(iter) ' iterations'])
     else
         disp(['did not converge in ' num2str(Nmax) ' iterations'])
     end
   
end

    
    
%output
bg = get_bg_fields(z,vbar,Pbar,p,flare);
bg.vbar = vbar;
bg.Pbar = Pbar;
bg.z = z; 
bg.p = p;




