
% Solve eigenvalue/vector problem for
% dU_dt + M(z)dU_dz = B(z)U, where U = [u P rho n]'
% diagonalizes to dW_dt + L(z)dW_dz = CW
% with L = inv(P)*M*P, C = inv(P)*B*P - L*inv(P)*Pz and
% boundary conditions on characteristics.


clear all
close all

% SET ORDER OF SPATIAL ACCURACY
order = 4; 

% SET IF CONSTANT COEFFICIENTS
constcoeff=false;

% Specify conduit geometry - to flare or not. 
flare = false;

% Get parameters - scaled or not
scale = false;
p = set_params(scale, flare);

% Set physical boundary conditions (will be translated to characteristic b.c.)
bc_L = 'pres';              % options: 'vel' or 'pres' - top boundary condition, data set to 0.
bc_0 = 'vel, dens, n';      % options: 'vel, dens, n' or 'press, dens, n' - bottom boundary condition, data set to 0.

bg_type =  'explosive';     %options 'explosive', 'magmastatic' or 'custom' 


% Set grid spacing (m) at which to compute background state, and color for
% plotting eigenvalues.
DZ=5; %grid resolution 
col = 'ko';
z = 0:DZ:p.L;
N = length(z);

%SOLVE FOR BACKGROUND STATE (vel, pressure), and obtain dependent fields (rho, n)
switch bg_type
    case 'custom'
         
        vdepth = 10; vsurf = 100;
        vbar = [vdepth + ((vsurf-vdepth)/2)*(8/p.L^3).*(z - p.L/2).^3 + (vsurf-vdepth)/2]'; 
        pdepth = 3e7; psurf = 3e6;
        Pbar = [pdepth + ((psurf-pdepth)/2)*(8/p.L^3).*(z - p.L/2).^3 + (psurf-pdepth)/2]';
        bg = get_bg_bg(z,vbar,Pbar,p,flare);
        bg.vbar = vbar;
        bg.Pbar = Pbar;
        
    case 'explosive'
        % uncomment one below
        bg = solve_bg_state(DZ,1e-3,0,100,10000,'choked_flow',p, flare); bc_L = 'vel';  % choked flow.
        % bg = solve_bg_state(DZ,1e-9,0,100,10000,25,p, flare);                            % a specified exit velocity.

%keyboard
    case 'magmastatic'
        bg = solve_bg_state(DZ,1e-4,0,2*(p.Patm+p.Pc+p.rho_tilde*p.g*p.L),10000,'magmastatic',p, flare);                            % specified exit velocity = 0.
end

% Evaluate if constant coefficient case.
if constcoeff
    % specify where to evaluate background state. E.g. Nev = (N-1)/2 +1 is
    % middle of conduit, Nev = N is Earth surface. 
    Nev = (N-1)/2 + 1;
    vbar = bg.vbar(Nev)*ones(N,1); 
    Pbar = bg.Pbar(Nev)*ones(N,1);
    bg = get_bg_bg(z(Nev)*ones(N,1),vbar,Pbar,p,flare);
    bg.vbar = vbar;
    bg.Pbar = Pbar;
    bg.z = z;

end

% Plot background state
figure(1)
subplot(1,3,1)        
hold on
plot(bg.vbar,z)
subplot(1,3,2)
hold on
plot(1./bg.Kbar .*bg.dPbar_dz - 1./bg.rhobar .*bg.drhobar_dz,z)
subplot(1,3,3)
hold on
plot(bg.cbar,z)
%plot(bg.bbar.*bg.dPbar_dz,z)
plot(bg.ceqbar,z)
%plot(z,bg.vbar,'bo',z,bg.ceqbar,'ro')

%pause

% GET/FORM OPERATORS
[z,dz,D,Hinv00,A,H] = grid_operators(N-1,p.L,order); 
Hinv = inv(H);
Hinv4 = blkdiag(Hinv,Hinv,Hinv,Hinv);
D4 = blkdiag(D,D,D,D);
H4 = kron(speye(4,4),H);

% Characteristic speeds/eigenvalues
lambda1 = bg.vbar; lambda2 = bg.vbar - bg.cbar; lambda3 = bg.vbar + bg.cbar; lambda4 = bg.vbar;

L4 = blkdiag(diag(lambda1),diag(lambda2),diag(lambda3),diag(lambda4));

[B, P, Pinv, Pz, Pinv_Pz] = form_matrices(bg, p, D, N);

% Quadrature for computing integrals used in continuous spectrum
Im1 = (dz/2)*(1/lambda1(1)) + sum(1./(lambda1(2:end-1)))*dz + (dz/2)*(1./lambda1(end));
Im2 = (dz/2)*(1/lambda2(1)) + sum(1./(lambda2(2:end-1)))*dz + (dz/2)*(1./lambda2(end));
Im3 = ((dz/2)*(1/lambda3(1)) + sum(1./(lambda3(2:end-1)))*dz + (dz/2)*(1./lambda3(end)));
Im4 = ((dz/2)*(1/lambda4(1)) + sum(1./(lambda4(2:end-1)))*dz + (dz/2)*(1./lambda4(end)));

% Penalty parameters 
sigma1 = -lambda1(1); sigma2 = lambda2(end); sigma3 = -lambda3(1); sigma4 = -lambda4(1);

% Determine reflection coefficients
switch bc_0
    case('vel, dens, n')
        R1 = -2;   
        R3 = 1;
        R4 = 0;
    case('pres, dens, n') 
        R1 = 0;
        R3 = -1;
        R4 = 0;
end


switch bc_L
    case('vel')
        R = 1;  %R such that v = 0    
    case('pres')
        R = -1;  %R such that p = 0    
end

% Build penalty matrix
S4 = 0.*speye(4*N,4*N);
S4(1,1) = sigma1*Hinv00;
S4(1,N+1) = -sigma1*Hinv00*R1;
S4(2*N+1,2*N+1) = sigma3*Hinv00;
S4(2*N+1,N+1) = -sigma3*Hinv00*R3;
S4(3*N+1,3*N+1) = sigma4*Hinv00;
S4(3*N+1,N+1) = -sigma4*Hinv00*R4;
S4(2*N,2*N) = sigma2*Hinv00;
S4(2*N,3*N) = -sigma2*Hinv00*R;

% Build numerical approximation to operator L(z)d_dz
DLAM = 0.5*(D4*L4 + L4*D4)  - 0.5*diag(D4*[lambda1;lambda2;lambda3;lambda4]);     % Skew-symmetric SBP
DLAM2 = L4*D4;                                                                    % alternative option ("standard" SBP) but doesn't have stability proof. 

% Build numerical operator C
C = 1.*Pinv*B*P - L4*Pinv*Pz;                                                 

% Build entire matrix for right hand side. 
MAT = -DLAM  + S4 + C;
e = eig(MAT);

% Plot numerical eigvenvalues
figure(2)
max(real(e))
hold on
plot(real(e),imag(e),'o')% col')

%if strcmp(col, 'ko')
 %   mm = -10.*(N-1)/1:1:10.*(N-1)/1;
 %   sn = (2*pi.*mm.*1i + log(R*R3))/(Im3-Im2);
 %   plot(real(sn),imag(sn),'b.')
 %   hold on
%end


return


