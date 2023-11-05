function [B, P, Pinv, Pz] = form_matrices(bg, p, D, N)
%[B, P, Pinv, Pz, Pinv_Pz, DP] = form_matrices(bg, p, D, N)

b11 = -diag(bg.dvbar_dz + bg.Fv./bg.rhobar);
b12 = -diag(bg.vbar.*bg.dvbar_dz + p.g)./bg.rhobar;
b13 = diag(-bg.Fp./bg.rhobar);
b14 = diag(-bg.Fn./bg.rhobar);

%b21 = -diag((1./bg.ceqbar.^2).*bg.dPbar_dz + bg.rhobar.*bg.dA_dz./bg.A);
b21 = -diag(bg.drhobar_dz + bg.rhobar.*bg.dA_dz./bg.A);
b22 =  -diag((bg.dvbar_dz.*bg.A + bg.vbar.*bg.dA_dz)./bg.A);
b23 = 0*diag(ones(N,1));
b24 = 0*diag(ones(N,1));

b31 = -diag(bg.dPbar_dz + bg.Kbar.*bg.dA_dz./bg.A);
b32 = -diag(bg.Kbar.*(bg.dvbar_dz.*bg.A + bg.vbar.*bg.dA_dz)./(bg.rhobar.*bg.A));
b33 = -diag(bg.Kbar.*bg.abar.*bg.bbar./p.tau);
b34 = -diag(bg.Kbar.*bg.abar./p.tau);

b41 = diag(bg.bbar.*bg.dPbar_dz);
b42 =  0*diag(ones(N,1)); 
b43 = -diag(bg.bbar./p.tau);
b44 = -diag(ones(N,1)./p.tau);

B = [b11 b12 b13 b14;b21 b22 b23 b24;b31 b32 b33 b34;b41 b42 b43 b44];
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p11 = 0*diag(ones(N,1));
p12 = diag(-1./(2.*bg.zbar));
p13 = diag(1./(2.*bg.zbar));
p14 = 0*diag(ones(N,1));

p21 = diag(1./(2.*bg.cbar.^2));
p22 = diag(1./(2.*bg.cbar.^2));
p23 = diag(1./(2.*bg.cbar.^2));
p24 =  0*diag(ones(N,1));

p31 = 0*diag(ones(N,1));
p32 = 1/2*diag(ones(N,1));
p33 = 1/2*diag(ones(N,1));
p34 = 0*diag(ones(N,1));

p41 = 0*diag(ones(N,1));
p42 = 0*diag(ones(N,1));
p43 = 0*diag(ones(N,1));
p44 = diag(sqrt(1./(2.*bg.Kbar)));

P = [p11 p12 p13 p14; p21 p22 p23 p24; p31 p32 p33 p34; p41 p42 p43 p44];

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pinv11 = 0*diag(ones(N,1));
pinv12 = diag(2*bg.cbar.^2);
pinv13 = -2*diag(ones(N,1));
pinv14 = 0*diag(ones(N,1));

pinv21 = -1.*diag(bg.zbar);
pinv22 =  0*diag(ones(N,1));
pinv23 = diag(ones(N,1));
pinv24 =  0*diag(ones(N,1));

pinv31 = diag(bg.zbar);
pinv32 = 0*diag(ones(N,1));
pinv33 = diag(ones(N,1));
pinv34 = 0*diag(ones(N,1));

pinv41 = 0*diag(ones(N,1));
pinv42 = 0*diag(ones(N,1));
pinv43 = 0*diag(ones(N,1));
pinv44 = diag(sqrt(2.*bg.Kbar));

Pinv = [pinv11 pinv12 pinv13 pinv14; pinv21 pinv22 pinv23 pinv24; pinv31 pinv32 pinv33 pinv34; pinv41 pinv42 pinv43 pinv44];
     
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical calculation of P_z:
pz11 = 0*diag(ones(N,1));
pz12 = diag(D*(-1./(2.*bg.zbar)));
pz13 = diag(D*(1./(2.*bg.zbar)));
pz14 = 0*diag(ones(N,1));

pz21 = diag(D*(1./(2.*bg.cbar.^2)));
pz22 = diag(D*(1./(2.*bg.cbar.^2)));
pz23 = diag(D*(1./(2.*bg.cbar.^2)));
pz24 =  0*diag(ones(N,1));

pz31 = 0*diag(ones(N,1));
pz32 = 0*diag(ones(N,1));
pz33 = 0*diag(ones(N,1));
pz34 = 0*diag(ones(N,1));

pz41 = 0*diag(ones(N,1));
pz42 = 0*diag(ones(N,1));
pz43 = 0*diag(ones(N,1));
pz44 = diag(D*(sqrt(1./(2.*bg.Kbar))));

DP = [pz11 pz12 pz13 pz14; pz21 pz22 pz23 pz24; pz31 pz32 pz33 pz34; pz41 pz42 pz43 pz44];

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analytic calculation of P_z
pz11 = 0*diag(ones(N,1));
%%old version (prior to Jasper's work!)
pz12 = diag((bg.dzbar_dz./(2*bg.zbar.^2)));
pz13 = diag((-bg.dzbar_dz./(2*bg.zbar.^2)));
%%new version
% pz12 = diag(-bg.d2Zinv_dz);
% pz13 = diag(bg.d2Zinv_dz);

pz14 = 0*diag(ones(N,1));

%%old version
pz21 = diag((-bg.dcbar_dz./(bg.cbar.^3)));
pz22 = diag((-bg.dcbar_dz./(bg.cbar.^3)));
pz23 = diag((-bg.dcbar_dz./(bg.cbar.^3)));
%%new version
% pz21 = diag(bg.d2C2inv_dz);
% pz22 = diag(bg.d2C2inv_dz);
% pz23 = diag(bg.d2C2inv_dz);

pz24 =  0*diag(ones(N,1));

pz31 = 0*diag(ones(N,1));
pz32 = 0*diag(ones(N,1));
pz33 = 0*diag(ones(N,1));
pz34 = 0*diag(ones(N,1));

pz41 = 0*diag(ones(N,1));
pz42 = 0*diag(ones(N,1));
pz43 = 0*diag(ones(N,1));
%%old version
pz44 = diag((0.5.*(1./(2.*bg.Kbar)).^(-0.5)) .* -2.*bg.dkbar_dz./(4.*bg.Kbar.^2));%diag(D*(sqrt(1./(2.*bg.Kbar))));
%%new version
% pz44 = diag(bg.dsqrt2K_dz);

Pz = [pz11 pz12 pz13 pz14; pz21 pz22 pz23 pz24; pz31 pz32 pz33 pz34; pz41 pz42 pz43 pz44];



% pinv_pz11 = diag((-2./bg.cbar).*(bg.dcbar_dz));
% pinv_pz12 = diag((-2./bg.cbar).*(bg.dcbar_dz));
% pinv_pz13 = diag((-2./bg.cbar).*(bg.dcbar_dz));
% pinv_pz14 = 0*diag(ones(N,1));
% 
% pinv_pz21 = 0*diag(ones(N,1));
% pinv_pz22 = (-1/2).*diag((bg.dzbar_dz)./bg.zbar);
% pinv_pz23 = (1/2).*diag((bg.dzbar_dz)./bg.zbar);
% pinv_pz24 =  0*diag(ones(N,1));
% 
% pinv_pz31 = 0*diag(ones(N,1));
% pinv_pz32 = (1/2).*diag((bg.dzbar_dz)./bg.zbar);
% pinv_pz33 = (-1/2).*diag((bg.dzbar_dz)./bg.zbar);
% pinv_pz34 = 0*diag(ones(N,1));
% 
% pinv_pz41 = 0*diag(ones(N,1));
% pinv_pz42 = 0*diag(ones(N,1));
% pinv_pz43 = 0*diag(ones(N,1));
% pinv_pz44 = (-1/2).*diag((sqrt(2.*bg.Kbar)./bg.Kbar.^2).*(bg.dkbar_dz));
% 
% Pinv_Pz = [pinv_pz11 pinv_pz12 pinv_pz13 pinv_pz14; pinv_pz21 pinv_pz22 pinv_pz23 pinv_pz24; pinv_pz31 pinv_pz32 pinv_pz33 pinv_pz34; pinv_pz41 pinv_pz42 pinv_pz43 pinv_pz44];

