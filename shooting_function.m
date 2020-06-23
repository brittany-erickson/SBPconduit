function dy = shooting_function(z0,y,p,flare)


vbar = y(1);
Pbar = y(2);

F = get_bg_fields(z0,vbar,Pbar,p,flare);

dy = zeros(2,1);

rat = (vbar/F.ceqbar)^2;


if rat == 0      
    dy(1) = 0;
    dy(2) = -F.rhobar*p.g ;%(-F.rhobar*p.g - F.fbar + (F.rhobar*vbar^2/F.A)*F.dA_dz);
else
    dy(1) = (1/(1-rat))*(rat*(F.rhobar*p.g + F.fbar)/(F.rhobar*vbar) - (vbar/F.A)*F.dA_dz + F.abar*vbar*F.dsigma_dz);
    dy(2) = (1/(1-rat))*(-F.rhobar*p.g - F.fbar + (F.rhobar*vbar^2/F.A)*F.dA_dz - F.rhobar*F.abar*vbar^2*F.dsigma_dz);
end




