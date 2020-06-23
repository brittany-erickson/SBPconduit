%Given a background state (vbar, Pbar), get_bg_fields() determines the remaining
%dependent variables in the background statevas well as descriptions of the
%conduit geometry (A, dA_dz).  


function fields = get_bg_fields(z,vbar,Pbar,p, flare)

N = length(vbar); 

%initialize outputs:
rhobar = zeros(N,1);
nbar = zeros(N,1);
Kbar = zeros(N,1);
Fv = zeros(N,1);
Fn = zeros(N,1);
Fp = zeros(N,1);
r = zeros(N,1);
rho_liq = zeros(N,1);
rhog = zeros(N,1);
dP_drho = zeros(N,1);
drho_dn = zeros(N,1);
cbar = zeros(N,1);
ceqbar = zeros(N,1);
abar = zeros(N,1);
phi = zeros(N,1);
fbar = zeros(N,1);
bbar = zeros(N,1); 
A = zeros(N,1);
dA_dz = zeros(N,1);

%initialize partial derivatives
dnbar_dP = zeros(N,1);
drho_dP = zeros(N,1);
drhobar_dPbar = zeros(N,1);
drho_liq_dP = zeros(N,1);
drhog_dP = zeros(N,1); 
dvbar_dz = zeros(N,1);
drhobar_dz = zeros(N,1);
dPbar_dz = zeros(N,1);
dnbar_dz = zeros(N,1);
dcbar_dz = zeros(N,1);
dkbar_dz = zeros(N,1);
dphi_dn = zeros(N,1);
dphi_dp = zeros(N,1);

dzbar_dz = zeros(N,1);

sigma = zeros(N,1);
dsigma_dz = zeros(N,1);


for i = 1:N
 
        sigma(i) = (p.nsigma2-p.nsigma1)*tanh((z(i)-p.L)/p.zstar_n) + p.nsigma2; %flare in external gas
        dsigma_dz(i) = (p.nsigma2-p.nsigma1)/p.zstar_n * sech((z(i)-p.L)/p.zstar_n)^2;
        
        if flare
           r(i) = (p.r2-p.r1)*tanh((z(i)-p.L)/p.zstar) + p.r2;
           A(i) = pi*r(i)^2;
           dA_dz(i) = 2*pi*r(i)*(p.r2-p.r1)*((sech((z(i)-p.L)/p.zstar))^2)*(1/p.zstar);   

        else
            r(i) = (p.rL-p.r0)*z(i)/p.L + p.r0;
            A(i) = pi*r(i)^2;
            dA_dz(i) = 2*pi*r(i)*(p.rL-p.r0)/p.L;           
        end
       
         
         if Pbar(i) < (p.nt/p.sw)^2
            dnbar_dP(i) = -p.sw*0.5*Pbar(i)^(-1/2);
            bbar(i) = -dnbar_dP(i); 
            nbar(i) = p.nt-p.sw*sqrt(Pbar(i))+ sigma(i);          %set nbar = neq(P) + sigma(z), where sigma is arbitrary fct of z

     
        else
            nbar(i) = sigma(i);%0;
            dnbar_dP(i) = 0;
            bbar(i) = 0;
            
         end

        rho_liq(i) = p.rho_tilde*(1 + Pbar(i)/p.Kliq); 
        drho_liq_dP(i) = p.rho_tilde/p.Kliq; 
        rhog(i) = Pbar(i)/(p.R*p.T);
        drhog_dP(i) = 1/(p.R*p.T);
        
        rhobar(i) = 1/((nbar(i)/rhog(i)) + (1-nbar(i))/rho_liq(i)); 
        
        drho_dn(i) = -rhobar(i)^(2)*(1/rhog(i) -1/rho_liq(i)); 
          
        drho_dP(i) = -rhobar(i)^(2)*(-nbar(i)*drhog_dP(i)/rhog(i)^2 ...
            + ((nbar(i)-1)*drho_liq_dP(i))/rho_liq(i)^2);
        
        cbar(i) = (drho_dP(i))^(-1/2);
        
         drhobar_dPbar(i) = -rhobar(i)^(2)*(dnbar_dP(i)/rhog(i) - nbar(i)*drhog_dP(i)/rhog(i)^2 + ...
           (-rho_liq(i)*dnbar_dP(i) + (nbar(i)-1)*drho_liq_dP(i))/rho_liq(i)^2);

    
         if Pbar(i) < (p.nt/p.sw)^2
            abar(i) = -(1/rhobar(i)) .* drho_dn(i);  %why is this line here? 
         else
             if dsigma_dz(i) == 0
             abar(i) = 0; %same as above. why is this line here?
             else
             abar(i) = -(1/rhobar(i)) .* drho_dn(i);
             end
         end
         
        ceqbar(i) = sqrt(1/(1/cbar(i)^2 + abar(i)*bbar(i)*rhobar(i)));
        
        phi(i) = nbar(i)/rhog(i)/(nbar(i)/rhog(i) + (1-nbar(i))/rho_liq(i));
        
        dphi_dn(i) = (1/rhog(i))*(nbar(i)/rhog(i) + (1-nbar(i))/rho_liq(i)) + (nbar(i)/rhog(i))*((1/rhog(i)) - 1/rho_liq(i));
        dphi_dp(i) = (-nbar(i)*drhog_dP(i)/rhog(i)^2)*(nbar(i)/rhog(i) + ...
            (1-nbar(i))/rho_liq(i)) + (nbar(i)/rhog(i))*(-nbar(i)*drhog_dP(i)/rhog(i)^2 - (1-nbar(i))*drho_liq_dP(i)/rho_liq(i)^2);
        

    if vbar(i) > 0
            f2 = p.f0;
            f1 = 8*p.eta*vbar(i)/r(i)^2;
            
            if p.phiscale ~= 0
                fbar(i) = (f2-f1)*(0.5 + 0.5*tanh((phi(i)-0.75)./p.phiscale)) + f1;
                Fv(i) = -8*p.eta/r(i)^2*(0.5 + 0.5*tanh((phi(i)-0.75)./p.phiscale)) + 8*p.eta/r(i)^2;
                Fn(i) = (f2-f1)*0.5*(1./p.phiscale)*dphi_dn(i)*(sech((phi(i)-0.75)./p.phiscale))^2;
                Fp(i) = (f2-f1)*0.5*(1./p.phiscale)*dphi_dp(i)*(sech((phi(i)-0.75)./p.phiscale))^2;
            else
                if phi(i) < 0.75
                    fbar(i) = f1;
                    Fv(i)   = 8*p.eta/r(i)^2;
                    Fn(i)   = 0; %Not sure here. Maybe (f1-f2) * dirac(phi(i)-0.75)*dphi_dn(i);
                    Fp(i)   = 0; %Not sure here. Maybe (f1-f2) * dirac(phi(i)-0.75)*dphi_dp(i);
                    
                else
                    fbar(i) = p.f0;
                    Fv(i) = 0;
                end

                
            end
            
            
         else
            fbar(i) = 0; 
            Fv(i)   = 0;
            Fn(i)   = 0;
            Fp(i)   = 0;
    end
   
        %define partial derivatives
        rat = (vbar(i)/ceqbar(i))^2;
        dvbar_dz(i) = (1/(1-rat))*(rat*(rhobar(i)*p.g + fbar(i))/(rhobar(i)*vbar(i)) - (vbar(i)/A(i))*dA_dz(i) + abar(i)*vbar(i)*dsigma_dz(i));
        dPbar_dz(i) = (1/(1-rat))*(-rhobar(i)*p.g - fbar(i) + (rhobar(i)*vbar(i)^2/A(i))*dA_dz(i) - rhobar(i)*abar(i)*vbar(i)^2*dsigma_dz(i));
        dnbar_dz(i) = -bbar(i)*dPbar_dz(i);
        drhobar_dz(i) = drhobar_dPbar(i)*dPbar_dz(i) - rhobar(i)*abar(i)*dsigma_dz(i); 
        Kbar(i)  = rhobar(i)*cbar(i)^2;
        
        b1 = -nbar(i)*drhog_dP(i)/rhog(i)^2 + (nbar(i)-1)*drho_liq_dP(i)/rho_liq(i)^2;
        b2 = ((-dnbar_dP(i)*drhog_dP(i))*rhog(i)^2 + nbar(i)*drhog_dP(i)*2*rhog(i)*drhog_dP(i))/rhog(i).^4;
        b3 = ((dnbar_dP(i)*drho_liq_dP(i))*rho_liq(i)^2 - (nbar(i)-1)*drho_liq_dP(i)*2*rho_liq(i)*drho_liq_dP(i))/rho_liq(i).^4;
        AA = -2*rhobar(i)*drhobar_dPbar(i)*b1 - rhobar(i)^2*(b2 + b3);
       
        dcbar_dz(i) = -0.5*(drho_dP(i))^(-3/2)*AA*dPbar_dz(i);
        dkbar_dz(i) = 2*cbar(i)*dcbar_dz(i)*rhobar(i) + cbar(i)^2*drhobar_dz(i);
        dzbar_dz(i) = cbar(i)*drhobar_dz(i) + rhobar(i)*dcbar_dz(i);

end


fields.rhobar = rhobar;
fields.nbar = nbar;
fields.rho_liq = rho_liq;
fields.rhog = rhog;
fields.dP_drho = dP_drho;
fields.drho_dn = drho_dn;
fields.cbar = cbar;
fields.dcbar_dz = dcbar_dz;
fields.ceqbar = ceqbar;
fields.abar = abar;
fields.phi = phi;
fields.fbar = fbar;
fields.bbar = bbar;
fields.dvbar_dz = dvbar_dz;
fields.drhobar_dz = drhobar_dz;
fields.dPbar_dz = dPbar_dz;
fields.dnbar_dz = dnbar_dz;
fields.Kbar = Kbar;
fields.dkbar_dz = dkbar_dz;
fields.Fv = Fv;
fields.Fn = Fn;
fields.Fp = Fp;
fields.r = r;
fields.A = A;
fields.dA_dz = dA_dz;
fields.zbar = rhobar.*cbar;
fields.dzbar_dz = dzbar_dz;
fields.sigma = sigma;
fields.dsigma_dz = dsigma_dz;

     