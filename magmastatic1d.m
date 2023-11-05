function [Out] = magmastatic1d(Param,N,F,x,flag)
%calculate steady state solution for flow in conduit for reduced equations
%evaluate at the grid points defined by N
%F contains a flag to determine wheter outputs are "frozen" (for testing)

options = odeset('RelTol',1e-13,'AbsTol',1e0*eps,'Events',@(t,y) event(t,y,Param));%,'Refine',100);%,'Events',@(t,y) events(t,y,Param));%'OutputFcn', @save_stream,,'OutputFcn',@odeplot,'OutputSel',[2]);

Pch=0e6; %Chamber overpressure in Pa

Param.rescale=0;

%initial value of Pressure at surface, 
if Param.rescale
InitP= 1e5/1e6;%Param.Rhol0*9.8*Param.L + Pch;
else
 InitP= 1e5;
end

%solve the conduit flow ODE    
sol=ode45(@(t,y) RHSMagmaStaticSolver(t,y,Param),[0 Param.L],[InitP],options);  
%keyboard
if flag
Out.Ex_above = Param.L-max(sol.xe(sol.ye<Param.GasExsolution));
Out.Ex_below = Param.L-min(sol.xe(sol.ye>Param.GasExsolution));
%Out.PEx = min(sol.ye);
return
end
x=Param.L-x;
%evaluate solution and derivatives at high order
[P,dPdz]=deval(x,sol); 

%keyboard
LEN=length(P);

 for ii=1:LEN    

         if(P(ii)<=Param.GasExsolution)
             if Param.GasLookup
     N(ii) = interp1(Param.GasPress,Param.GasTot,P(ii));
             else
     N(ii) = (Param.nt-Param.sw .*P(ii).^(1/2)); %equilibrium gas exsolution for water 
             end
     R(ii) = (N(ii) .*Param.R .*Param.T./P(ii) + (1-N(ii))./(Param.rho_tilde.*(1+ P(ii)./Param.Kliq))).^(-1);
     
     a0(ii) = -1./R(ii) .* (P(ii).*(Param.Kliq+P(ii)).*Param.rho_tilde.*(Param.Kliq.*P(ii)+(-1).*(Param.Kliq+P(ii)).*Param.R.*Param.rho_tilde.*Param.T).*( ...
            Param.Kliq.*((-1)+N(ii)).*P(ii)+(-1).*N(ii).*(Param.Kliq+P(ii)).*Param.R.*Param.rho_tilde.*Param.T).^(-2));
     b0(ii) = Param.sw./(2.*sqrt(P(ii)));    
        
         else
     N(ii)=0;
     R(ii) = (Param.rho_tilde.*(1+ P(ii)./Param.Kliq));
     
     a0(ii) = 0;
     b0(ii) = 0;
         end

 end

C0 = (Param.rho_tilde.^(-1).*(Param.Kliq.*((-1)+N).*P+(-1).*N.*(Param.Kliq+P).*Param.R.*Param.rho_tilde.* ...
    Param.T).^2.*((-1).*Param.Kliq.*((-1)+N).*P.^2+N.*(Param.Kliq+P).^2.*Param.R.*Param.rho_tilde.* ...
    Param.T).^(-1)).^(1/2);  

dNdz=R.*Param.g.*b0;
    
dRdz=-R.*Param.g./C0.^2 - R.^2.*Param.g.*a0.*b0; 
    
Ceq = 1./sqrt(1./C0.^2 + a0/Param.g .*dNdz);   

Nfreq= Param.g .*sqrt(1./Ceq.^2 - 1./C0.^2);

%smooth transition
%[vt tt]=min(abs(Param.GasExsolution - P));
%ex = 1-tanh((x/x(tt)).^20); % smoothed transition
%N=ex.*N;
%a0=ex.*a0;
%b0=ex.*b0;

% if(F.Fixed)
%  CD=.8; %depth at which coefficients frozen, <1 is below exsolution
%         [MP IPx]=min(abs(P-Param.GasExsolution/CD));    
%      %frozen coefficients, flip outputs so that z is positive up
%     Out.P = fliplr(P(IPx))*ones(1,length(x))';
%     Out.R = fliplr(R(IPx))*ones(1,length(x))'; 
%     Out.N = fliplr(N(IPx))*ones(1,length(x))'; 
%     Out.dRdz = dRdz(IPx)*ones(size(Out.P))'; 
%     Out.dNdz = dNdz(IPx)*ones(size(Out.P))'; 
%     Out.C0 = fliplr(C0(IPx))*ones(1,length(x))'; 
%     Out.a0 = fliplr(a0(IPx))*ones(1,length(x))'; 
%     Out.b0 = fliplr(b0(IPx))*ones(1,length(x))'; 
%     Out.Ceq = fliplr(Ceq(IPx))*ones(1,length(x))'; 
%     Out.Nfreq = fliplr(Nfreq(IPx))*ones(1,length(x))'; 
%     Out.b0 = fliplr(Nfreq(IPx))*ones(1,length(x))'; 
%     Out.x = x';
% else
    %flip outputs 
% Out.P=fliplr(P)';
% Out.R=fliplr(R)';
% Out.N=fliplr(N)';
% Out.dRdz=flipud(dRdz');
% Out.dNdz=flipud(dNdz');
% Out.C0=fliplr(C0)';
% Out.a0=fliplr(a0)';
% Out.Ceq = fliplr(Ceq)';
% Out.Nfreq = fliplr(Nfreq)';
% Out.b0 = fliplr(b0)';
% Out.x=x';
Out.P=P';
Out.rhobar=R';
Out.nbar=N';
Out.drhobar_dz=dRdz';
Out.dPbar_dz = -Param.g*R';
Out.dNdz=dNdz';
Out.cbar=C0';
Out.abar=a0';
Out.ceqbar = Ceq';
Out.Nfreq = Nfreq';
Out.bbar = b0';
Out.x=Param.L-x';
Out.Kbar = R' .* transpose(C0).^2;

Out.vbar=0*C0'; %zero vector
Out.dvbar_dz=0*C0';
Out.Fv=0*C0'; %zero vector
Out.Fn=0*C0'; %zero vector
Out.Fp=0*C0'; %zero vector
Out.dA_dz=0*C0'; %zero vector
Out.A = pi*Param.r0*ones(size(transpose(C0)));

Out.zbar = R' .* transpose(C0);


rho_liq = Param.rho_tilde*(1 + P./Param.Kliq); 
drhol_dP = Param.rho_tilde./Param.Kliq; 
rhog = P./(Param.R*Param.T);
drhog_dP = 1/(Param.R*Param.T);

drho_dn = -R.^2.*(1./rhog -1./rho_liq); 
          
drho_dP = R.^2.*(N.*drhog_dP./rhog.^2 ...
            + ((1-N).*drhol_dP)./rho_liq.^2);

d2rho_dP2 = -2.*R.^3 .*(N./rhog.^2.*drhog_dP + (1-N)./rho_liq.^2.*(1-N).*drhol_dP).^2+...
    R.^2 .*(-2.*N./rhog.^3 .*drhog_dP.^2 - 2.*(1-N)./rho_liq.^3 .*drhol_dP.^2);

dCdP = -.5 .*drho_dP.^(-3/2) .* d2rho_dP2;

dCdN = (.5 .*(N.*rho_liq+(1-N).*rhog).*(rhog.^2*drhol_dP-rho_liq.^2 .*drhog_dP) - ...
    (N.*rho_liq.^2 .*drhog_dP + (1-N).*rhog.^2 .*drhol_dP).*(rhog-rho_liq)) .* ...
    ((N.*rho_liq.^2 .*drhog_dP + (1-N).*rhog.^3 .*drhol_dP)./(N.*rho_liq+(1-N).*rhog).^2).^(-3/2) .* ...
    (N.*rho_liq+(1-N).*rhog).^(-3);

% drhobar_dPbar = -R.^(2).*(b0./rhog - N.*drhog_dP./rhog.^2 + ...
%            (-rho_liq.*b0 + (N-1).*drho_liq_dP)./rho_liq.^2);
% 
% b1 = -N.*drhog_dP./rhog.^2 + (N-1).*drho_liq_dP./rho_liq.^2;
% b2 = ((-b0.*drhog_dP).*rhog.^2 + N.*drhog_dP.*2.*rhog.*drhog_dP)./rhog.^4;
% b3 = ((b0.*drho_liq_dP).*rho_liq.^2 - (N-1)*drho_liq_dP.*2.*rho_liq.*drho_liq_dP)./rho_liq.^4;
% AA = -2.*R.*drhobar_dPbar.*b1 - R.^2.*(b2 + b3);
% 
% dcbar_dz = -0.5.*(drho_dP).^(-3/2).*AA.*Out.dPbar_dz';
% dkbar_dz = 2.*C0.*dcbar_dz.*R + C0.^2.*dRdz;

% Out.dzbar_dz = C0'.*dRdz' + R'.*dcbar_dz';
% Out.dcbar_dz=dcbar_dz';
% Out.dkbar_dz = dkbar_dz';


d2Zinv_dz = (Param.g ./(2.*C0) .*(1./R .*(drho_dP - Param.sw./(2*sqrt(P)).*drho_dn) + ...
    1./C0.*(dCdP - Param.sw./(2.*sqrt(P)).*dCdN))); %term 1 from Jasper

d2C2inv_dz = R.*Param.g./C0 .*(dCdP - Param.sw./(2*sqrt(P)).*dCdN); %term 2

dsqrt2K_dz = Param.g./(sqrt(8*R).*C0.^2) .* (C0.*(drho_dP - Param.sw./(2*sqrt(P)).*drho_dn) +...
    2*R.*(dCdP - Param.sw./(2*sqrt(P)).*dCdN)); %term 3 from Jasper

Out.d2Zinv_dz = d2Zinv_dz';
Out.d2C2inv_dz = d2C2inv_dz';
Out.dsqrt2K_dz = dsqrt2K_dz';

Out.dnbar_dz = Out.dNdz;
%end

%keyboard
end

function RHS = RHSMagmaStaticSolver(t,y,Param)
%right hand side of ODE in pressure and velocity for steady conduit flow

P=y;%(1);
         if(P(1)<=Param.GasExsolution)
             if Param.GasLookup
                 %if(P>1e6)
                 %    keyboard
                 %end
     N = interp1(Param.GasPress,Param.GasTot,P);    
             else
     N = (Param.nt-Param.sw .*P.^(1/2)); %equilibrium gas exsolution  
             end
         else
     N=0;      
         end
%mixture density 
    Rho = (N .*Param.R .*Param.T./P + (1-N)./(Param.rho_tilde.*(1+ P./Param.Kliq))).^(-1);   

RHS = (Rho.*Param.g );

end

function [value,isterminal,direction] = event(t,y,Param)
%keyboard
if Param.GasLookup
    if Param.rescale
threshold  =  1e6/1e6;%3e2;%2e2  
    else
threshold  =  1e6;%3e2;%2e2
    end
else
     if Param.rescale
threshold  =  3e2/1e6;%3e2;%2e2  
    else
threshold  =  3e2;%3e2;%2e2
    end
end

value      = abs(y - Param.GasExsolution)-threshold;
% if abs(y - Param.GasExsolution)>threshold
%     keyboard
% end
isterminal = 0;
direction  = 0;
end