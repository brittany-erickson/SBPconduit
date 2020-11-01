% script to process eigenmodes from conduit code
function eigenmodes(A,p,N) 

%% solve the eigenproblem
%unpack the physical variables
%    IV=1:N+1; IRho=N+2:2*(N+1); IP=2*(N+1)+1:3*(N+1); IN=3*(N+1)+1:4*(N+1); 
    IV=1:N; IRho=N+1:2*N; IP=2*N+1:3*N; IN=3*N+1:4*N; 

    [Q,Sm] = eig(A); 
    S = diag(Sm); % solve eigenproblem
    
    %focus on lower frequency modes
   [s1,k] = sort(abs(real(S)),1,'ascend'); 
%    k=find(abs(imag(S))<25);%find(imag(S)~=0&abs(imag(S))<50);
%    Sv=S(k);
%    Qv = Q(:,k);
%    [s1,k1] = sort(abs(real(Sv)),1,'ascend');

for ii=1:length(S)
  %      Hp = ctranspose(Q(IP,ii))*H*Q(IP,ii);Hv = ctranspose(Q(IV,ii))*H*Q(IV,ii);
        pp(:,ii) = Q(IP,ii);%/Hp;
        vv(:,ii) = Q(IV,ii);
        nn(:,ii) = Q(IN,ii);
        rr(:,ii) = Q(IRho,ii);
end



%% Filter out modes with power near the Nyquist

% spatial sampling frequency (1/m)
fs=2*(N+1)/p.L;
% length of signal
Nt=2*(N+1);
% number of FFT points
nfft=2^nextpow2(Nt);
%spatial frequency plotting vector
fr=fs/2*[-1:2/nfft:1-2/nfft];
%define thresholds abovewhich to cut off
thr1=15;
HiFrq = find(fr>1/(2/fs +thr1));
thr2 = -4; %threshold in power
% Analysis
Good=ones(size(S));
% analyze spectrum
for ii=1:length(S)
    FS(:,ii)=fftshift(fft(Q(:,ii),nfft)/Nt);
    % power spectrum (actually sqrt of)
    PS(:,ii)=(log10(abs(FS(:,ii))));
    if(any(PS(HiFrq,ii)>thr2))    
        Good(ii)=0;
    end
end
%indexes of spatially well resolved modes
GdId = find(Good==1);

IQ = Q\eye(size(Q));
obs=zeros(1,length(S));
obs(2*(N+1))=1; %pick out one mode to study

HIQ = obs*IQ;

numr=1;
numi=length(S(GdId));%200;
ks=logspace(-1.5,1,numi/2);
ims=imag(S(GdId));%[-ks ks];

for jj=1:numi
    sv = 1i*ims(jj); %test values
    
    Ft= exp(sv*(sv*M.pT.T^2-M.pT.t))*(1-erfz(.5*(1/M.pT.T)*(2*sv*M.pT.T^2-M.pT.t)));
    Tsval = HIQ(GdId)*diag(1./(sv-S(GdId)))*R(GdId)*Ft;
    
    sl((ii-1)*numi+jj)=sv;
    Ts((ii-1)*numi+jj) = Tsval;
end

%% Energy integrals

 KE=zeros(size(v1));
  AE=zeros(size(v1));

 for i=1:length(S); %do this in the SBP H norm
    KE(:,i) = 1/2.*(M.rho.*v1(:,i)).*v1(:,i);%(ctranspose(v(:,i))*Hnorm*v(:,i));               
    %Kinetic(i) = trapz(x,KE(:,i));%(ww*KE(:,i));  
    Kinetic(i) = real(1/2.*ctranspose(M.rho.*v1(:,i))*H*v1(:,i));
    
    AE(:,i) = ((1./(2.*M.rho.*M.c.^2).*p1(:,i)).*p1(:,i));      
    %Acoustic(i) = trapz(x,AE(:,i));%(ww*AE(:,i));
    Acoustic(i) =  real(ctranspose(1./(2.*M.rho.*M.c.^2).*p1(:,i))*H*p1(:,i)); 
 end 
%% plotting
figure(1)
plot(S,'o')
xlabel('real(s)')
ylabel('imag(s)')

% figure(2)
%scatter(2*pi./abs(imag(S(k))),abs(R),40,abs(real(S(k))),'filled');
%set(gca,'xscale','log');set(gca,'yscale','log');colorbar; 
%xlabel('Mode period (sec)'); ylabel('abs(Q*f)'); title('color by damping rate')


figure(2)
IS=(imag(S));
scatter(2*pi./IS(GdId),abs(R(GdId)),40,log10(abs(real(S(GdId)))),'filled'); xlim([0 65]); box on;colorbar;set(gca,'xscale','log');

%scatter(2*pi./IS(GdId),abs(R(GdId)),40,abs(real(S(GdId))),'filled'); xlim([0 65])

figure(3)
scatter(2*pi./abs(imag(sl)),abs(Ts)/max(abs(Ts)),40,sign(imag(sl)),'filled'); xlim([0 65]); box on;set(gca,'xscale','log');colorbar;

end