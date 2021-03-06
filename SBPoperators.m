function [D,Hinv00,A,H] = SBPoperators(N,h,order)

  if nargin==1, h=1; order=2; end
  if nargin==2, order=2; end
  
  m = N+1;
  
  % build differentiation matrix D and artificial dissipation A
  switch order
    
   case 2
    
    Pinv = [2 ones(1,N-1) 2];
    Q = 0.5*sparse(diag([-1;zeros(N-1,1);1])...
                   +diag(ones(N,1),1)-diag(ones(N,1),-1));
    
    DD = buildAD(m,1);
    B = spdiags([1; ones(N,1)],0,m,m)/2;
    
   case 4
    
    if N<8,disp('N>=8 for 4th order'),return,end
    
    Pinv = [48/17 48/59 48/43 48/49 ...
            ones(1,N-7) 48/49 48/43 48/59 48/17];
    diagQ = ones(N+1,1)*[1/12 -2/3 0 2/3 -1/12];
    Q = spdiags(diagQ,(-2:2),N+1,N+1);
    Q(1:4,1:4)=[-1/2 59/96 -1/12 -1/32;-59/96 0 59/96 0;...
                1/12 -59/96 0 59/96; 1/32 0 -59/96 0];
    Q(N-2:N+1,N-2:N+1) = -rot90(Q(1:4,1:4),2);
    
    DD = buildAD(m,2);
    
    B = spdiags(ones(N+1,1),0,m,m)/12;
    
   case 6
    
    if N<12,disp('N>=12 for 6th order'),return,end
    
    Pinv = [43200/13649 8640/12013 4320/2711 4320/5359 8640/7877 43200/43801 ...
            ones(1,N+1-12) ...
            43200/43801 8640/7877 4320/5359 4320/2711 8640/12013 43200/13649];
    diagQ = ones(N+1,1)*[-1 9 -45 0 45 -9 1]/60;
    Q = spdiags(diagQ,(-3:3),N+1,N+1);   
    
    Q(1,1) = -1/2;
    Q(1,2) =  104009/172800;
    Q(1,3) =  30443/259200;
    Q(1,4) = -33311/86400;
    Q(1,5) =  16863/86400;
    Q(1,6) = -15025/518400;
    
    Q(2,1) = -104009/172800;
    Q(2,2) =  0;
    Q(2,3) = -311/51840;
    Q(2,4) =  20229/17280;
    Q(2,5) = -24337/34560;
    Q(2,6) =  36661/259200;
    
    Q(3,1) = -30443/259200;
    Q(3,2) =  311/51840;
    Q(3,3) =  0;
    Q(3,4) = -11155/25920;
    Q(3,5) =  41287/51840;
    Q(3,6) = -21999/86400;
    
    Q(4,1) =  33311/86400;
    Q(4,2) = -20229/17280;
    Q(4,3) =  11155/25920;
    Q(4,4) =  0;
    Q(4,5) =  4147/17280;
    Q(4,6) =  25427/259200;

    Q(5,1) = -16863/86400;
    Q(5,2) =  24337/34560;
    Q(5,3) = -41287/51840;
    Q(5,4) = -4147/17280;
    Q(5,5) =  0;
    Q(5,6) =  342523/518400;

    Q(6,1) =  15025/518400;
    Q(6,2) = -36661/259200;
    Q(6,3) =  21999/86400;
    Q(6,4) = -25427/259200;
    Q(6,5) = -342523/518400;
    Q(6,6) =  0;

    Q(N-4:N+1,N-4:N+1) = -rot90(Q(1:6,1:6),2);

    DD = buildAD(m,3);
    B = spdiags([1; ones(N,1)],0,m,m)/60;
    
   case 8
    % Optimal 8th order, nonequidistant grid    
       
    mmm = 8; m = N+1;
    Hinv00 = 1/h*1/1.0758368078310e-01; 
    
    P=spdiags(ones(m,1),0,m,m);
    P_U =[0.10758368078309999999e0 0 0 0 0 0 0 0; 0 0.61909685107891000003e0 0 0 0 0 0 0; 0 0 0.96971176519117000000e0 0 0 0 0 0; 0 0 0 0.11023441350947000001e1 0 0 0 0; 0 0 0 0 0.10244688965833000000e1 0 0 0; 0 0 0 0 0 0.99533550116830999995e0 0 0; 0 0 0 0 0 0 0.10008236941028000000e1 0; 0 0 0 0 0 0 0 0.99992060631812000004e0;];
    P(1:mmm,1:mmm)=P_U;
    P(m-mmm+1:m,m-mmm+1:m)=rot90( P_U(1:mmm,1:mmm) ,2 );
    
    DD = buildAD(m,40);
    B = spdiags([1; ones(N,1)],0,m,m)/280;

    D_i1 = 0.800000000000000; 
    D_i2 = -0.200000000000000;  
    D_i3 = 0.038095238095238;  
    D_i4 = -0.003571428571429; 

    diagD1 = ones(N+1,1)*[-D_i4 , -D_i3, -D_i2, -D_i1, 0, D_i1, D_i2, D_i3, D_i4];
    D1 = spdiags(diagD1,(-4:4),N+1,N+1);   

    D1(1,1) = -4.6475450213316e+00;
    D1(1,2) =  6.2541786625636e+00;
    D1(1,3) = -2.4139100510440e+00;
    D1(1,4) =  1.2566395095719e+00;
    D1(1,5) = -6.4766769665992e-01;
    D1(1,6) =  2.4570663393331e-01;
    D1(1,7) = -5.2045357676976e-02;
    D1(1,8) =  4.6433206436931e-03;

    D1(2,1) = -1.0868211647678e+00;
    D1(2,2) =  0.0000000000000e+00;
    D1(2,3) =  1.5195364183857e+00;
    D1(2,4) = -6.5436679827908e-01;
    D1(2,5) =  3.1286206957080e-01;
    D1(2,6) = -1.1086808101812e-01;
    D1(2,7) =  2.1234895991755e-02;
    D1(2,8) = -1.5773398832359e-03;

    D1(3,1) =  2.6780878369504e-01;
    D1(3,2) = -9.7012354133588e-01;
    D1(3,3) =  0.0000000000000e+00;
    D1(3,4) =  9.7262296639767e-01;
    D1(3,5) = -3.6843978324227e-01;
    D1(3,6) =  1.1618593544435e-01;
    D1(3,7) = -1.8907619882972e-02;
    D1(3,8) =  8.5325892405388e-04;

    D1(4,1) = -1.2264219453174e-01;
    D1(4,2) =  3.6750449461988e-01;
    D1(4,3) = -8.5559845023345e-01;
    D1(4,4) =  0.0000000000000e+00;
    D1(4,5) =  7.9552641570539e-01;
    D1(4,6) = -2.2405034810100e-01;
    D1(4,7) =  4.2900978549802e-02;
    D1(4,8) = -3.6408960088886e-03;

    D1(5,1) =  6.8014241294556e-02;
    D1(5,2) = -1.8906569319898e-01;
    D1(5,3) =  3.4874693977150e-01;
    D1(5,4) = -8.5599853894095e-01;
    D1(5,5) =  0.0000000000000e+00;
    D1(5,6) =  7.9186344382307e-01;
    D1(5,7) = -1.9783080393206e-01;
    D1(5,8) =  3.7756538075870e-02;
    D1(5,9) = -3.4861268930073e-03;

    D1(6,1) = -2.6557903380662e-02;
    D1(6,2) =  6.8959742481719e-02;
    D1(6,3) = -1.1319486586973e-01;
    D1(6,4) =  2.4813802673085e-01;
    D1(6,5) = -8.1504122739101e-01;
    D1(6,6) =  0.0000000000000e+00;
    D1(6,7) =  8.0483962089931e-01;
    D1(6,8) = -2.0182899393713e-01;
    D1(6,9) =  3.8273766032180e-02;
    D1(6,10) = -3.5881655655168e-03;

    D1(7,1) =  5.5946228886809e-03;
    D1(7,2) = -1.3135637494343e-02;
    D1(7,3) =  1.8319851498636e-02;
    D1(7,4) = -4.7252720307139e-02;
    D1(7,5) =  2.0250470348441e-01;
    D1(7,6) = -8.0042614113573e-01;
    D1(7,7) =  0.0000000000000e+00;
    D1(7,8) =  7.9973532195798e-01;
    D1(7,9) = -1.9983539676215e-01;
    D1(7,10) =  3.8063885097553e-02;
    D1(7,11) = -3.5684892278956e-03;

    D1(8,1) = -4.9958518981227e-04;
    D1(8,2) =  9.7660369095530e-04;
    D1(8,3) = -8.2748091416587e-04;
    D1(8,4) =  4.0138390353475e-03;
    D1(8,5) = -3.8683470124513e-02;
    D1(8,6) =  2.0090351330033e-01;
    D1(8,7) = -8.0045761050338e-01;
    D1(8,8) =  0.0000000000000e+00;
    D1(8,9) =  8.0006351998859e-01;
    D1(8,10) = -2.0001587999715e-01;
    D1(8,11) =  3.8098262856599e-02;
    D1(8,12) = -3.5717121428062e-03;

    D1(N-6:N+1,N-10:N+1) = -rot90(D1(1:8,1:12),2);
    D = 1/h*D1;

  end
  
  if(order == 8)
      H = h*P;
      A = -inv(P)*transpose(DD)*B*DD;
  else
  
      D = spdiags(Pinv',0,N+1,N+1)*Q/h;
      Hinv00 = Pinv(1,1)/h;
      

      A = -spdiags(Pinv',0,N+1,N+1)*transpose(DD)*B*DD;
      H = h*spdiags(1./Pinv',0,m,m);

  end
  
end

  
