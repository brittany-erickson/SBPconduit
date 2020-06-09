function [x,dx,D,Hinv00,A,H] = grid_operators(N,L,order)
% Special distribution of grid points for the 8th order operator
if(order == 8)
    x1 = 3.8118550247622e-01; x2 = 1.1899550868338e+00;
    x3 = 2.2476300175641e+00; x4 = 3.3192851303204e+00;
    
    h = L/(N-8+2*x4); dx = h;

    x = [0,h*x1,h*x2,h*x3,h*x4,...
        linspace(h*x4+h,L-(h*x4+h),N-9),...
        L - h*x4,L - h*x3,L - h*x2,L - h*x1,L]';
else 
    dx = L/N; % length of domain, grid spacing
    x = (0:dx:L)'; % coordinates (at N+1 grid points)
end

% SBP differentiation matrix
% D = first derivative
% Hinv = corner entry of diagonal norm (for SAT penalties)
% A = artificial dissipation operator
% order = interior order of accuracy

[D,Hinv00,A,H] = SBPoperators(N,dx,order); % sparse arrays

end

