% % 
% % L = 2*pi;
% % N = 100;
% % x = linspace(0,L,N)'; 
% % dx = x(2) - x(1);
% % % Interval Length
% % % No. of coordinate points % Coordinate vector
% % % Coordinate step
% % % Two-point finite-difference representation of Derivative 
% % D=(diag(ones((N-1),1),1)-diag(ones((N-1),1),-1))/(2*dx);
% % % Next modify D so that it is consistent with f(0) = f(L) = 0 D(1,1) = 0; D(1,2) = 0; D(2,1) = 0; % So that f(0) = 0 D(N,N-1) = 0; D(N-1,N) = 0; D(N,N) = 0; % So that f(L) = 0
% % % Three-point finite-difference representation of Laplacian 
% % Lap = (-2*diag(ones(N,1),0) + diag(ones((N-1),1),1)+ diag(ones((N-1),1),-1))/(dx^2);
% % % Next modify Lap so that it is consistent with f(0) = f(L) = 0 Lap(1,1) = 0; Lap(1,2) = 0; Lap(2,1) = 0; % So that f(0) = 0 Lap(N,N-1) = 0; Lap(N-1,N) = 0; Lap(N,N) = 0;% So that f(L) = 0
% % % To verify that D*f corresponds to taking the derivative of f % and Lap*f corresponds to taking a second derviative of f,
% % % define f = sin(x) or choose your own f
% % f = sin(x);
% % % And try the following:
% % Df = D*f; Lapf = Lap*f;
% % plot(x,f,'b',x,Df,'r', x,Lapf,'g');
% % axis([0 5 -1.1 1.1]); % Optimized axis parameters
% % % To display the matrix D on screen, simply type D and return ... D % Displays the matrix D in the workspace
% % Lap % Displays the matrix Lap
% % 
% % % Total Hamiltonian where hbar=1 and m=1
% % hbar = 1; m = 1; 
% % H = -(1/2)*(hbar^2/m)*Lap;
% % % Solve for eigenvector matrix V and eigenvalue matrix E of H
% % [V,E] = eig(H);
% % % Plot lowest 3 eigenfunctions 
% % plot(x,V(:,3),'r',x,V(:,4),'b',x,V(:,5),'k'); shg;
% % E % display eigenvalue matrix
% % diag(E) % display a vector containing the eigenvaluesdisplay a vector containing the eigenvalues


L = 5;
N = 1000;
x = linspace(-L,L,N)';
dx = x(2) - x(1);
% POTENTIAL, choose one or make your own
U = 1/2*10*x.^(2); % quadratic harmonic oscillator potential %
%U = 1/2*10*x.^(4); % quartic potential
% Finite square well of width 2w and depth given %
w = L/50;
trig = 1:400;
U = [100*ones(1,100) trig-400 zeros(1,500)]';
% Two finite square wells of width 2w and distance 2a apart %
w = L/50; 
a=3*w;
%U = -200*(heaviside(x+w-a) - heaviside(x-w-a) + heaviside(x+w+a) - heaviside(x-w+a));
% Three-point finite-difference representation of Laplacian
% using sparse matrices, where you save memory by only
% storing non-zero matrix elements
e = ones(N,1); Lap = spdiags([e -2*e e],[-1 0 1],N,N)/dx^2;
% Total Hamiltonian
hbar = 1; m = 1; % constants for Hamiltonian 
H = -1/2*(hbar^2/m)*Lap + spdiags(U,0,N,N);
% Find lowest nmodes eigenvectors and eigenvalues of sparse matrix 
nmodes = 10; options.disp = 0;
[V,E] = eigs(H,nmodes,'sa',options); % find eigs
[E,ind] = sort(diag(E));% convert E to vector and sort low to high 
V = V(:,ind); % rearrange corresponding eigenvectors
% Generate plot of lowest energy eigenvectors V(x) and U(x)
Usc = U*max(abs(V(:)))/max(abs(U)); % rescale U for plotting 
plot(x,V,x,Usc,'--k'); % plot V(x) and rescaled U(x)
% Add legend showing Energy of plotted V(x)
lgnd_str = [repmat('E = ',nmodes,1),num2str(E)];
legend(lgnd_str) % place lengend string on plot
shg

plot(x,100*V(:,1:10)+E(1:10)',x,U)
hold on 
surface(A)