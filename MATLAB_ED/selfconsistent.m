L = 5;
N = 1000;
x = linspace(-L,L,N)';
dx = x(2) - x(1);

trig = 1:200;
U = [100*ones(1,500) trig-200 zeros(1,300)]';
%U = [100*ones(1,100) zeros(1,900)]';

% Three-point finite-difference representation of Laplacian
% using sparse matrices, where you save memory by only
% storing non-zero matrix elements
e = ones(N,1); 
Lap = spdiags([e -2*e e],[-1 0 1],N,N)/dx^2;

for i=1:5
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

%plot(x,V,x,Usc,'--k'); % plot V(x) and rescaled U(x)
% Add legend showing Energy of plotted V(x)

%plot(x,100*V(:,1:10)+E(1:10)',x,U)

e = ones(N,1); 
Lap = spdiags([e -2*e e],[-1 0 1],N,N)/dx^2;
Poisson = Lap;

dielcost = [1*ones(1,500) 20*ones(1,500)]';
%second term
n = 2000*V(:,1).^2 + 1000*V(:,2).^2;
%n = n+ 0.01;
xx = [1:500];
norm = 0.01+0.7*normpdf(xx,500,50);
Nd  = [1*norm*sum(n) zeros(1,500)]';
Tot = -1000*(Nd-40*n)./dielcost;
%solver
[X, R] = linsolve(zeros(N,N)+Poisson, Tot);

DeltaE =  [100000*ones(500,1)' -10000*ones(500,1)']; %zeros(1000,1); 
U = -X+DeltaE';
end
plot(U)