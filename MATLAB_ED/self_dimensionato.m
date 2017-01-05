mo = 9.11e-31;   %mass of electron
h = 1.05e-34;    %Js 
hbar = h/(2*pi);
m = 0.22*mo;     %effective mass cit: http://journals.aps.org/prb/pdf/10.1103/PhysRevB.57.1374
q = 1.60e-19;
Kb = 8.61673e-5;

L = 50e-9; %2L = 100nm
N = 100000;
x = linspace(-L,L,N)';
dx = x(2) - x(1);

trig = -2.7:0.00027/2:-0.0001;
U = [3*ones(1,50000) trig zeros(1,30000)]'; %guess potential

% Three-point finite-difference representation of Laplacian
% using sparse matrices, where you save memory by only
% storing non-zero matrix elements
e = ones(N,1); 
Lap = spdiags([e -2*e e],[-1 0 1],N,N)/dx^2;

for i=1:1
% Total Hamiltonian
% constants for Hamiltonian 
H = -1/2*(hbar^2/m)*Lap + spdiags(U,0,N,N);
%H = -1/2*Lap + spdiags(U,0,N,N);

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

dielcost = [9*ones(1,500) 11*ones(1,500)]';
%second term
n=0;
for ei=1:nmodes
Efermi = 0.403; %V
n = n + m*Kb*300/(pi*hbar^2)*log(1 + exp((Efermi - E(ei))/(Kb*300)))*V(:,ei).^2 ;
end

xx = [1:50000];
norm = 0.00000001+2*sum(n)*normpdf(xx,50000,50);
Nd  = [norm zeros(1,50000)]';
Tot = -q*(Nd-n)./dielcost;
%solver
[X, R] = linsolve(zeros(N,N)+Poisson, Tot);
% 0.195 eV
DeltaE =  -[19.5*ones(500,1)' 0*ones(500,1)']; %zeros(1000,1); 
%U = -X-DeltaE';
end
plot(U)