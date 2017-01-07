
L = 50; %2L = 100nm
N = 1000;
x = linspace(-L,L,N)';
dx = x(2) - x(1);

Eltrig = -2.256463756872821e-3*51.42208619083232; %eV/51 = ATOMIC UNIT : HA/BOHR
DeltaPot = Eltrig*dx;
trig = 1:500; %-0.2:0.2/900:-0.00002;

U = [50*ones(1,500) -DeltaPot*trig]';


% trig = -2.7:0.027/2:-0.01;
% U = [3*ones(1,500) trig zeros(1,300)]'; %guess potential

% Three-point finite-difference representation of Laplacian
% using sparse matrices, where you save memory by only
% storing non-zero matrix elements
e = ones(N,1); 
Lap = spdiags([e -2*e e],[-1 0 1],N,N)/dx^2;

for i=1:1
    %% 
% Total Hamiltonian
% constants for Hamiltonian 
hbar = 1; m = 0.24; % constants for Hamiltonian 
H = -1/2*(hbar^2/m)*Lap + spdiags(U,0,N,N);
%H = -1/2*Lap + spdiags(U,0,N,N);

% Find lowest nmodes eigenvectors and eigenvalues of sparse matrix 
nmodes = 4; options.disp = 0;
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

dielcost = [9*ones(1,500)*8e-12 11*ones(1,500)*8e-12]';
%second term
n=0;
for ei=1:nmodes
Efermi = 0.403; %V
mo = 9.11e-31;   %mass of electron
hbar = 1.05e-34;    %Js 
h = hbar*(2*pi);
m = 0.24*mo;
Kb1 = 1.3806485e-23;
Kb = 8.61673e-5;

n = n + m*Kb1*300/(pi*hbar^2)*log(1 + exp((Efermi - E(ei))/(Kb*300)))*V(:,ei).^2 ;
end

xx = [1:500];
norm = 0.00001+2*sum(n)*normpdf(xx,500,50);
Nd  = [norm zeros(1,500)]';
Tot = -q*(Nd-n)./dielcost;
%solver
[X, R] = linsolve(zeros(N,N)+Poisson, Tot);
% 0.195 eV
X = X; %max(abs(q*X))*10
DeltaE =  -[1.95*ones(500,1)' 0*ones(500,1)']; %zeros(1000,1); 
%U = -q*X-DeltaE';

end
plot(U)
plot(x,50*V(:,1:1)+E(1:1)',x,U)