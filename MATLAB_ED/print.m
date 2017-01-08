%  Airy functions coeffients [C0 C1 .. ]
c = [2.338 4.088 5.521 6.787 7.944 9.022]; 
for i=7:30
    c(i) = (3*pi/2 * (i-0.25))^(2/3);
end 

global d T

%using the energies we obtain from [] we can estimate the electric field
%at equilibrium (NO BIAS APPLIED)

%Obtain electric field from E1 and E2
El0 = El_fiel_trig_well(E0(V0),c(1));
El1 = El_fiel_trig_well(E1(V0),c(2));

%Relative errors (thereticcaly we should obtain the same value)
Accuracy_EL = abs((El1 - El0)/El0); 

%quantum well depth 
Emax= El0*d;

%using the Electric field estimated below we can compute more bands.
Airy_energies = 1:length(c);
P_occupation = 1:length(c);

for i=1:length(c)
    Airy_energies(i) = Energy_trig_well(El0, c(i))
    P_occupation(i) = fermi(Airy_energies(i), Efermi(V0), T);
end 

function El =  El_fiel_trig_well(En,c)
    global q h m 
    eq = @(El) En - 1./q*c*((q*El*h).^2./(2*m)).^(1/3); %expressed in eV
    El = fzero(eq, 1);
end 

function En =  Energy_trig_well(El,c)
    global q h m 
    eq = @(En) En - 1/q*c*((q*El*h).^2./(2*m)).^(1/3);
    En = fzero(eq, 1);
end 

function P = fermi(E, Ef, T)
   global Kb
   P=1/(1+ exp((E-Ef)/(T*Kb)));
end

