function [El0, El1, Emax, Airy_energies, P_occupation] = airy_fuctions(Efermi, E0, E1, V0, f)

global d T m q hbar 


%Energies and Elecctric field in an idea triangular well. 
%Airy functions
%En = 1/q*c(1)*((q*5e6*h)^2/(2*m))^(1/3) %expressed in eV


%  Airy functions coeffients [C0 C1 .. ]
c = [2.338 4.088 5.521 6.787 7.944 9.022]; 
for i=7:30
    c(i) = (3*pi/2 * (i-0.25))^(2/3);
end 

%using the energies we obtain from [] we can estimate the electric field
%at equilibrium (NO BIAS APPLIED)

%Obtain electric field from E1 and E2
El0 = El_fiel_trig_well(E0(V0),c(1));
El1 = El_fiel_trig_well(E1(V0),c(2));

%Relative errors (thereticcaly we should obtain the same value)
Accuracy_EL = abs((El1 - El0)/El0)  

%using the Electric field estimated below we can compute more bands.
Airy_energies = 1:length(c);
P_occupation = 1:length(c);
for i=1:length(c)
    Airy_energies(i) = Energy_trig_well(El0, c(i))
    P_occupation(i) = fermi(Airy_energies(i), Efermi(V0), T);
end 


Emax= El0*d;

%number of bounded states
%1 + round(sqrt(-(2*m*Emax*q*(d.^2))/(pi.^2*hbar.^2 )))

clf
plot(Airy_energies,'o')
hold on 
plot([0 30],[-Emax -Emax],'k-')

    title('Bound states in triangular well','FontSize',15)
    xlabel('Bound States E_n','FontSize',15) 
    ylabel('Energy [eV]','FontSize',15)
    h_legend=legend('E_n Bound States in the well','Maximum Allowed Energy');
    set(h_legend,'Location','southeast','FontSize',14);
    saveas(f,'Trig_well','svg');





%fare analisi dimensionale, è interessante si ottiene N/C

end 

%  forse Efermi AlGaN 3.96eV
%                 Gan 3.39eV

function El =  El_fiel_trig_well(En,c)
    global q h m 
    eq = @(El) En - 1./q*c*((q*El*h).^2./(2*m)).^(1/3);
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

