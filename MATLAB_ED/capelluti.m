clear all

Vmin = -2.95;
Vmax = 3;   %gate voltage

T = 300; %temperature

str = 'temp400.png';


%parameters
global Voff l0 l1 ep D q Vth P2 P1 mo h m Kb hbar d

f = figure()

kk = 0;
% for T=20:20:500
    
Voff= -3;
d = 20e-9;
l0 = 2.12e-12;
l1 = 3.73e-12;   %experimental parameters
ep = 9*8.85e-12; %permittivity of AlGaN
D = 1.001e18;    %density of states
q = 1.60e-19;    %charge of electron in Coulomb
Kb = 8.61673324e-5;
Vth = T*Kb; %bolzamn costant * T
P2 = D*Vth;
P1 = ep/(q*d);
mo = 9.11e-31; %mass of electron
h = 1.05e-34;  %Js 
hbar = h/(2*pi);
m = 0.22*mo;    %effective mass cit: http://journals.aps.org/prb/pdf/10.1103/PhysRevB.57.1374

DeltaV= 100; %resolution - number of points in the plot

Efermi= 1:DeltaV;
Volts = 1:DeltaV;
E1 = 1:DeltaV;
E0 = 1:DeltaV;


    clf
    kk = kk +1;
    i=0;
    for Vg = Vmin:(abs(Vmin-Vmax)/DeltaV):Vmax
        i = i+1;
        Vgo = Vg - Voff;
        eq_1 = @(Ef) -P1*(Vgo - Ef) + P2*(log(exp((Ef - l0*(P1*(Vgo - Ef)).^(2/3))/Vth) + 1) + log(exp((Ef - l1*(P1*(Vgo - Ef)^(2/3)))/(Vth)) + 1 )    );                        
        Efermi(i) = fzero(eq_1,0);
        Volts(i) = Vg;

        ns = @(Vg, Ef) P1*(Vg - Ef);

        %two bands in juction well as define in the paper
        E1(i) = l1*ns(Vgo,Efermi(i))^(2/3);
        E0(i) = l0*ns(Vgo,Efermi(i))^(2/3);
        
        
        %first ragion
        eq_1_approx = @(Ef) -P1*(Vgo - Ef) + P2*(log(exp((Ef - l0*(P1*(Vgo - Ef)).^(2/3))/Vth) + 1)    );                        
        Ef_approx_reg_1(i) = fzero(eq_1_approx,0); 

    end

    

    %figure [1] - Energies vs. Voltages
    fig1_EvsVolts(Volts,Efermi, E0, E1, T ,f)
    hold on
    plot(Volts,Ef_approx_reg_1)
    %ok, il secondo termine è inifluente
    
% end
%find where E0 intersect Ef
V0 = find(Efermi - E0 >= 0);
V0 = V0(1); %index for Vgate = 0 
Region1_2_volts  = Volts(V0);
Region1_2_Efermi = Efermi(V0);

clf

%approssimazione eq 3?
apprx1 = @(Vgo, Ef) exp((Ef - l0*(P1*(Vgo - Ef)).^(2/3))/Vth);
apprx = apprx1(Volts-Voff,Efermi);
rel_err_eq3 = (log(apprx+1)-log(apprx))./log(apprx+1);
plot(Volts,rel_err_eq3) %(1:V0)
hold on 
%axis([-3 -1.5 -40 0])


clf
ns = ns_eq_2(Volts- Voff, Efermi);
Ef_eq_3 = eq_3(Volts - Voff,ns);
plot(Volts,Ef_eq_3)
hold on
plot(Volts,Ef_approx_reg_1,'o')
plot(Volts,Efermi,'-x')
(Ef_eq_3-Ef_approx_reg_1)./Ef_approx_reg_1
plot(Volts,-(Ef_eq_3-Ef_approx_reg_1)./abs(Ef_approx_reg_1))
yL = get(gca,'YLim');
line([Volts(V0) Volts(V0)],yL,'Color','k','LineStyle','-.');

clf

% figure [2] - Delta Energies


% plot(Volts, Efermi - E1)
% hold on 
% plot(Volts, Efermi - E0)
%     
%       
%     xL = get(gca,'XLim');
%     yL = get(gca,'YLim');
%     line([Volts(V0) Volts(V0)],yL,'Color','k','LineStyle','-.');
%     line(xL,[0 0],'Color','k','LineStyle','-.');
% str = ['Diff_Ef_vs_E0_E1_vs_Volts' num2str(T) 'K.png'];
% saveas(f,str);









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

V0 = find(Volts >= 0);
V0 = V0(1); %index for Vgate = 0

%Electric field from E1 and E2
El0 = El_fiel_trig_well(E0(V0),c(1))
El1 = El_fiel_trig_well(E1(V0),c(2))

%relative errors (thereticcaly we should obtain the same value)
Accuracy_EL = abs((El1 - El0)/El0)  

%using the Electric field estimated below we can compute more bands.


Airy_energies = 1:length(c);
P_occupation = 1:length(c);
for i=1:length(c)
    Airy_energies(i) = Energy_trig_well(El0, c(i));
    
    P_occupation(i) = fermi(Airy_energies(i), Efermi(V0), T);
end 

clf

%number of bounded states
Emax= El0*d;
1 + round(sqrt(-(2*m*Emax*q*(d.^2))/(pi.^2*hbar.^2 )))

%plot(Airy_energies(1:length(c)),P_occupation(1:length(c)),'-o')
plot(Airy_energies,'o')
hold on 
plot([0 30],[-Emax -Emax],'k-')

%fare analisi dimensionale, è figa si ottiene N/C


%  forse Efermi AlGaN 3.96eV
%                 Gan 3.39eV
% % 
% % 
function El =  El_fiel_trig_well(En,c)
    global q h m 
    eq = @(El) En - 1/q*c*((q*El*h)^2/(2*m))^(1/3);
    El = fzero(eq, 1);
end 

function En =  Energy_trig_well(El,c)
    global q h m 
    eq = @(En) En - 1/q*c*((q*El*h)^2/(2*m))^(1/3);
    En = fzero(eq, 1);
end 

function P = fermi(E, Ef, T)
   global Kb
   P=1/(1+ exp((E-Ef)/(T*Kb)));
end

function ns = ns_eq_2(Volt, Ef)
   global q ep d Voff 
   ns = ep/(d*q)*(Volt - Ef);
end

function Ef = eq_3(Vg, ns)
 global l0 Vth D
 Ef = l0*(ns).^(2/3) + Vth*log(ns./(D*Vth));
end
