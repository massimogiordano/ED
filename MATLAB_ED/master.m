%clear all

Vmin = -2.95;
Vmax = 3;   %gate voltage


%parameters
global Voff l0 l1 ep D q mo h m Kb hbar d DeltaV T Reg1_2

T = 300; %temperature
% for T=100:20:500
f = figure()
    
Voff= -3;

d = 20e-9;
l0 = 2.12e-12;
l1 = 3.73e-12;   %experimental parameters
ep = 9*8.85e-12; %permittivity of AlGaN
D = 1.001e18;    %density of states
q = 1.60e-19;    %charge of electron in Coulomb
Kb = 8.61673e-5;
mo = 9.11e-31;   %mass of electron
hbar = 1.05e-34;    %Js 
h = hbar*(2*pi);
m = 0.24*mo;     %effective mass cit: http://journals.aps.org/prb/pdf/10.1103/PhysRevB.57.1374

DeltaV= 100;     %resolution - number of points in the plot

%fig1
[Efermi, E0, E1, Volts ] = Ef_E1_E2_vs_GATE_V(T, Vmin, Vmax, f);  

%fi2: Ef - E0 and Ef - E1
[Reg1_2] =  Efermi_E0_E1(Volts, Efermi, E0, E1, f);
            Region1_2_volts  = Volts(Reg1_2);
            Region1_2_Efermi = Efermi(Reg1_2);
% end        
%voltage at which BIAS is not applyied        
V0 = find(Volts >= 0);
V0 = V0(1); %index for Vgate = 0

%bound states into trig_well AIRY functions        
[El0, El1, Emax, Airy_energies, P_occupation] = airy_fuctions(Efermi, E0, E1, V0, f);
        
%approximations

approximations(Volts, Efermi, T,f);