function [Efermi, E0, E1, Volts ] = Ef_E1_E2_vs_GATE_V(T, Vmin, Vmax,  f)

global Voff l0 l1 Vth P2 P1 DeltaV Vgo Kb ep q d D 

Vth = T*Kb; %bolzamn costant * T
P2 = D*Vth;
P1 = ep/(q*d);



Efermi= 1:DeltaV;
Volts = 1:DeltaV;
E1 = 1:DeltaV;
E0 = 1:DeltaV;

clf
    
    i=0;
    for Vg = Vmin:(abs(Vmin-Vmax)/DeltaV):Vmax
        i = i+1;
        
        Vgo = Vg - Voff;
        Volts(i) = Vg;
        
        eq_1 = @(Ef) -P1*(Vgo - Ef) + P2*(log(exp((Ef - l0*(P1*(Vgo - Ef)).^(2/3))./Vth) + 1) + log(exp((Ef - l1*(P1*(Vgo - Ef).^(2/3)))./(Vth)) + 1 )    );                        
                               

        Efermi(i) = fzero(eq_1,0);
        

        ns = @(Vg, Ef) P1*(Vg - Ef);

        %two bands in juction well as define in the paper
        E1(i) = l1*ns(Vgo,Efermi(i))^(2/3);
        E0(i) = l0*ns(Vgo,Efermi(i))^(2/3);
         

    end


    %figure [1] - Energies vs. Voltages
    fig1_EvsVolts(Volts,Efermi, E0, E1, T ,f)
end
