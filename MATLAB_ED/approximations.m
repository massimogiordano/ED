function approximations(Volts, Efermi, T, f)

global P1 l0 Voff Vth Kb Reg1_2

Vth = T*Kb;
clf
hold on 
%approssimazione eq 3
apprx1 = @(Vgo, Ef) exp((Ef - l0*(P1*(Vgo - Ef)).^(2/3))/Vth);
apprx = apprx1(Volts-Voff,Efermi);

rel_err_eq3 = (log(apprx+1)-log(apprx))./log(apprx+1);






ns = ns_eq_2(Volts- Voff, Efermi);
Ef_eq_3 = eq_3(Volts - Voff,ns);
plot(Volts,Ef_eq_3,'LineWidth',1)
hold on
plot(Volts,Efermi,'-','LineWidth',2)

plot(Volts,-(Ef_eq_3-Efermi)./abs(Efermi),'LineWidth',1)
yL = get(gca,'YLim');
line([Volts(Reg1_2) Volts(Reg1_2)],yL,'Color','k','LineStyle','-.');
%axis([-3 -1 -0.1 0.3])
    %%IMAGE SET UP
    title(['Relative Error on Fermi Energy due to Approximation'],'FontSize',15)
    xlabel('Gate Voltage [V]','FontSize',15) 
    ylabel('Energy [V]','FontSize',15)
    h_legend=legend('E approximated','E Fermi','Relative Error');
    set(h_legend,'Location','northwest','FontSize',14);

    saveas(f,['RelativeError_eq3'],'svg');
    

output_args = 1;

end

function ns = ns_eq_2(Volt, Ef)
   global q ep d
   ns = ep/(d*q)*(Volt - Ef);
end

function Ef = eq_3(Vg, ns)
 global l0 Vth D
 Ef = l0*(ns).^(2/3) + Vth*log(ns./(D*Vth));
end

