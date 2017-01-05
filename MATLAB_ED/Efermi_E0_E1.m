function [Reg1_2] = Efermi_E0_E1(Volts, Efermi, E0, E1, f)

global T


Reg1_2 = find(Efermi - E0 >= 0);
Reg1_2 = Reg1_2(1); %index for Vgate = 0 

% figure [2] - Delta Energies

clf

plot(Volts, Efermi - E1,'LineWidth',2)
hold on 
plot(Volts, Efermi - E0,'LineWidth',2)
    
      
    xL = get(gca,'XLim');
    yL = get(gca,'YLim');
    line([Volts(Reg1_2) Volts(Reg1_2)],yL,'Color','k','LineStyle','-.');
    line(xL,[0 0],'Color','k','LineStyle','-.');
    
    %%IMAGE SET UP
    title(['Sub-bands E_i relative to Fermi Energy ' ' at ' num2str(T) 'K'],'FontSize',15)
    xlabel('Gate Voltage [V]','FontSize',14) 
    ylabel(' \Delta E [V]','FontSize',14)
    h_legend=legend('E_1','E_0');
    set(h_legend,'Location','southeast','FontSize',11);
    %axis([-3 3 -0.15 1])
    
    saveas(f,['Diff_Ef_vs_E0_E1_vs_Volts' num2str(T) '.svg'],'svg');
    
end 