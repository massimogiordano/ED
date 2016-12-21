function fig1_EvsVolts(Volts, Efermi, E0, E1, T,f)
plot(Volts, Efermi,'LineWidth',2)
hold on
plot(Volts, E1,'LineWidth',2)
plot(Volts, E0,'LineWidth',2)
    title(['Fermi energy, E_0 and E_1 subbands' ' at ' num2str(T) 'K'],'FontSize',15)
    xlabel('Gate Voltage [V]','FontSize',15) 
    ylabel('Energy [eV]','FontSize',15)
    h_legend=legend('E_{fermi}','E_1','E_0');
    set(h_legend,'Location','southeast','FontSize',14);
    axis([-3 3 -0.15 1])
saveas(f,['Energies_vs_Volts' num2str(T) '.svg'],'svg');

