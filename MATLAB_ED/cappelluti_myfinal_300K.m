%parametri
Vmin = -2.95;
Vmax = 3;   %gate voltage

T = 300; %temperature
d = 20e-9;

fake=2.0105e-12;

Voff= -3;
l0 = 2.12e-12;
hbar= 1.055e-34
l1 = 3.73e-12; %experimental parameters
ep = 9*8.85e-12; %permittivity of AlGaN
D = 1.001e18; %density of states
q = 1.60e-19; %charge of electron in Coulomb
Vth = T*(8.61673324e-5); %bolzamn costant * T
P2 = D*Vth;
P1 = ep/(q*d);
Cg=ep/d;
beta=Cg/(q*D*Vth);
alfad= 1/beta;
alfan= exp(1)/beta; %SONO RINCOGLIONITO


DeltaV= 100;

Efermi= 1:DeltaV;
Volts = 1:DeltaV;
E1 = 1:DeltaV;
E0 = 1:DeltaV;
NS = 1:DeltaV;
%err1 = 1:DeltaV;
%EF1 = 1:DeltaV;
EF1_corr = 1:DeltaV;
EF1_1 = 1:DeltaV;
NS1 = 1:DeltaV;
EF2 = 1:DeltaV;
EF2_corr = 1:DeltaV;
%err2 = 1:DeltaV;
NSu = 1:DeltaV;
diff1 = 1:DeltaV;
diff2 = 1:DeltaV;
diffU = 1:DeltaV;
FERMI=1:DeltaV;
subb0=1:DeltaV;
subb1=1:DeltaV;
deltaf=1:DeltaV;
delta1=1:DeltaV;
delta0=1:DeltaV;
gate_bias=1:DeltaV;
NS_noVTH=1:DeltaV;
EF_noVTH=1:DeltaV;
regione1_2=1:DeltaV;
i = 0;
for Vg = Vmin:(abs(Vmin-Vmax)/DeltaV):Vmax
    i = i+1;
    Vgo = Vg - Voff;
    gate_bias(i)=Vg;
    Vgod= Vgo*alfad./(sqrt((Vgo.^2)+(alfad^2)));
Vgon= Vgo*alfan./(sqrt((Vgo.^2)+(alfan^2)));
    eq_1 = @(Ef) -P1*(Vgo - Ef) + P2*(log(exp((Ef - l0*(P1*(Vgo - Ef)).^(2/3))/Vth) + 1) + log(exp((Ef - l1*(P1*(Vgo - Ef)^(2/3)))/(Vth)) + 1 )    );                        
    Efermi(i) = fzero(eq_1,0);
    Volts(i) = Vg; 
    ns = @(Vg, Ef) P1*(Vg - Ef);
    E1(i) = l1*ns(Vgo,Efermi(i))^(2/3);
    E0(i) = l0*ns(Vgo,Efermi(i))^(2/3);
      NS(i) = ns(Vgo,Efermi(i));
      regione1_2(i)=Efermi(i)-E0(i);
      %regione1
   %EF1(i) = l0*(NS(i)^(2/3))+Vth*log(NS(i)/(D*Vth)) ; %formula per regione 1 ma valori NS generale!! poi confronto
    EF1_1(i)  = Vgo*(Vth*log(beta*Vth)+l0*((Cg*Vgo/q)^(2/3)))/(Vgo+Vth+2*l0*((Cg*Vgo/q)^(2/3))) ;
    EF_noVTH(i)  = Vgo*(l0*((Cg*Vgo/q)^(2/3)))/(Vgo+2*l0*((Cg*Vgo/q)^(2/3))) ;
    NS1(i) = (Cg*Vgo/q)*(Vgo+Vth-Vth*log(beta*Vgo)-l0*((Cg*Vgo/q)^(2/3))/3)/(Vgo+Vth+2*l0*((Cg*Vgo/q)^(2/3))/3);
    NS_noVTH(i) = (Cg*Vgo/q)*(Vgo-l0*((Cg*Vgo/q)^(2/3))/3)/(Vgo+2*l0*((Cg*Vgo/q)^(2/3))/3);
   
    %calcolo E_fermi con 'corretta' NS del modello1
    EF1_corr(i) = l0*(NS1(i)^(2/3))+Vth*log(NS1(i)/(D*Vth))  ;
    %err1(i)=abs((EF1_corr-EF1)/EF1_corr);
    diff1(i)=100*abs((NS1(i)-NS(i))/NS(i));
    %====regione2
   %EF2(i)= ; %formula regione 2 con valori stimati da caso generale! poi confronto
   NS2(i) = ((Cg*Vgo/q)*(Vgo-l0*((Cg*Vgo/q)^(2/3))/3))/(Vgo+(Vgo*beta*Vth)+2*l0*((Cg*Vgo/q)^(2/3))/3) ;
   EF2_corr(i) = Vgo*(beta*Vth*Vgo+l0*((Cg*Vgo/q)^(2/3)))/(Vgo+Vgo*beta*Vth+2*l0*((Cg*Vgo/q)^(2/3))/3) ;
   %err2(i)=;
   diff2(i)= 100*abs((NS2(i)-NS(i))/NS(i));
   %=====UNIFIED MODEL
  NSu(i) = (Cg*Vgo/q)*(Vgo+Vth-Vth*log(beta*Vgon)-l0*(((Cg*Vgo/q)^(2/3))/3))/(Vgo+(Vgo*Vth/Vgod)+2*l0*(((Cg*Vgo/q)^(2/3))/3)) ;
% NSu(i) = (Cg*Vgo/q)*(Vgo-l0*(((Cg*Vgo/q)^(2/3))/3))/(Vgo+2*l0*(((Cg*Vgo/q)^(2/3))/3))  
diffU(i)= 100*abs((NSu(i)-NS(i))/NS(i));

FERMI(i)= Vgo - NSu(i)/P1;
subb0(i)=l0*(NSu(i)^(2/3));
subb1(i)=l1*(NSu(i)^(2/3));
deltaf(i)= 100*abs((Efermi(i)-FERMI(i))/Efermi(i));
delta0(i)= 100*abs((E0(i)-subb0(i))/E0(i));
delta1(i)= 100*abs((E1(i)-subb1(i))/E1(i));
end
%--------------
Eturn=find(Efermi-E0>=0);
Eturn=Eturn(1);
tmp=Volts(Eturn) %stampa valore gate bia per cui Ef-E0 diventa positivo
%--------------
%--------
fixed=find(Volts-Voff>=0);
fixed=fixed(1);
tmp1=NSu(fixed);
tmp2=NS(fixed);
tt= (tmp1-tmp2)/ tmp2
save k300.txt tmp1 tmp2 tt -ASCII;
NSut=NSu.'
save nsU_300k.txt NSut -ASCII; 
%--------------
% plot(Volts, Efermi)
% 
% hold on
% title('Fermi energy, first and second subbands vs. Vg')
% xlabel('Vg(V)') % x-axis label
% ylabel('Ef, E0, E1 (V) ') % y-axis label
% plot(Volts, E1)
% plot(Volts, E0)
% grid on
% figure %stampa in nuova finestra
% 
% plot(Volts,Efermi-E0)
% hold on
% title('energy difference vs. Vg')
% xlabel('Vg(V)') % x-axis label
% ylabel('Ef-E0,Ef-E1 (V) ') % y-axis label
% 
% plot(Volts,Efermi-E1)
% grid on
% %legend('y = sin(x)','y = cos(x)','Location','southwest') METTERE LEGENDA
% % plot regione 1
% %figure
% %plot(NS,EF1)
% %grid on
% % figure
% % plot(NS1,log10(EF1_corr))
% % grid on
% % %semilogy(NS1,log10(EF1_corr))
% % figure
% % plot(Volts,log10(EF1_1))
% % %semilogy(Volts,EF1_1)
% % grid on
% % %figure
% % %plot(Volts,err1)
% % %grid on
% % %figure
% % %plot(Volts,log10(NS1))
% % 
% % %grid on
% % %plot regione 2
% % %figure
% % %plot(NS,EF2)
% % %grid on
% % figure
% % plot(NS2,log10(EF2_corr))
% % %semilogy(NS2,log10(EF2_corr))
% % grid on
% % figure
% % plot(Volts,log10(EF2_corr))
% % %semilogy(Volts,EF2_corr)
% % grid on
% 
% 
% 
% %grid on
% figure
% plot(Volts,log10(NS1))
% grid on
% hold on
% plot(Volts,log10(NS2))
% plot(Volts,log10(NSu))
% plot(Volts,log10(NS),'o')
% title('Ns vs.  ')
% xlabel('Vg(V)') % x-axis label
% ylabel('Ns / m^{-2}') % y-axis label
% figure
% plot(Volts,diff1, 'LineWidth',2)
% grid on
% hold on
% plot(Volts,diff2, 'LineWidth',2)
% 
% plot(Volts,diffU, 'LineWidth',2)
% title(['Relative error modulus of regional and unified model n_s vs V_g'],'FontSize',15)
%     xlabel('Vg (V)','FontSize',15) 
%     ylabel('Relative error modulus (%)','FontSize',15)
%     h_legend=legend('n_s^I','n_s^{II}','n_s^u');
%     set(h_legend,'Location','southeast','FontSize',14);
%     axis([-3 3 0 10])
% %fplot(eq_1)
% figure
% plot(Volts,log10(NSu),'square')
% hold on
% plot(Volts,log10(NS1),'o')
% 
% plot(Volts,log10(NS2),'*')
% plot(Volts,log10(NS),'Linewidth',2)
% grid on
% title(['n_s vs V_g'],'FontSize',15)
%     xlabel('Vg (V)','FontSize',15) 
%     ylabel('n_s(m^{-2})','FontSize',15)
%     h_legend=legend('n_s^u','n_s^{I}','n_s^{II}','n_s^N');
%     set(h_legend,'Location','southeast','FontSize',14);
%     
% grid on
% figure
% plot(Volts, Efermi./(Volts-Voff))
% grid on
% figure
% plot(Volts,log10(NS_noVTH))
% hold on
% plot(Volts,log10(NS1))
% 
% plot(Volts,log10(NS2))
% grid on
% figure
% plot(Volts,abs((NS1-NS_noVTH)./NS_noVTH))
% grid on
% figure
% plot(Volts,abs((NS2-NS_noVTH)./NS_noVTH))
% 
% grid on
% figure
% plot(Volts,(NS1-NS),'LineWidth',2)
% hold on
% plot(Volts,(NS2-NS),'LineWidth',2)
% grid on
% title(['Difference between regional and numerical n_s vs V_g'],'FontSize',15)
%     xlabel('Vg (V)','FontSize',15) 
%     ylabel('(n_s^I - n_s^u), (n_s^{II} - n_s^u)      (m^{-2})','FontSize',15)
%     h_legend=legend('n_s^I-n_s^u','n_s^{II}-n_s^u');
%     set(h_legend,'Location','southeast','FontSize',14);
%     
% figure
% plot(Volts,(NS-NSu))
% grid on
% figure
% plot(Volts,100*abs((NS-NSu)./NS))
% grid on
% % figure
% % 
% % plot(Volts, FERMI)
% % 
% % hold on
% % title('Vg(V) vs. Fermi energy, first and second subbands')
% % xlabel('Vg(V)') % x-axis label
% % ylabel('Ef, E0, E1 (V) ') % y-axis label
% % plot(Volts, subb0)
% % plot(Volts, subb1)
% % grid on
% % figure
% % plot(Volts,deltaf)
% % grid on
% % hold on
% % plot(Volts,delta0)
% % 
% % plot(Volts,delta1)
% % NS_t=NS.';
% % Volts_t=Volts.';
% % Efermi_t=Efermi.';
% % discr=abs(NSu-NS);
% % discr_t=discr.';
% % media_err=sum(discr)/size(NSu,2);
% % i = 0;
% % veterr=1:DeltaV;
% % for Vg = Vmin:(abs(Vmin-Vmax)/DeltaV):Vmax
% %     i = i+1;
% %     veterr(i)=media_err;
% % end
% 
%     
E0nsfake=fake.*NS.^(2/3);
figure
plot(Volts,100*abs(E0nsfake-E0)./E0)
fake
E1nsfake=fake2.*NS.^(2/3);
figure
plot(Volts,100*abs(E1nsfake-E1)./E1)
