clearvars
close all

%conparing the three LPT trb maps below
%load Results_TwoStage_TwoStage.mat
%load Results_TwoStage_HighPqP.mat
%load Results_TwoStage_LPT01.mat
%load Results_TwoStage_HighPqP_REPEAT.mat
%load  Results_2Stage_HighPqPTrb_4_kg_s.mat
%load EXP_Results_2Stage_HighPqPTrb_4p00_kg_s.mat

load Results_2Stage_HighPqPTrb_4p00_kg_s.mat
load Haglind_data.mat

load_vector = 10:10:100;

figure(1),plot(generator_load*100, efficiency_vector,':o',load_vector,thermal_efficiency_haglind(:,2),'d');
xlabel('Percent Load');ylabel('Thermal Efficiency');legend('Modeled','GE Data');

efficiency_vector_std = interp1(generator_load*100, efficiency_vector,load_vector,'pchip','extrap');
RMS_efficiency =  sqrt(   mean(    (efficiency_vector_std'-thermal_efficiency_haglind(:,2)).^2      )        )
Percent_error_efficiency =      mean(    100*abs((efficiency_vector_std'-thermal_efficiency_haglind(:,2))./thermal_efficiency_haglind(:,2))    )

figure(2),plot(generator_load*100, mdot_vector,':',load_vector,mdot_haglind(:,2),'d');
xlabel('Percent Load');ylabel('Mass Flow Rate [kg/s]');legend('Modeled','GE Data');

mdot_vector_std = interp1(generator_load*100, mdot_vector,load_vector,'pchip','extrap');
RMS_mdot =  sqrt(mean((mdot_vector_std'-mdot_haglind(:,2)).^2))
Percent_error_mdot =      mean(    100*abs((mdot_vector_std'-mdot_haglind(:,2))./mdot_haglind(:,2))    )

figure(3),plot(generator_load*100, PR_vector,':',load_vector,PR_haglind(:,2),'d');
xlabel('Percent Load');ylabel('Pressure Ratio [None]');legend('Modeled','GE Data');

PR_vector_std = interp1(generator_load*100, PR_vector,load_vector,'pchip','extrap');
RMS_PR =  sqrt(mean((PR_vector_std'-PR_haglind(:,2)).^2))
Percent_error_PR =      mean(    100*abs((PR_vector_std'-PR_haglind(:,2))./PR_haglind(:,2))    )

figure(4),plot(generator_load*100, exh_tmptr_vector,':',load_vector,exh_tmptr_haglind(:,2),'d');
xlabel('Percent Load');ylabel('Exhaust Temperature [Deg C]');legend('Modeled','GE Data');

exh_tmptr_vector_std = interp1(generator_load*100, exh_tmptr_vector,load_vector,'pchip','extrap');
RMS_exh_tmptr =  sqrt(mean((exh_tmptr_vector_std'-exh_tmptr_haglind(:,2)).^2))
Percent_error_exh_tmptr =      mean(    100*abs((exh_tmptr_vector_std'-exh_tmptr_haglind(:,2))./exh_tmptr_haglind(:,2))    )

fuel_rate_vector = design_fuel_rate*[0.1:0.1:1]';
AFR_vector = (mdot_vector-3.5)./fuel_rate_vector; 
AFR_vector(1)=0;
figure(5),plot(generator_load*100,AFR_vector,':o');xlabel('Percent Load');ylabel('AFR [None]');

figure(6),plot(generator_load(2:end)*100,RPM_vector(2:end)-3600,'.');xlabel('Percent Load');ylabel('RPM Deviation [RPM]');



