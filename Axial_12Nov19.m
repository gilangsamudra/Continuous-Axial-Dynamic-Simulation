clear all;
% close all;
clc;
global Ff FFtrack
FFtrack = [];
% F_hook = 6.8*222411;
load('Input Accel.mat');
Hook = Trajectory(2:end);

%% DRAG FORCES
Loaded_Data = readtable('Drag 1KM 4 Segments.xlsx'); % in table format
drag        = Loaded_Data{:,5};                      % drag forces (lbf)
drag        = drag(~isnan(drag));

% Panjang in meter
Total = 1000;
ns    = 4;
LengthS = Total/ns;
mod_el = 206842718795.3;          % Elastic modulus, N/m^2

% String Table, inches
ODm = 5 * 0.0254; %in meter now
ODj = 6.5 * 0.0254;
ODc = 6.5 * 0.0254;
ODb = 6.5 * 0.0254;
IDm = 4.276 * 0.0254;
IDj = 3.75 * 0.0254;
IDc = 2.813 * 0.0254;
IDb = 3.89 * 0.0254;
holeD = (7+(0/8)) * 0.0254;

% in SI unit
pipe_den = 7800.001722079;        %Pipe density, in kg/m^3
rho  = 1840;                      %Mud density, in kg/m^3

alp = ODm/holeD;
Kc = (alp^2 - sqrt((alp^4 +alp)/(1+alp)))/(1-alp^2);
n = 0.55; 
k = 0.13;
h1 = (holeD-ODm)/2; %The width of the annular space
h2 = (holeD-ODj)/2;
hC = (holeD-ODc)/2;
hB = (holeD-ODb)/2;
aa = 0.138; %Rheology parameter for annulus
bb = 0.312;
Au1 = 2*pi*0.95*LengthS*(ODm/2 + holeD/2);
Au2 = 2*pi*0.05*LengthS*(ODj/2 + holeD/2);
AuC = 2*pi*0.911*LengthS*(ODc/2 + holeD/2);
AuB = 2*pi*0.088*LengthS*(ODb/2 + holeD/2);

% Drillstring Parameters
Stiff_bha = mod_el*((ODb^2 - IDb^2)*pi/4)/(0.088*LengthS);
Stiff_coll= mod_el*((ODc^2 - IDc^2)*pi/4)/(0.911*LengthS);
Stiff_main= mod_el*((ODm^2 - IDm^2)*pi/4)/(0.95*LengthS);
Stiff_sec = mod_el*((ODj^2 - IDj^2)*pi/4)/(0.05*LengthS);

stiffness1 = Stiff_main*Stiff_sec/(Stiff_main+Stiff_sec);     %stiffness of drill pipe
stiffness2 = stiffness1;
stiffness3 = stiffness1;
stiffness4 = Stiff_bha*Stiff_coll/(Stiff_bha+Stiff_coll);
stiffness = [stiffness1 stiffness2 stiffness3 stiffness4];

mass_top  = 20000;
mass_bha  = ((ODb^2 - IDb^2)*pi/4)*(0.088*LengthS)*pipe_den;        
mass_coll = ((ODc^2 - IDc^2)*pi/4)*(0.911*LengthS)*pipe_den;
mass_main = ((ODm^2 - IDm^2)*pi/4)*(0.95*LengthS)*pipe_den;
mass_sec  = ((ODj^2 - IDj^2)*pi/4)*(0.05*LengthS)*pipe_den;

mass1 = mass_main+mass_sec; %mass of drill pipe
mass2 = mass1;
mass3 = mass1;
mass4 = mass_bha+mass_coll;
mass  = [mass_top mass1 mass2 mass3 mass4];

y0    = ones(10,1)*0.001;
kecepatan = [];
waktu = [];
for i = 1:30
% run simulation ==> solve the EOM
tspan = [0 1];
F_hook = Hook(i);
% options = odeset('OutputFcn',@outputfcn);
[t,y] = ode45(@(t,y)fun_axial12nov19(t,y,stiffness,mass,drag,ODm,holeD,Kc,rho,k,h1,n,ODj,h2,...
    ODc,hC,ODb,hB,aa,bb,Au1,Au2,AuC,AuB,F_hook), tspan, y0);
y0 = y(end,:);
kecepatan = [kecepatan;y];
end
waktu = linspace(0,30,size(kecepatan,1));

figure()
xlabel('time (s)');
ylabel('velocity (m/s)');
title('BHA');
set(gca,'fontsize',15);
hold on
plot(waktu,kecepatan(:,2));
plot(waktu,kecepatan(:,4));
plot(waktu,kecepatan(:,6));
plot(waktu,kecepatan(:,8));
plot(waktu,kecepatan(:,10));
hold off




v1 = [];
v2 = [];
Re1= [];
Re2= [];
for i=4:2:10
    if i < 10
        ve = (-(ODm^2)/(holeD^2 - ODm^2))*(kecepatan(:,i)) + Kc*(kecepatan(:,i));
        v1 = [v1 ve];
        Rem = (8*rho/k)*(((h1*n)/(4*n+2))^n)*(abs(ve).^(2-n));
        Re1 = [Re1 Rem];
        
        vse= (-(ODj^2)/(holeD^2 - ODj^2))*(kecepatan(:,i)) + Kc*(kecepatan(:,i));
        v2 = [v2 vse];
        Res = (8*rho/k)*(((h2*n)/(4*n+2))^n)*(abs(vse).^(2-n));
        Re2 = [Re2 Res];

    else
        ve = (-(ODc^2)/(holeD^2 - ODc^2))*(kecepatan(:,i)) + Kc*(kecepatan(:,i));
        v1 = [v1 ve];
        Rem = (8*rho/k)*(((hC*n)/(4*n+2))^n)*(abs(ve).^(2-n));
        Re1 = [Re1 Rem];

        vse= (-(ODb^2)/(holeD^2 - ODb^2))*(kecepatan(:,i)) + Kc*(kecepatan(:,i));
        v2 = [v2 vse];
        Res = (8*rho/k)*(((hB*n)/(4*n+2))^n)*(abs(vse).^(2-n));
        Re2 = [Re2 Res];
    end
end

f1 = aa./(Re1.^bb);
f2 = aa./(Re2.^bb);

%Model Turbulent in annulus
Tau1 = -0.5*f1*rho.*abs(v1).*v1;
Tau2 = -0.5*f2*rho.*abs(v2).*v2;

Ff = Au1*Tau1 + Au2*Tau2;
Ff(end) = AuC*Tau1(end) + AuB*Tau2(end);%It must be different 

% figure()
% curve = animatedline('Color','r','marker','o');
% for i = 1:size(Ff,1)
%     addpoints(curve,kecepatan(i,8),Ff(i,3));
%     drawnow limitrate
%     pause(0.01)
% end
% 
% 
% figure()
% for i = 1:size(Ff,1)
%     clf
%     plot(kecepatan(i,8),Ff(i,3),'ro');
%     hold on
%     plot(kecepatan(:,8),Ff(:,3));
%     drawnow limitrate
%     pause(0.01)
%     movec(i)=getframe;
% end
% title('Viscous Force')
% xlabel('Speed')
% ylabel('Force (N)')
% 



pp1 = polyfit(kecepatan(:,4),Ff(:,1),1);
pp2 = polyfit(kecepatan(:,6),Ff(:,2),1);
pp3 = polyfit(kecepatan(:,8),Ff(:,3),1);
pp4 = polyfit(kecepatan(:,10),Ff(:,4),1);

save('Viscous Force Approximation','pp1','pp2','pp3','pp4');

popo = polyval(pp3,kecepatan(:,8));
manualpopo = pp3(1) * kecepatan(:,8) + pp3(2);

figure()
plot(kecepatan(:,8),popo,'-','linewidth',2)
hold on
plot(kecepatan(:,8),Ff(:,3),'-','linewidth',2)
title('Viscous Force','FontSize', 15)
xlabel('Velocity (m/s)','FontSize', 15)
ylabel('Force (N)','FontSize', 15)
set(gca,'FontSize', 15)
legend('Linearized Viscous Force','True Viscous Force')
hold off
% 
% mywriter = VideoWriter('curve');
% mywriter.FrameRate = 20;
% open(mywriter);
% writeVideo(mywriter,movec);
% close(mywriter);