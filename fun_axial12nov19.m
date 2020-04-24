function dydt = fun_axial12nov19(~,y,stiffness,mass,drag,ODm,holeD,Kc,rho,k,h1,n,ODj,h2,...
    ODc,hC,ODb,hB,aa,bb,Au1,Au2,AuC,AuB,F_hook)
global Ff FFtrack
v1 = [];
v2 = [];
Re1= [];
Re2= [];
for i=4:2:numel(y)
    if i < numel(y)
        ve = (-(ODm^2)/(holeD^2 - ODm^2))*(y(i)) + Kc*(y(i));
        v1 = [v1 ve];
        Rem = (8*rho/k)*(((h1*n)/(4*n+2))^n)*(abs(ve)^(2-n));
        Re1 = [Re1 Rem];
        
        vse= (-(ODj^2)/(holeD^2 - ODj^2))*(y(i)) + Kc*(y(i));
        v2 = [v2 vse];
        Res = (8*rho/k)*(((h2*n)/(4*n+2))^n)*(abs(vse)^(2-n));
        Re2 = [Re2 Res];

    else
        ve = (-(ODc^2)/(holeD^2 - ODc^2))*(y(i)) + Kc*(y(i));
        v1 = [v1 ve];
        Rem = (8*rho/k)*(((hC*n)/(4*n+2))^n)*(abs(ve)^(2-n));
        Re1 = [Re1 Rem];

        vse= (-(ODb^2)/(holeD^2 - ODb^2))*(y(i)) + Kc*(y(i));
        v2 = [v2 vse];
        Res = (8*rho/k)*(((hB*n)/(4*n+2))^n)*(abs(vse)^(2-n));
        Re2 = [Re2 Res];
    end
end

f1 = aa./(Re1.^bb);
f2 = aa./(Re2.^bb);

%Model Turbulent in annulus
Tau1 = -0.5*f1*rho.*abs(v1).*v1;
Tau2 = -0.5*f2*rho.*abs(v2).*v2;

% Laminar in Drillstring
% TauIn1 = fun_laminarflow(n,k,tauY,vIn1,IDm);
% TauIn2 = fun_laminarflow(n,k,tauY,vIn2,IDrat2);

Ff = Au1*Tau1 + Au2*Tau2;
Ff(end) = AuC*Tau1(end) + AuB*Tau2(end);%It must be different 
FFtrack = [FFtrack;Ff];

dydt =[ y(2)
    F_hook %(F_hook - stiffness(1)*(y(1)-y(3)))/mass(1)
    y(4)
    (stiffness(1)*(y(1) -y(3)) - stiffness(2)*(y(3)-y(5)) - (drag(1)+Ff(1)))/mass(2)
    y(6)
    (stiffness(2)*(y(3)-y(5)) - stiffness(3)*(y(5)-y(7)) - (drag(2)+Ff(2)))/mass(3)
    y(8)
    (stiffness(3)*(y(5)-y(7)) - stiffness(4)*(y(7)-y(9)) - (drag(3)+Ff(3)))/mass(4)
    y(10)
    (stiffness(4)*(y(7)-y(9)) - (drag(4)+Ff(4)))/mass(5)];
