clear
clc
Sbase=10000;%VA
fbase=50;
Vbase=381;%V
Zbase=Vbase^2/Sbase;
L1=0.0262;
R1=0.0103;
L2=0.0262;
R2=0.0358;
% 
%  Z1=R1+i*L1*w;
%  Z2=R2+i*L2*w;
 %parameter
 Dw=100;
 Dv=20;
 f=50;

 P4=0.42;
 P5=0.42;
 Q4=0;
 Q5=0;
 
 w0=2*pi*f;



%theta1 as reference angle
% f=@(V1,V2,V3,theta1,theta2,theta3,w) [V1*V2*(real(-1/(R1+i*L1*w))*cos(theta1-theta2)+imag(-1/(R1+i*L1*w))*sin(theta1-theta2))+V1*V1*real(1/(R1+i*L1*w))+P4-Dw*(1-w);
%     V1*V2*(real(-1/(R1+i*L1*w))*cos(theta2-theta1)+imag(-1/(R1+i*L1*w))*sin(theta2-theta1))+V2^2*real(1/(R1+i*L1*w)+1/(R2+i*L2*w))+V2*V3*(real(-1/(R2+i*L2*w))*cos(theta2-theta3)+imag(-1/(R2+i*L2*w))*sin(theta2-theta3))-Dw*(1-w);
%     V2*V3*(real(-1/(R2+i*L2*w))*cos(theta3-theta2)+imag(-1/(R2+i*L2*w))*sin(theta3-theta2))+V3^2*real(1/(R2+i*L2*w))+P5-Dw*(1-w);
%     theta1;
%     V1*V2*(real(-1/(R1+i*L1*w))*sin(theta1-theta2)-imag(-1/(R1+i*L1*w))*cos(theta1-theta2))-V1^2*imag(1/(R1+i*L1*w))+Q4-Dv*(1-V1);
%     V1*V2*(real(-1/(R1+i*L1*w))*sin(theta2-theta1)-imag(-1/(R1+i*L1*w))*cos(theta2-theta1))-V2^2*imag(1/(R1+i*L1*w)+1/(R2+i*L2*w))+V2*V3*(real(-1/(R2+i*L2*w))*sin(theta2-theta3)-imag(-1/(R2+i*L2*w))*cos(theta2-theta3))-Dv*(1-V2);
%     V2*V3*(real(-1/(R2+i*L2*w))*sin(theta3-theta2)-imag(-1/(R2+i*L2*w))*cos(theta3-theta2))-V3^2*imag(1/(R2+i*L2*w))+Q5-Dv*(1-V3) ];
%% after
f=@(V1,V2,V3,theta1,theta2,theta3,w) [V1*V2*(real(-1/(R1+i*L1*w))*cos(theta1-theta2)+imag(-1/(R1+i*L1*w))*sin(theta1-theta2))+V1*V1*real(1/(R1+i*L1*w))+P4;
    V1*V2*(real(-1/(R1+i*L1*w))*cos(theta2-theta1)+imag(-1/(R1+i*L1*w))*sin(theta2-theta1))+V2^2*real(1/(R1+i*L1*w)+1/(R2+i*L2*w))+V2*V3*(real(-1/(R2+i*L2*w))*cos(theta2-theta3)+imag(-1/(R2+i*L2*w))*sin(theta2-theta3))-Dw*(1-w);
    V2*V3*(real(-1/(R2+i*L2*w))*cos(theta3-theta2)+imag(-1/(R2+i*L2*w))*sin(theta3-theta2))+V3^2*real(1/(R2+i*L2*w))+P5-Dw*(1-w);
    theta1;
    V1*V2*(real(-1/(R1+i*L1*w))*sin(theta1-theta2)-imag(-1/(R1+i*L1*w))*cos(theta1-theta2))-V1^2*imag(1/(R1+i*L1*w))+Q4;
    V1*V2*(real(-1/(R1+i*L1*w))*sin(theta2-theta1)-imag(-1/(R1+i*L1*w))*cos(theta2-theta1))-V2^2*imag(1/(R1+i*L1*w)+1/(R2+i*L2*w))+V2*V3*(real(-1/(R2+i*L2*w))*sin(theta2-theta3)-imag(-1/(R2+i*L2*w))*cos(theta2-theta3))-Dv*(1-V2);
    V2*V3*(real(-1/(R2+i*L2*w))*sin(theta3-theta2)-imag(-1/(R2+i*L2*w))*cos(theta3-theta2))-V3^2*imag(1/(R2+i*L2*w))+Q5-Dv*(1-V3) ];
%  V1,V2,V3,theta1,theta2,theta3,w

fp=@(x) f(x(1),x(2),x(3),x(4),x(5),x(6),x(7));
 
%[x, fval, info] = fsolve (fp, [381;381;381;0;0;0;2*pi*50]);
 [x, fval, info] = fsolve (fp, [1;1;1;0;0;0;1]);
%disp(x);

 Ibase=Sbase/Vbase/3^0.5;
%  disp(x(1)*381);
% disp(x(2)*381);
% disp(x(3)*381);
% % disp("theta1: "+x(4)+"    theta1 in angle: "+x(4)*180/pi);
% % disp("theta2: "+x(5)+"    theta2 in angle: "+x(5)*180/pi);
% % disp("theta3: "+x(6)+"    theta3 in angle: "+x(6)*180/pi);
% disp(x(7)*w0/2/pi);
% disp(Dw*(1-x(7))*Sbase);
% 
% disp(Dv*(1-x(1))*Sbase);
% disp(Dv*(1-x(2))*Sbase);
% disp(Dv*(1-x(3))*Sbase);

 Z1=R1+i*L1*x(7);
 Z2=R2+i*L2*x(7);
disp("V1: "+x(1)+" V"+x(1)*381);
disp("V2: "+x(2)+" V"+x(2)*381);
disp("V3: "+x(3)+" V"+x(3)*381);
disp("theta1: "+x(4)+"    theta1 in angle: "+x(4)*180/pi);
disp("theta2: "+x(5)+"    theta2 in angle: "+x(5)*180/pi);
disp("theta3: "+x(6)+"    theta3 in angle: "+x(6)*180/pi);
disp("omega: "+x(7)*w0+"   omega infrequency: "+x(7)*w0/2/pi);
disp("P1 setpoint "+Dw*(1-x(7))*Sbase+" W");
disp("P2 setpoint "+Dw*(1-x(7))*Sbase+" W");
disp("P3 setpoint: "+Dw*(1-x(7))*Sbase+" W");
disp("Q1 setpoint: "+Dv*(1-x(1))*Sbase+" Var");
disp("Q2 setpoint: "+Dv*(1-x(2))*Sbase+" Var");
disp("Q3 setpoint: "+Dv*(1-x(3))*Sbase+" Var");
disp("I12 current: "+(x(1)*exp(i*x(4))-x(2)*exp(i*x(5)))/Z1+" pu   A"+(x(1)*exp(i*x(4))-x(2)*exp(i*x(5)))/Z1*Ibase);
disp("I23 current: "+(x(2)*exp(i*x(5))-x(3)*exp(i*x(6)))/Z2+" pu   A"+(x(2)*exp(i*x(5))-x(3)*exp(i*x(6)))/Z2*Ibase);
disp("I23 current: "+(x(2)*exp(i*x(5))-x(3)*exp(i*x(6)))/Z2+" A");