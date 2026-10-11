function [f,A,B]=fw_motion(x,u)
% 7-state coordinated point-mass model; SI units and radians.
% x=[X Y Z V gamma chi mu]'; u=[a_t n omega_mu]'.
V=x(4);ga=x(5);ch=x(6);mu=x(7);at=u(1);n=u(2);om=u(3);g=9.81;
f=[V*cos(ga)*cos(ch);V*cos(ga)*sin(ch);V*sin(ga);at-g*sin(ga);...
   g/V*(n*cos(mu)-cos(ga));g*n*sin(mu)/(V*cos(ga));om];
A=zeros(7);B=zeros(7,3);
A(1,4)=cos(ga)*cos(ch);A(1,5)=-V*sin(ga)*cos(ch);A(1,6)=-V*cos(ga)*sin(ch);
A(2,4)=cos(ga)*sin(ch);A(2,5)=-V*sin(ga)*sin(ch);A(2,6)=V*cos(ga)*cos(ch);
A(3,4)=sin(ga);A(3,5)=V*cos(ga);A(4,5)=-g*cos(ga);B(4,1)=1;
A(5,4)=-f(5)/V;A(5,5)=g/V*sin(ga);A(5,7)=-g*n/V*sin(mu);B(5,2)=g/V*cos(mu);
A(6,4)=-f(6)/V;A(6,5)=f(6)*tan(ga);A(6,7)=g*n*cos(mu)/(V*cos(ga));B(6,2)=g*sin(mu)/(V*cos(ga));B(7,3)=1;
end
