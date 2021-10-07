function DF = DF_realified(in1,in2)
%DF_REALIFIED
%    DF = DF_REALIFIED(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    06-Oct-2021 20:26:00

X1 = in1(1,:);
X2 = in1(2,:);
Xprime1 = in1(3,:);
Xprime2 = in1(4,:);
cspeed = in1(5,:);
lamprime1 = in1(9,:);
lamprime2 = in1(10,:);
omega = in1(6,:);
para1 = in2(1,:);
para2 = in2(2,:);
para3 = in2(3,:);
v1 = in1(7,:);
v2 = in1(8,:);
vprime1 = in1(11,:);
vprime2 = in1(12,:);
DF = reshape([0.0,X1.*para1+para1.*(X1-1.0),0.0,Xprime1.*para1.*2.0,0.0,para1+para1.*((para2.*sin(omega))./omega-(sin(omega.*para3).*(para2-1.0))./(omega.*para3)),0.0,-para1.*((para2.*(cos(omega)-1.0))./omega+((cos(omega.*para3)-1.0).*(para2-1.0))./(omega.*para3)),0.0,lamprime1.*(para1+para1.*(-1.0./omega.^2.*para2.*(cos(omega)+omega.*sin(omega)-1.0)+(1.0./omega.^2.*cos(omega.*para3).*(cos(omega.*para3)-1.0).*(para2-1.0))./para3+(1.0./omega.^2.*sin(omega.*para3).*(sin(omega.*para3)-omega.*para3).*(para2-1.0))./para3))+lamprime2.*para1.*(1.0./omega.^2.*para2.*(sin(omega)-omega.*cos(omega))+(1.0./omega.^2.*sin(omega.*para3).*(cos(omega.*para3)-1.0).*(para2-1.0))./para3-(1.0./omega.^2.*cos(omega.*para3).*(sin(omega.*para3)-omega.*para3).*(para2-1.0))./para3),0.0,lamprime2.*(para1+para1.*(-1.0./omega.^2.*para2.*(cos(omega)+omega.*sin(omega)-1.0)+(1.0./omega.^2.*cos(omega.*para3).*(cos(omega.*para3)-1.0).*(para2-1.0))./para3+(1.0./omega.^2.*sin(omega.*para3).*(sin(omega.*para3)-omega.*para3).*(para2-1.0))./para3))-lamprime1.*para1.*(1.0./omega.^2.*para2.*(sin(omega)-omega.*cos(omega))+(1.0./omega.^2.*sin(omega.*para3).*(cos(omega.*para3)-1.0).*(para2-1.0))./para3-(1.0./omega.^2.*cos(omega.*para3).*(sin(omega.*para3)-omega.*para3).*(para2-1.0))./para3),1.0,cspeed,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-para1+X1.*para1.*2.0,0.0,0.0,0.0,0.0,0.0,para1.*((para2.*sin(omega))./omega+omega.*para3.*sin(omega.*para3).*(para2-1.0)+1.0),0.0,-para1.*((para2.*(cos(omega)-1.0))./omega-omega.*para3.*(cos(omega.*para3)-1.0).*(para2-1.0)),0.0,0.0,1.0,cspeed,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,X2,0.0,Xprime2,0.0,v1,0.0,v2,0.0,vprime1+lamprime1.*v1-lamprime2.*v2,0.0,vprime2+lamprime1.*v2+lamprime2.*v1,0.0,0.0,0.0,0.0,0.0,v2+X1.*para1.*((para2.*cos(omega))./omega-1.0./omega.^2.*para2.*sin(omega)-(cos(omega.*para3).*(para2-1.0))./omega+(1.0./omega.^2.*sin(omega.*para3).*(para2-1.0))./para3),-1.0,-v1+X1.*para1.*((para2.*sin(omega))./omega+(sin(omega.*para3).*(para2-1.0))./omega+1.0./omega.^2.*para2.*(cos(omega)-1.0)+(1.0./omega.^2.*(cos(omega.*para3)-1.0).*(para2-1.0))./para3),0.0,vprime2+Xprime1.*para1.*(para3.*sin(omega.*para3).*(para2-1.0)+(para2.*cos(omega))./omega-1.0./omega.^2.*para2.*sin(omega)+omega.*para3.^2.*cos(omega.*para3).*(para2-1.0))+X1.*lamprime2.*para1.*(-1.0./omega.^2.*sin(omega.*para3).^2.*(para2-1.0)+(para2.*sin(omega))./omega-1.0./omega.^3.*para2.*(sin(omega)-omega.*cos(omega)).*2.0+1.0./omega.^2.*sin(omega.*para3).*(sin(omega.*para3)-omega.*para3).*(para2-1.0)+1.0./omega.^2.*cos(omega.*para3).*(cos(omega.*para3)-1.0).*(para2-1.0)-(1.0./omega.^3.*sin(omega.*para3).*(cos(omega.*para3)-1.0).*(para2-1.0).*2.0)./para3+(1.0./omega.^2.*cos(omega.*para3).*(para3-para3.*cos(omega.*para3)).*(para2-1.0))./para3+(1.0./omega.^3.*cos(omega.*para3).*(sin(omega.*para3)-omega.*para3).*(para2-1.0).*2.0)./para3)-X1.*lamprime1.*para1.*((para2.*cos(omega))./omega-1.0./omega.^3.*para2.*(cos(omega)+omega.*sin(omega)-1.0).*2.0+1.0./omega.^2.*cos(omega.*para3).*sin(omega.*para3).*(para2-1.0)+1.0./omega.^2.*sin(omega.*para3).*(cos(omega.*para3)-1.0).*(para2-1.0)-1.0./omega.^2.*cos(omega.*para3).*(sin(omega.*para3)-omega.*para3).*(para2-1.0)+(1.0./omega.^3.*cos(omega.*para3).*(cos(omega.*para3)-1.0).*(para2-1.0).*2.0)./para3+(1.0./omega.^2.*sin(omega.*para3).*(para3-para3.*cos(omega.*para3)).*(para2-1.0))./para3+(1.0./omega.^3.*sin(omega.*para3).*(sin(omega.*para3)-omega.*para3).*(para2-1.0).*2.0)./para3),0.0,-vprime1+Xprime1.*para1.*((para2.*sin(omega))./omega+para3.*(cos(omega.*para3)-1.0).*(para2-1.0)+1.0./omega.^2.*para2.*(cos(omega)-1.0)-omega.*para3.^2.*sin(omega.*para3).*(para2-1.0))-X1.*lamprime1.*para1.*(-1.0./omega.^2.*sin(omega.*para3).^2.*(para2-1.0)+(para2.*sin(omega))./omega-1.0./omega.^3.*para2.*(sin(omega)-omega.*cos(omega)).*2.0+1.0./omega.^2.*sin(omega.*para3).*(sin(omega.*para3)-omega.*para3).*(para2-1.0)+1.0./omega.^2.*cos(omega.*para3).*(cos(omega.*para3)-1.0).*(para2-1.0)-(1.0./omega.^3.*sin(omega.*para3).*(cos(omega.*para3)-1.0).*(para2-1.0).*2.0)./para3+(1.0./omega.^2.*cos(omega.*para3).*(para3-para3.*cos(omega.*para3)).*(para2-1.0))./para3+(1.0./omega.^3.*cos(omega.*para3).*(sin(omega.*para3)-omega.*para3).*(para2-1.0).*2.0)./para3)-X1.*lamprime2.*para1.*((para2.*cos(omega))./omega-1.0./omega.^3.*para2.*(cos(omega)+omega.*sin(omega)-1.0).*2.0+1.0./omega.^2.*cos(omega.*para3).*sin(omega.*para3).*(para2-1.0)+1.0./omega.^2.*sin(omega.*para3).*(cos(omega.*para3)-1.0).*(para2-1.0)-1.0./omega.^2.*cos(omega.*para3).*(sin(omega.*para3)-omega.*para3).*(para2-1.0)+(1.0./omega.^3.*cos(omega.*para3).*(cos(omega.*para3)-1.0).*(para2-1.0).*2.0)./para3+(1.0./omega.^2.*sin(omega.*para3).*(para3-para3.*cos(omega.*para3)).*(para2-1.0))./para3+(1.0./omega.^3.*sin(omega.*para3).*(sin(omega.*para3)-omega.*para3).*(para2-1.0).*2.0)./para3),0.0,0.0,0.0,0.0,1.0,cspeed,0.0,-omega,lamprime1,-lamprime1+cspeed.*lamprime1+1.0,lamprime2,-lamprime2+cspeed.*lamprime2,0.0,0.0,0.0,0.0,0.0,omega,1.0,cspeed,-lamprime2,lamprime2-cspeed.*lamprime2,lamprime1,-lamprime1+cspeed.*lamprime1+1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,v1-1.0,-v1+cspeed.*v1+para1.*(X1-1.0)+X1.*para1.*(-1.0./omega.^2.*para2.*(cos(omega)+omega.*sin(omega)-1.0)+(1.0./omega.^2.*cos(omega.*para3).*(cos(omega.*para3)-1.0).*(para2-1.0))./para3+(1.0./omega.^2.*sin(omega.*para3).*(sin(omega.*para3)-omega.*para3).*(para2-1.0))./para3),v2,-v2+cspeed.*v2-X1.*para1.*(1.0./omega.^2.*para2.*(sin(omega)-omega.*cos(omega))+(1.0./omega.^2.*sin(omega.*para3).*(cos(omega.*para3)-1.0).*(para2-1.0))./para3-(1.0./omega.^2.*cos(omega.*para3).*(sin(omega.*para3)-omega.*para3).*(para2-1.0))./para3),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-v2,v2-cspeed.*v2+X1.*para1.*(1.0./omega.^2.*para2.*(sin(omega)-omega.*cos(omega))+(1.0./omega.^2.*sin(omega.*para3).*(cos(omega.*para3)-1.0).*(para2-1.0))./para3-(1.0./omega.^2.*cos(omega.*para3).*(sin(omega.*para3)-omega.*para3).*(para2-1.0))./para3),v1-1.0,-v1+cspeed.*v1+para1.*(X1-1.0)+X1.*para1.*(-1.0./omega.^2.*para2.*(cos(omega)+omega.*sin(omega)-1.0)+(1.0./omega.^2.*cos(omega.*para3).*(cos(omega.*para3)-1.0).*(para2-1.0))./para3+(1.0./omega.^2.*sin(omega.*para3).*(sin(omega.*para3)-omega.*para3).*(para2-1.0))./para3),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,cspeed,0.0,-omega,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,omega,1.0,cspeed],[12,12]);