%%
clc;
close all;
%%
%prompt2='enter the frequency of source in MHz';
%f=input(prompt2)*1000000;%freqency
f=850000000;
disp("frequency used:"+f+"Hz");
c=300000000;
lambda = c/f;
disp("wavelength:"+lambda+"m");


prompt='enter the co-ordinates of source:';
source= input(prompt);
source=source.*lambda;
%disp(source);

s=30*lambda;%distance between edge and the reference point
disp("distance between field of interest and edge:"+s);
psi=lambda;
E0=10000;%energy at the source 
k=2*pi/lambda;%k=2pi/lamda
beta0=pi/2;%considering normal incidence to edge.
edge = [0 0];%edge co-ordinates
edge2 = [-3 0];%edge taken to calculate edge unit vector
s1=norm(edge-source);%distance between source and edge
disp("distance between source and edge:"+s1);
p1=norm(edge-source);%considering rho1 and rho2 as same and equal to s1

%%
%angles of boundaries
src_unit_vector = (edge-source)/norm(edge-source);%source unit vector
edge_unit_vector = (edge-edge2)/norm(edge-edge2);%edge unit vector
incidenceangle = dot(src_unit_vector,edge_unit_vector)/(norm(src_unit_vector)*norm(edge_unit_vector));
incidenceangle = acosd(incidenceangle);%measuring the incident angle phi-
disp("incidence angle:"+incidenceangle);
RB=180-incidenceangle;%measuring the angle for reflection boundary
disp("reflection shadow boundary:"+RB);
ISB=180+incidenceangle;%measuring the angle for incident shadow boundary
disp("incidence shadow boundary:"+ISB);

%%
%filed calculation
R=[1 0;0 -1];
L=(s*s1/(s+s1))*sin(beta0)*sin(beta0);
amplitude=sqrt(s1/(s*(s+s1)));%calculation of amplitude for difracted rays
for theta=1:1:360
    if(theta<=ISB)
        x=s*cosd(theta);
        y=s*sind(theta);
        distance = sqrt((x-source(1))^2+(y-source(2))^2);%calculation of the distance of reference point from the source 
         dir(1,theta)=E0*(exp(-1i*k*distance))*(exp(-j*k*psi))*sqrt(p1*p1/((p1+distance)*(p1+distance)));
        directrays(1,theta) = 20*log10(norm(E0*(exp(-1i*k*distance))*(exp(-1i*k*psi))*sqrt(p1*p1/((p1+distance)*(p1+distance)))));%calculation of field of direct rays at point
    else
        directrays(1,theta)=0;
        dir(1,theta)=0;
    end
    
    if(theta<RB)
        
        incident = E0*(exp(-1i*k*s1))*(exp(-1i*k*psi))*sqrt(p1*p1/((p1+s1)*(p1+s1)));%field incidented at the edge(0,0) 
        x=s*cosd(theta);
        y=s*sind(theta);
        distance = sqrt((x)^2+(y)^2);%calculating the distance of the point from the source
        ref(1,theta)=conj(incident*exp(-1i*k*distance)*(exp(-j*k*psi))*sqrt(p1*p1/((p1+distance)*(p1+distance))));
        reflectedrays(1,theta)=20*log10(norm(conj(incident*exp(-1i*k*distance)*(exp(-1i*k*psi))*sqrt(p1*p1/((p1+distance)*(p1+distance))))));%calculation of field of refracted rays at point
    else
        reflectedrays(1,theta)=0;
        ref(1,theta)=0;
        
    end
    
    
    betaplus=degtorad(incidenceangle+theta);
    betaminus=degtorad(theta-incidenceangle);
    aplus=2*cos(betaplus/2)*cos(betaplus/2);
    aminus=2*cos(betaminus/2)*cos(betaminus/2);
    diffcoef=[-exp(-1i*pi/4)/(2*sqrt(2*pi*k)*sin(degtorad(beta0)))]*[[FresnelIntegral(k*L*aminus)/(cos(degtorad((theta-incidenceangle)/2)))]-[FresnelIntegral(k*L*aminus)/(cos(degtorad((theta+incidenceangle)/2)))]];
    Ei=E0*(exp(-1i*k*s1))*(exp(-1i*k*psi))*sqrt(p1*p1/((p1+s1)*(p1+s1)));%field incidented at (0,0) edge 
    diffractedrays(1,theta)=20*log10(norm(Ei*amplitude*exp(-1i*k*s)*diffcoef));%calculation of the field at 3 units diatance from edge  
    diff(1,theta)=Ei*amplitude*exp(-1i*k*s)*diffcoef;
    
    
    total(1,theta)=20*log10(norm(diff(1,theta)+dir(1,theta)+ref(1,theta)));
end

%%
%graphs ploting
figure(1)
plot(directrays);
title('Direct Incident field');
xlabel('angle(0-360)');
ylabel('field energy(dB)');
figure(2)
plot(reflectedrays);
title('Reflected field');
xlabel('angle(0-360)');
ylabel('field energy(dB)');
figure(3)
plot(diffractedrays);
title('Diffracted field');
xlabel('angle(0-360)');
ylabel('field energy(dB)');
figure(4)
plot(total);
title('Total field');
xlabel('angle(0-360)');
ylabel('field energy(dB)');

figure(5)
subplot(2,2,1);
plot(directrays);
title('Direct Incident field');
xlabel('angle(0-360)');
ylabel('field energy(dB)');
subplot(2,2,2);
plot(reflectedrays);
title('Reflected field');
xlabel('angle(0-360)');
ylabel('field energy(dB)');
subplot(2,2,3);
plot(diffractedrays);
title('Diffracted field');
xlabel('angle(0-360)');
ylabel('field energy(dB)');
subplot(2,2,4);
plot(total);
title('Total field');
xlabel('angle(0-360)');
ylabel('field energy(dB)');

figure(6);
plot(directrays);
hold on;
plot(reflectedrays);
hold on;
plot(diffractedrays);
hold on;
plot(total);
hold off;
legend('incident field','reflectedrays field','diffractrays field','totalfield');
xlabel('angle(0-360)');
ylabel('field energy');