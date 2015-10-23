%%FOR LOCAL ATTITUDE STABILITY, WE CONSIDER LINEARIZED MODEL


longstr1='1 25635U 99008B   12341.91858393  .00000229  00000-0  63384-4 0  6835';                           %Two Line Element Data for Satellite
longstr2='2 25635 096.4809 248.3902 0143174 115.1607 246.4523 14.46530052727099 ';
time='22-June-2014'

%STATE VARIABLES CORRESPONDING TO SMALL PERTURBATION IN ANGULAR VELOCITY
%AND ATTITUDE

dwx_0=0;                          %small perturbation in x component of angular velocity of Control CS wrt Orbital CS
dwy_0=0;                          %small perturbation in y component
dwz_0=0;                          %small perturbation in z component
dq1_0=0;                          %small perturbation in quaternion vector component along 'i' direction 
dq2_0=0;                          %along 'j' direction
dq3_0=-1;                          %along 'k' direction

Jx=181.78;                           %moment of inertia about principal x axis in kgm2
Jy=181.25;                           %moment of inertia about principal y axis in kgm2
Jz=1.28;                             %moment of inertia about principal z axis in kgm2

sigx=(Jy-Jz)/Jx;                     %few constants which are late required in computation of state matrix
sigy=(Jz-Jx)/Jy;
sigz=(Jx-Jy)/Jz;

%%STATE MODEL FOR LINEAR SYSTEM

% A=[0,0,0,(-2*k*sigx),0,0;
%     0,0,(wo*sigy),0,(2*k*sigy),0;
%     0,(wo*sigz),0,0,0,0;
%     0.5,0,0,0,0,0;
%     0,0.5,0,0,0,wo;
%     0,0,0.5,0,-wo,0];                       %6*6 state matrix

I=[Jx,0,0;0,Jy,0;0,0,Jz];                 %moment of inertia matrix resolved along principle body axes



wo=0.001047;              %orbital rate of satellite in rad/s(assuming that orbit is circular)

A22_init=0.001;                           %initial anlge of Orbital CS wrt World CS

% q1=0.01;                                  %breaking quaternion into its component, q4 is scalar
% q2=0;                                  %supposing that this quaternion contains information of orientation of Control CS wrt Orbital CS
% q3=0;
% q4=0.011;
 
% q1=0.0850;                                  %breaking quaternion into its component, q4 is scalar
% q2=0.0396;                                  %supposing that this quaternion contains information of orientation of Control CS wrt Orbital CS
% q3=-0.0472;                                 %quaternion initial value corresponding to 5 degrees(heading), -5(attitude), 10(banking)
% q4=0.9944;

q1=dq1_0;                                  %breaking quaternion into its component, q4 is scalar
q2=dq2_0;                                  %supposing that this quaternion contains information of orientation of Control CS wrt Orbital CS
q3=dq3_0;
q4=0;                                  %quaternion value corresponding to 40,-40,80


% Q1=1;                                       %supposing that this quaternion contains information of orientation of CS wrt Orbital CS
% Q2=0;
% Q3=0;
% Q4=0;

hg=100000000;                                %gain for angular velocity feedback control law
eg=300000;                                   %gain for quaternion control
k=3*wo*wo;                                   %this value is derived from theory of linearization

%%%feeding the dynamic model of satellite inside the simulation, so that it
%%%would generate response. We need to have code for system of differential
%%%equation

clear Q;                      %these are arrays used to plot graph, clear is used to clear their previous data
clear R;
clear S;
clear Tt;
clear W;

satrec = twoline2rvMOD(longstr1,longstr2);
tsince_offset=810415.5;                                       %tsince_offset is the number of mins elapsed since epoch time

for t0=1:1:200                                          %for loop for continuous evaluation of system for one complete orbital motion(time=100*60 seconds)

    
 tsince=tsince_offset+(t0/60);                                
[satrec, r, v] = sgp4(satrec,tsince);                       %sgp4 will give position of satellite
% 
r;
% Q(t)=r(1);
% R(t)=r(2);
% S(t)=r(3);
                                                            %converting r
                                                            %into metres
ecf=[-r(1)*1000,-r(2)*1000,-r(3)*1000];                     %sign changed of ecf vector after comparison with online available data

llh=ecf2llhT(ecf);

lat=llh(1)*(180/3.14);                                        %latitude and longitude in degrees
long=llh(2)*(180/3.14);
alt1=llh(3)/1000;                                             %altitude in kms

% Q(t0)=lat;
% R(t0)=long;
% S(t0)=alt1;

[Bx,By,Bz]=igrf(time,lat,long,alt1);                     %latitude always in [-90 to 90] range, longitude in [0 to 360] range
B_temp=[Bx;By;Bz];                                       %magnetic field in ECEF frame
[B_ECI] = ECEFtoECI(2456830.9851,B_temp,0,0);            %magnetic field in ECI frame

[B_ECI_new]=[B_ECI(1)*0.000000001;B_ECI(2)*0.000000001;B_ECI(3)*0.000000001];            %magnetic field in nano Tesla
% Q(t0)=B_ECI(1);
% R(t0)=B_ECI(2);
% S(t0)=B_ECI(3);

    distance=7469.76*t0;                          %since we want distance travelled by satellite in metres, so wo is in m/s
    A22=A22_init+wo*t0;                           %this is the angle(in radians) at which orbital CS is inclined wrt World CS(about x axis)
    
   

     Aoc=[((q1*q1)-(q2*q2)-(q3*q3)+(q4*q4)),2*((q1*q2)-(q3*q4)),2*((q1*q3)+(q2*q4));
        2*((q1*q2)+(q3*q4)),(-(q1*q1)+(q2*q2)-(q3*q3)+(q4*q4)),2*((q2*q3)-(q1*q4));
        2*((q1*q3)-(q2*q4)),2*((q2*q3)+(q1*q4)),(-(q1*q1)-(q2*q2)+(q3*q3)+(q4*q4))];          %quaternion matrix for converting B from orbital CS to control CS
    

%      Aoc=[((q1*q1)-(q2*q2)-(q3*q3)+(q4*q4)),2*((q1*q2)+(q3*q4)),2*((q1*q3)-(q2*q4));
%         2*((q1*q2)-(q3*q4)),(-(q1*q1)+(q2*q2)-(q3*q3)+(q4*q4)),2*((q2*q3)+(q1*q4));
%         2*((q1*q3)+(q2*q4)),2*((q2*q3)-(q1*q4)),(-(q1*q1)-(q2*q2)+(q3*q3)+(q4*q4))];          %quaternion matrix for converting B from orbital CS to control CS
%    
    
    Awo=[1,0,0;0,-sin(A22),cos(A22);0,cos(A22),sin(A22)];                           %rotation matrix to trasfer magnetic field from World CS(ECI) to Orbital CS
    A22_init=A22;
    

    
    B_t=(Awo*B_ECI_new);                               %magnetic field in Orbital CS

    Bx=B_t(1)*0.000000001;
    By=B_t(2)*0.000000001;
    Bz=B_t(3)*0.000000001;
    

dwx=dwx_0;
dwy=dwy_0;
dwz=dwz_0;
dq1=dq1_0;
dq2=dq2_0;
dq3=dq3_0;


h=1;                                         %step size for RK2 method
%for i=0:1:t0                                %RK2 method to solve system of dynamical equations 
    %t=t+(i*h);
		k1=(-2*k*sigx*dq1)-((Bz/Jx)*((hg*dwy)+(eg*dq2)))+((By/Jx)*((hg*dwz)+(eg*dq3)))
		k3=(wo*sigy*dwy)+(2*k*sigy*dq2)+((Bz/Jy)*((hg*dwx)+(eg*dq1)))-((Bx/Jy)*((hg*dwz)+(eg*dq2)));
		k5=(wo*sigz*dwy)-((By/Jz)*((hg*dwx)+(eg*dq1)))+((Bx/Jz)*((hg*dwy)+(eg*dq2)));
        k7=(0.5*dwx)-(wo*0.5);
        k9=(0.5*dwy)+(wo*dq3);
        k11=(0.5*dwz)-(wo*dq2);
		
		k2=(-2*k*sigx*(dq1+(k7*h))-((Bz/Jx)*((hg*(dwy+(k3*h)))+(eg*(dq2+(k9*h)))))+((By/Jx)*((hg*(dwz+(k5*h)))+(eg*(dq3+(k11*h))))))
		k4=((wo*sigy*(dwy+(k3*h)))+(2*k*sigy*(dq2+(k9*h)))+((Bz/Jy)*((hg*(dwx+(k1*h)))+(eg*(dq1+(k7*h)))))-((Bx/Jy)*((hg*(dwz+(k5*h)))+(eg*(dq2+(k9*h))))));
		k6=(wo*sigz*(dwy+(k3*h)))-((By/Jz)*((hg*(dwx+(k1*h)))+(eg*(dq1+(k7*h)))))+((Bx/Jz)*((hg*(dwy+(k3*h)))+(eg*(dq2+(k9*h)))));
		k8=(0.5*(dwx+(k1*h)))-(wo*0.5);
        k10=(0.5*(dwy+(k3*h)))+(wo*(dq3+(k11*h)));
        k12=(0.5*(dwz+(k5*h)))-(wo*(dq2+(k9*h)));
        

%         k1=(-2*k*sigx*dq1)-((Bz/Jx)*((hg*dwy)-(eg*dq2)))+((By/Jx)*((hg*dwz)-(eg*dq3)));
% 		k3=(wo*sigy*dwy)+(2*k*sigy*dq2)+((Bz/Jy)*((hg*dwx)-(eg*dq1)))-((Bx/Jy)*((hg*dwz)-(eg*dq2)));
% 		k5=(wo*sigz*dwy)-((By/Jz)*((hg*dwx)-(eg*dq1)))+((Bx/Jz)*((hg*dwy)-(eg*dq2)));
%         k7=(0.5*dwx)-(wo*0.5);
%         k9=(0.5*dwy)+(wo*dq3);
%         k11=(0.5*dwz)-(wo*dq2);
% 		
% 		k2=(-2*k*sigx*(dq1+(k7*h))-((Bz/Jx)*((hg*(dwy+(k3*h)))-(eg*(dq2+(k9*h)))))+((By/Jx)*((hg*(dwz+(k5*h)))-(eg*(dq3+(k11*h))))));
% 		k4=((wo*sigy*(dwy+(k3*h)))+(2*k*sigy*(dq2+(k9*h)))+((Bz/Jy)*((hg*(dwx+(k1*h)))-(eg*(dq1+(k7*h)))))-((Bx/Jy)*((hg*(dwz+(k5*h)))-(eg*(dq2+(k9*h))))));
% 		k6=(wo*sigz*(dwy+(k3*h)))-((By/Jz)*((hg*(dwx+(k1*h)))-(eg*(dq1+(k7*h)))))+((Bx/Jz)*((hg*(dwy+(k3*h)))-(eg*(dq2+(k9*h)))));
% 		k8=(0.5*(dwx+(k1*h)))-(wo*0.5);
%         k10=(0.5*(dwy+(k3*h)))+(wo*(dq3+(k11*h)));
%         k12=(0.5*(dwz+(k5*h)))-(wo*(dq2+(k9*h)));
        
        
        
        dwx=dwx+((k1+k2)*0.5*h);
		dwy=dwy+((k3+k4)*0.5*h);
		dwz=dwz+((k5+k6)*0.5*h);
        dq1=dq1+((k7+k8)*0.5*h);
        dq2=dq2+((k9+k10)*0.5*h);
        dq3=dq3+((k11+k12)*0.5*h);
        
        
        k13=-0.5*(((dwx-wo)*dq1)+((dwy+(2*dq3*wo))*dq2)+((dwz-(2*dq2*wo))*dq3));             %for calculating q4
        k14=-0.5*((((dwx+(k1*h))-wo)*(dq1+(k7*h)))+(((dwy+(k3*h))+(2*(dq3+(k11*h))*wo))*(dq2+(k9*h)))+(((dwz+(k5*h))-(2*(dq2+(k9*h))*wo))*(dq3+(k11*h))));
        
        q4=q4+((k13+k14)*0.5*h);
%end    

q_temp=[dq1,dq2,dq3,q4]

q1=q_temp(1);
q2=q_temp(2);
q3=q_temp(3);

if(q4<0)
    q1=-q1;
    q2=-q2;
    q3=-q3;
    q4=-q4;
end

q=[q1,q2,q3,q4]

%%%CALCULATION OF EULER ANGLES
   
% OUTPUT=SpinCalc('QtoEA321',q,0.01,0);
% 
% Q(t0)=OUTPUT(1)-180;
% R(t0)=OUTPUT(2)-180;
% S(t0)=OUTPUT(3)-180;

m11=(2*(q4*q4))+(2*(q1*q1))-1;
m12=(2*q1*q2)+(2*q4*q3);
m13=2*((q1*q3)-(q4*q2));
m23=((2*q3*q2)+(2*q4*q1));
m33=((2*(q4*q4))+(2*(q3*q3))-1);

  

Sa=((atan(m12/m11))*(180/3.14));
Q(t0)=Sa;
The=((asin(-m13))*(180/3.14));
R(t0)=The;
Phi=((atan2(m23,m33))*(180/3.14));
S(t0)=Phi;

% qx=q1;
% qy=q2;
% qz=q3;
% qw=q4;
% 
% heading = (atan2(2*qy*qw-2*qx*qz , 1 - 2*qy*qy - 2*qz*qz))*(180/3.14)
% Q(t0)=heading;
% attitude = (asin(2*qx*qy + 2*qz*qw))*(180/3.14)
% R(t0)=attitude;
% bank = (atan2(2*qx*qw-2*qy*qz , 1 - 2*qx*qx - 2*qz*qz))*(180/3.14)
% S(t0)=bank;


Ww=(Aoc*[wo;0;0])+[dwx;dwy;dwz];
% Wo=Ww-(Aoc*[wo;0;0]);                    %to get orbital from world
dwx_0=dwx;
dwy_0=dwy;
dwz_0=dwz;
dq1_0=dq1;
dq2_0=dq2;
dq3_0=dq3;
 
W(t0)=Ww(1);                               %array W for plotting graph




end


plot(Q);
hold on;
plot(R);
hold on;
plot(S);
hold on;
% plot(Tt);
% hold on;

% plot(W)
% hold on;



