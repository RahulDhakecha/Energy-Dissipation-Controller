longstr1='1 25635U 99008B   12341.91858393  .00000229  00000-0  63384-4 0  6835';                           %Two Line Element Data for Satellite
longstr2='2 25635 096.4809 248.3902 0143174 115.1607 246.4523 14.46530052727099 ';
time='22-June-2014'
wo=0.001047;              %orbital rate of satellite in rad/s(assuming that orbit is circular)
alt=755000;              %altitude of satellite form point on earth in metres
rad=6378100;             %radius of earth in metres
flag=0;

prev_long=0;             %used for longitude measurement
prev_time=0;             %used for longitude measurement

A22_init=0.001;                           %initial anlge of Orbital CS wrt World CS
 
% q1=0.0850;                                  %breaking quaternion into its component, q4 is scalar
% q2=0.0396;                                  %supposing that this quaternion contains information of orientation of Control CS wrt Orbital CS
% q3=-0.0472;                                 %quaternion initial value corresponding to 5 degrees(heading), -5(attitude), 10(banking)
% q4=0.9944;

q1=0.6570;
q2=0.0398;
q3=-0.4526;
q4=0.6016;

Wo0=[0.005;0.003;-0.003];            %initial angular velocity of Control CS with respect to Orbit CS(in rad/s)

%h=5;                                %gain for feedback control law
Jx=181.78;                           %moment of inertia about principal x axis in kgm2
Jy=181.25;                           %moment of inertia about principal y axis in kgm2
Jz=1.28;                             %moment of inertia about principal z axis in kgm2

I=[Jx,0,0;0,Jy,0;0,0,Jz];            %moment of inertia matrix

% E_kin=0.5*(transpose(Wo0)*(I*Wo0));
% E_gg=(3/2)*wo*wo*(((transpose(cko))*(I*cko))-Jz);
% E_gyro=(1/2)*wo*wo*(Jx-((transpose(cio))*(I*cio)));
E_gg_max=(3/2)*wo*wo*(Jx-Jz);
E_gyro_max=(1/2)*wo*wo*(Jx-Jz);

%%%feeding the dynamic model of satellite inside the simulation, so that it
%%%would generate response. We need to have code for system of differential
%%%equation

clear Q;                      %these are arrays used to plot graph, clear is used to clear their previous data
clear R;
clear S;
clear Tt;
clear W;

satrec = twoline2rvMOD(longstr1,longstr2);
tsince_offset=810465.5;                                       %tsince_offset is the number of mins elapsed since epoch time

for t0=1:1:20                                          %for loop for continuous evaluation of system for one complete orbital motion(time=100*60 seconds)

  if(q4>0)                                     %whether this should be placed here or not?????
            flag=1;
            q1=-q1;
            q2=-q2;
            q3=-q3;
            q4=-q4;
        end    
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
B_temp=[Bx;By;Bz];
[B_ECI] = ECEFtoECI(2456830.9851,B_temp,0,0);  
[B_ECI_new]=[B_ECI(1)*0.000000001;B_ECI(2)*0.000000001;B_ECI(3)*0.000000001];            %magnetic field in nanotachnology
% Q(t0)=B_ECI(1);
% R(t0)=B_ECI(2);
% S(t0)=B_ECI(3);

    distance=7469.76*t0;                          %since we want distance travelled by satellite in metres, so wo is in m/s
    if(lat>80)
    A22_init=0;
    end

A22=A22_init+(wo*1);                        %angular rate multiplied by 1 sec because we are computing dynamics at every one second


    Aoc=[((q1*q1)-(q2*q2)-(q3*q3)+(q4*q4)),2*((q1*q2)-(q3*q4)),2*((q1*q3)+(q2*q4));
        2*((q1*q2)+(q3*q4)),(-(q1*q1)+(q2*q2)-(q3*q3)+(q4*q4)),2*((q2*q3)-(q1*q4));
        2*((q1*q3)-(q2*q4)),2*((q2*q3)+(q1*q4)),(-(q1*q1)-(q2*q2)+(q3*q3)+(q4*q4))];          %quaternion matrix for converting B from orbital CS to control CS
  

%      Aoc=[((q1*q1)-(q2*q2)-(q3*q3)+(q4*q4)),2*((q1*q2)+(q3*q4)),2*((q1*q3)-(q2*q4));
%         2*((q1*q2)-(q3*q4)),(-(q1*q1)+(q2*q2)-(q3*q3)+(q4*q4)),2*((q2*q3)+(q1*q4));
%         2*((q1*q3)+(q2*q4)),2*((q2*q3)-(q1*q4)),(-(q1*q1)-(q2*q2)+(q3*q3)+(q4*q4))];          %quaternion matrix for converting B from orbital CS to control CS
%    
    Awo=[-1,0,0;0,-cos(A22),-sin(A22);0,-sin(A22),cos(A22)];
   % Awo=[1,0,0;0,-sin(A22),cos(A22);0,cos(A22),sin(A22)];
    A22_init=A22;
    
%     Awo=[((Q1*Q1)-(Q2*Q2)-(Q3*Q3)+(Q4*Q4)),2*((Q1*Q2)-(Q3*Q4)),2*((Q1*Q3)+(Q2*Q4));
%         2*((Q1*Q2)+(Q3*Q4)),(-(Q1*Q1)+(Q2*Q2)-(Q3*Q3)+(Q4*Q4)),2*((Q2*Q3)-(Q1*Q4));
%         2*((Q1*Q3)-(Q2*Q4)),2*((Q2*Q3)+(Q1*Q4)),(-(Q1*Q1)-(Q2*Q2)+(Q3*Q3)+(Q4*Q4))];          %quaternion matrix for converting B from world CS to orbital CS
    
    %R=[q4,q3,-q2,q1;-q3,q4,q1,q2;-q2,-q1,q4,q3;-q1,-q2,-q3,q4];
    
    B=Aoc*(Awo*B_ECI_new);

cko=[2*((q1*q3)-(q2*q4));2*((q2*q3)+(q1*q4));(-(q1*q1)-(q2*q2)+(q3*q3)+(q4*q4))];              %unit vector along z direction in orbital CS projected along Control CS
Ngg=3*wo*wo*cross(cko,(I*cko));                              %torque due to gravity gradient
ckoz=transpose(cko)*[0;0;1];   
cio=[((q1*q1)-(q2*q2)-(q3*q3)+(q4*q4));2*((q1*q2)-(q3*q4));2*((q1*q3)+(q2*q4))];

Q_matrix=[0,-q3,q2;q3,0,-q1;-q2,q1,0];                            %assuming quaternion multiplication and ignoring scalar part,this matrix is used to represent quaternion multiplication with vector    

E_kin=0.5*(transpose(Wo0)*(I*Wo0))
E_gg=(3/2)*wo*wo*(((transpose(cko))*(I*cko))-Jz);
E_gyro=(1/2)*wo*wo*(Jx-((transpose(cio))*(I*cio)));

Q(t0)=E_kin;
R(t0)=E_gg;
S(t0)=E_gyro;

E_tot=E_kin+E_gg+E_gyro;
if(E_tot>(E_gg_max+E_gyro_max))
    m=100000000*(cross(Wo0,B));                                           %angular rate control law
else
    if(ckoz>0)
        m=((100000000*(cross(Wo0,B)))-(300000*(Q_matrix*B)));                                           %attitude control law
    else
        m=[0;0;0];
    end
end


T=cross(m,B);                                                %torque to be applied on each magnetic coil,this torque would be applied to body CS
Tx=T(1);
Ty=T(2);
Tz=T(3);


Ww0=Wo0+(Aoc*[wo;0;0]);                         %initial condition for angular velocity wrt to World, because dynamical eqns consider that velocity

Wx=Ww0(1);
Wy=Ww0(2);
Wz=Ww0(3);

%     if(flag==1)
%     flag=0;
%     q1=-q1;
%             q2=-q2;
%             q3=-q3;
%             q4=-q4;
% end


h=1;
%t=2;                       %time at which we want to evaluate the system
n=t0/h;

%for i=0:1:t0                                %RK2 method to solve system of dynamical equations 
    %t=t+(i*h);
		k1=((Tx+Ngg(1))-((Jz-Jy)*Wy*Wz))/Jx;
		k3=((Ty+Ngg(2))-((Jx-Jz)*Wx*Wz))/Jy;
		k5=((Tz+Ngg(3))-((Jy-Jx)*Wy*Wx))/Jz;
		
		k2=((Tx+Ngg(1))-((Jz-Jy)*(Wy+(k3*h))*(Wz+(k5*h))))/Jx;
		k4=((Ty+Ngg(2))-((Jx-Jz)*(Wx+(k1*h))*(Wz+(k5*h))))/Jy;
		k6=((Tz+Ngg(3))-((Jy-Jx)*(Wy+(k3*h))*(Wx+(k1*h))))/Jz;
		
		Wx=Wx+((k1+k2)*0.5*h);
		Wy=Wy+((k3+k4)*0.5*h);
		Wz=Wz+((k5+k6)*0.5*h);
        
        Ww=[Wx;Wy;Wz];
        
        Wo=Ww-(Aoc*[wo;0;0]);                                     %to transfer rate wrt world to rate wrt orbit
        
        Wx1=Wo(1);
        Wy1=Wo(2);
        Wz1=Wo(3);
        
        
        K1=0.5*((Wx1*q4)+(Wz1*q2)-(Wy1*q3));             % for q1
		K3=0.5*((Wy1*q4)-(Wz1*q1)+(Wx1*q3));            % for q2
		K5=0.5*((Wz1*q4)+(Wy1*q1)-(Wx1*q2));             % for q3
        K7=0.5*(-(Wx1*q1)-(Wy1*q2)-(Wz1*q3));            % for q4
		
		K2=0.5*((Wx1*(q4+(K7*h)))+(Wz1*(q2+(K3*h)))-(Wy1*(q3+(K5*h))));              % for q1
		K4=0.5*((Wy1*(q4+(K7*h)))-(Wz1*(q1+(K1*h)))+(Wx1*(q3+(K5*h))));              % for q2
		K6=0.5*((Wz1*(q4+(K7*h)))+(Wy1*(q1+(K1*h)))-(Wx1*(q2+(K3*h))));              % for q3
        K8=0.5*(-(Wx1*(q1+(K1*h)))-(Wy1*(q2+(K3*h)))-(Wz1*(q3+(K5*h))));             % for q4
		
		q1=q1+((K1+K2)*0.5*h);
		q2=q2+((K3+K4)*0.5*h);
		q3=q3+((K5+K6)*0.5*h);
        q4=q4+((K7+K8)*0.5*h);
        
      
       
        
%         if(q4<0)                                     %whether this should be placed here or not?????
%             q1=-q1;
%             q2=-q2;
%             q3=-q3;
%             q4=-q4;
%         end
%         
%         Q(t0)=q1;
%         R(t0)=q2;
%         S(t0)=q3;
%         Tt(t0)=q4;
        

        
        
%        Wx=Wx+wo;                                         % to change rate wrt Orbit to rate wrt World
%end


q=[q1,q2,q3,q4]
% Q(t0)=q(1);
% R(t0)=q(2);
% S(t0)=q(3);
% Tt(t0)=q(4);

cko=[2*((q1*q3)-(q2*q4));2*((q2*q3)+(q1*q4));(-(q1*q1)-(q2*q2)+(q3*q3)+(q4*q4))];
ckoz=transpose(cko)*[0;0;1]
% Tt(t0)=ckoz;
cio=[((q1*q1)-(q2*q2)-(q3*q3)+(q4*q4));2*((q1*q2)-(q3*q4));2*((q1*q3)+(q2*q4))];
ciox=transpose(cio)*[1;0;0];
% R(t0)=ciox;

%%%CALCULATION OF EULER ANGLES
   

m11=(2*(q4*q4))+(2*(q1*q1))-1;
m12=(2*q1*q2)+(2*q4*q3);
m13=2*((q1*q3)-(q4*q2));
m23=(2*q3*q2)+(2*q4*q1);
m33=(2*(q4*q4))+(2*(q3*q3))-1;

  

% Sa=(atan2(m12,m11)*(180/3.14));
% Q(t0)=Sa;
% The=(asin(-m13)*(180/3.14));
% R(t0)=The;
% Phi=(atan2(m23,m33)*(180/3.14));
% S(t0)=Phi;

% Ww=[Wx,Wy,Wz];
% Wo=Ww-(Aoc*[wo;0;0]);                    %to get orbital from world
Wo0=Wo
 
W(t0)=Wo(1);


v1=((q1*q1)-(q2*q2)-(q3*q3)+(q4*q4));
v2=2*((q1*q2)-(q3*q4));
v3=2*((q1*q3)+(q2*q4));

end

plot(Q);
hold on;
% plot(R);
% hold on;
% plot(S);
% hold on;
% plot(Tt);
% hold on;

% plot(W)
% hold on;


%%this system of differential equation will give us angular velocity with
%%respect to world CS. We need to convert this into orbital CS. This should
%%be our simulation output. This conversion can be carried by formula 2.20
%%(Pg 15)Probably it will fluctuate around 0. 

% Wxo=Wx-wo;
% Wyo=Wy;
% Wzo=Wz;

%%If we want to check for attitude, then simulation output should be
%%projection of control CS on orbital CS

%%What if we want to check for angular velocity??



