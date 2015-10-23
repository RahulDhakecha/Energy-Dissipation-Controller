wo=0.001047;              %orbital rate of satellite in rad/s(assuming that orbit is circular)
alt=755000;              %altitude of satellite form point on earth
rad=6378100;             %radius of earth in metres
time='17-April-1997';       %date on which satellite position is to be calculated


Wo0=[0.005,0.003,-0.003];            %initial angular velocity of Control CS with respect to Orbit CS

%h=5;                                %gain for feedback control law
Jx=181.78;                           %moment of inertia about principal x axis in kgm2
Jx=181.25;                           %moment of inertia about principal y axis in kgm2
Jx=1.28;                             %moment of inertia about principal z axis in kgm2



%%%feeding the dynamic model of satellite inside the simulation, so that it
%%%would generate response. We need to have code for system of differential
%%%equation


for t0=1:1:60                                           %for loop for continuous evaluation of system for one complete orbital motion(time=100*60 seconds)

    
    if(t0>6000)
        t1=floor(t0/6000);
        t=t0-(t1*6000);
    else
        t=t0;
    end
    distance=wo*t;                          %distance travelled by satellite


                        %LATITUDE MEASUREMENT


                        if(t<=1500)
                            lat=(distance/(rad+alt))*(180/3.14);          %latitude of point on earth from where altitude is measured
                        end
                        if(t>1500&&t<=3000)
                            lat=((22409294.56-distance)/(rad+alt))*(180/3.14);
                        end
                        if(t>3000&&t<=4500)
                            lat=((22409294.56-distance)/(rad+alt))*(180/3.14);
                        end
                        if(t>4500&&t<=6000)
                            lat=((distance-44818589.12)/(rad+alt))*(180/3.14);
                        end


                       %LONGITUDE MEASUREMENT


                         if(t0>86400)
                            t2=floor(t0/86400);
                            t3=t0-(t2*86400);
                         else
                            t3=t0;
                         end
                            long=(t3*360)/86400;
                            
                            
                           alt=alt/1000;                               %altitude converted in kms
    
    
    [Bx,By,Bz]=igrf(time,lat,long,alt);                                %magnetic field of earth obtained by IGRF function
    B=[Bx,By,Bz];


m=100000000*(cross(Wo0,B));                                           %magnetic moment
T=cross(m,B);                                                %torque to be applied on each magnetic coil,this torque would be applied to body CS
Tx=T(1);
Ty=T(2);
Tz=T(3);
Ww0=Wo0+[wo,0,0];                         %initial condition for angular velocity wrt to World, because dynamical eqns consider that velocity

inits='Wx(0)=Ww0(1),Wy(0)=Ww0(2),Wz(0)=Ww0(3)';
[Wx,Wy,Wz]=dsolve('Dx=(Tx-(Jz-Jy)*Wy*Wz)/Jx','Dy=(Ty-(Jx-Jz)*Wx*Wz)/Jy','Dz=(Tz-(Jy-Jx)*Wy*Wx)/Jz',inits)
Ww=[Wx,Wy,Wz];
Wo=Ww-[wo,0,0];
Wo0=Wo
 




% t = linspace(0,0.75,20);
% z = eval(vectorize(Wx));
% plot(t,z)



end


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



