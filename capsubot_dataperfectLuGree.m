
clc;
clear all;

a=8;
b=25;
x =-a:1:a ;
y=sqrt((1-x.^2/a.^2)*b.^2);
figure(10),plot(x+8,y), grid on, title('desired trajectory of capsubot'),xlabel('time(s)'),ylabel('capsubot displacement(cm)')
ylim([0 28])
j=30


for l=2:(2*a+1)
   
d_desired(l-1) = y(l)-y(l-1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%Abas optimisation%%%%%%%%%%%%%%%%%
% 
% %%%%%%%utorque plofile scenario 1%%%%%%%
% y11 = [7.5 7.5 2 2 0.5 0.5 8 8];
% x = [0 3 3 5 5 6.5 6.5 9.5];
% figure(1),plot(x,y11, 'LineWidth', 2, 'MarkerSize', 20);
% ylim([0,8.5]);
% xlim([0,10]);
% ylabel('accelaration'),title('utorque plofile scenario 1');
% xlabel('time');
% grid on;
% yticks([0]);
% xticks([0]);
% % %%%%%%%utorque plofile scenario 2 %%%%%%%
% y12 = [0.5 0.5 6 6 8 8 0.2 0.2];
% x = [0 3 3 5 5 6.5 6.5 9];
% figure(2),plot(x,y12,'LineWidth', 2, 'MarkerSize', 20);
% ylim([0,8.5]);
% xlim([0,9.5]);
% ylabel('accelaration'),title('utorque plofile scenario 2');
% xlabel('time');
% grid on;
% yticks([0]);
% xticks([0]);
% % 
% % %%%%%%%contrarium plofile scenario 1 %%%%%%%
% y2=[0.2 0.2 4 4 ];
% x2=[0 2 2 6];
% figure(3),plot(x2,y2,'LineWidth', 2, 'MarkerSize', 20),title('contrarium plofile scenario 1');
% ylim([0,4.5]);
% xlim([0,7]);
% ylabel('accelaration');
% xlabel('time');
% grid on;
% yticks([0]);
% xticks([0]);
% % 
% % %%%%%%%contrarium plofile scenario 2 %%%%%%%
% y22=[4 4 0.2 0.2 ];
% x22=[0 2 2 6];
% figure(4),plot(x22,y22,'LineWidth', 2, 'MarkerSize', 20),title('contrarium plofile scenario 2');
% ylim([0,4.5]);
% xlim([0,7]);
% ylabel('accelaration');
% xlabel('time');
% grid on;
% yticks([0]);
% xticks([0]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 0
% %%%%%%%%%%% Desired capsubot position %%%%%%%%%%%%%%%%%
% %%%%IM accelaration utorque profile %%%%%%
% yIM=[6 6  0.2 0.2 6 6 6];
% xIM=[0 3 3 5 5 8 9];
% figure(5),plot(xIM,yIM,'LineWidth', 2, 'MarkerSize', 20');
% ylim([0,7])
% xlim([0,9.5])
% ylabel('accelaration');
% xlabel('time')
% yticks([0]);
% xticks([0]);
% grid on;
% hold on;
% % %%%% capsubot accelaration utorque profile %%%%%%
% ycapsubot=[5 5 6.5 6.5 4 4 5 5 5];
% xcapsubot=[0 3 3 5 5 7 7 8 9 ];
% plot(xcapsubot,ycapsubot,'LineWidth', 2, 'MarkerSize', 20'),title('IM and capsubot accelaration for utorque profile');
% legend({'IM accelaration','capsubot accelaration ,'},'Location','southeast');
% yticks([0])
% xticks([0])
% grid off 
% 
%%%IM Velocity optimisation utorque profile %%%%%%
% yIMve=[0 1 -1 0];
% xIMve=[0 4 8 16];
% figure(6),plot(xIMve,yIMve,'LineWidth', 2, 'MarkerSize', 20');
% yticks([0]);
% xticks([0]);
% xlim([0 13])
% ylim([-1.2 1.2])
% hold on
% %%%%acpsubot Velocity optimisation utorque profile %%%%%%
% ycabuve=[0 0 0.5 0 0];
% xcabuve=[0 4 8 14 16];
% plot(xcabuve,ycabuve,'LineWidth', 2, 'MarkerSize', 20');
% yticks([0]);
% xticks([0]);
% xlim([0 17])
% ylim([-1.2 1.2])
% legend({'IM velocity','Capsubot velocity ,'},'Location','northeast');
% title('IM and Caspsubot Velocity for Utorque profile')
% ylabel('Velocity'),xlabel('time')
% hold off

%%IM position optimisation contrarium profile %%%%%%
% yIMa=[2 2 0 0 0.8 0.8];
% xIMa=[0 2 2 6 6 14];
% figure(6),plot(xIMa,yIMa,'LineWidth', 2, 'MarkerSize', 20');
% yticks([0]);
% xticks([0]);
% hold on
% %%%%acpsubot accelaration optimisation contrarium profile %%%%%%
% ycapa=[-2 -2 1.8 1.8];
% xcapa=[0 2 2 14];
% figure(6),plot(xcapa,ycapa,'LineWidth', 2, 'MarkerSize', 20');
% yticks([0]);
% xticks([0]);
% xlim([0 14.2])
% ylim([-2.2 2.2])
% legend({'IM accelaration','Capsubot accelaration ,'},'Location','southeast');
% title('IM and Caspsubot accelaration for contrarium profile')
% ylabel('accelaration'),xlabel('time')
% hold off

%%%IM velocity optimisation contrarium profile %%%%%%
% yIMvc=[0 1 0 0 0];
% xIMvc=[0 3 7 12 20];
% figure(6),plot(xIMvc,yIMvc,'LineWidth', 2, 'MarkerSize', 20');
% yticks([0]);
% xticks([0]);
% xlim([0 20.5])
% ylim([-1.7 1.2])
% hold on
% %%%%acpsubot velocity optimisation contrarium profile %%%%%%
% ycapvc=[0 -1.5 0 ];
% xcapvc=[0 3 20];
% plot(xcapvc,ycapvc,'LineWidth', 2, 'MarkerSize', 20');
% yticks([0]);
% xticks([0]);
% legend({'Capsubot Velocity','IM Velocity ,'},'Location','northeast');
% title('IM and Caspsubot Velocity for contrarium profile')
% ylabel('Velocity'),xlabel('time')
% hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%Abas%%%%%%%%%%%%%%%%%
constrains_data= [.396 0.05 .1 0.2 .009 0.5 0.001 0.3 0.1 1 0.3];

a_last=-30;
a_first=-8;
a_min=5;

data=constrains_data;

%%%%%%%%%%%%%Abas%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


M=data(1);             %Mass of the capsubot in Kg
m=data(2);             % Mass of the Inner Mass in Kg
m1=data(2);            % Mass of the Inner Mass in Kg
mu=data(3);            %friction coefficient of capsule
mu_1=data(4);          %friction coefficient of IM
mu_m=data(4);          %friction coefficient of IM
a_stoke=data(5);
x_6=data(6);
z_cap=data(7);
sig_IM0=data(8);
sig_IM1=data(9);
sig_IM2=data(9)*2;
sig_cap0=data(8);
sig_cap1=data(11);
f_sm=data(8);
fv_IM=data(9);
f_sM=0.49;
fv_Cap=data(9);
z_IM=data(9);



k=a_stoke;
close all;
a1=-a_stoke; % Integrator in the First Subsystem (Signal Generator)


g=9.8;         % Acceleration due to Gravity = 9.8 m/sec^2
k1=300;        %aero controller coefficient1   
k2=100;        %aero controller coefficient1


n=(abs(a_last)- abs(a_first))/.5 + 1;       %% n=45
a_all1=a_first:-.5:a_last;                  %%create a vector from -8:0.5:-30 
a_all2=a_min*ones(1,n);                     %% create column of 5 at n=45                         
a_c_desired=[a_all1' a_all2'];              %% desired accelaration 
a_c_desired_forward=[a_all1' a_all2'];      %%desired forward accelaration
a_c_desired_backward=[-a_all1' -a_all2'];   %%desired backward accelaration 

a_u_desired=[a_all2' a_all1' a_all1' a_all2'];
a_u_desired_forward=[a_all2' a_all1' a_all1' a_all2'];
a_u_desired_backward=[-a_all2' -a_all1' -a_all1' -a_all2'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Contrarium Cycle %%%%%%%%%%%%%%%%%%


tc = zeros(n,2);
d_c=zeros(n,1);
for i=1:n
        a_c = a_c_desired_forward(i,:);  
         
        a_c1= a_c(1);
        a_c2= a_c(2);


       aMc1=(-m*a_c1-mu*M*g)/M;                          %%equation 3.15 for i1
       aMc2=(-m*a_c2-mu*M*g)/M;                          %%equation 3.15 for i2
 

        Vmc=-sqrt((-4*k*(a_c1)^2*a_c2*aMc2)/(a_c1*aMc2*(a_c2-a_c1)- aMc1*a_c2*(aMc2- aMc1)));  %%equation 3.14

        tc1=Vmc/a_c1;                                    %%equation 3.11 for tc1                                    
        tc2=tc1-Vmc/a_c2;                                %%equation 3.11 for tc2

         
       VMc=Vmc* aMc1/a_c1;                               %%equation 3.13                                                   
        d_c(i)= VMc^2/(2* aMc1)- VMc^2/(2*aMc2);
     
        tc(i,:)=[tc1 tc2];
      
        
end




%%%%%%%%% Contrarium Cycle %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%% Utroque Cycle %%%%%%%%%%%%%%%%%%
% 
tu =zeros(n,4);                         %%utorque time
d_u=zeros(n,1);
d_u_individual=zeros(n,3);
for i=1:n
        a_u = a_u_desired_forward(i,:);  
         
        a_u1= a_u(1);
        a_u2= a_u(2);

        a_u3= a_u(3);
        a_u4= a_u(4);

        a_u2p=(-m*a_u2-mu*M*g)/M;
        
 a_u3p=(-m*a_u3-mu*M*g)/M;
 a_u4p=(-m*a_u4-mu*M*g)/M;
 
 
 v_u12=sqrt(4*k*a_u1*a_u2^2/(a_u2^2-a_u1*a_u2-a_u1*a_u2p));

 vp_u12=- v_u12*a_u2p/a_u2;
 
 

A=1/a_u3 - 1/a_u4 + (1/a_u4p-1/a_u3p)*a_u3p^2/a_u3^2;
B=2*vp_u12*(a_u3p/a_u3)*(1/a_u4p-1/a_u3p);
C=4*k+vp_u12^2/a_u4p;

 
 v_u34_approx1=(-B+sqrt(B^2-4*A*C))/(2*A) ;

 v_u34_approx2=(-B-sqrt(B^2-4*A*C))/(2*A) ;
 
 v_u34=v_u34_approx1;


tu1=v_u12/a_u1;
tu2=tu1-v_u12/a_u2;

tu3=tu2+v_u34/a_u3;
tu4=tu3-v_u34/a_u4;


v_u12p=-v_u12*a_u2p/a_u2;
v_u34p=v_u34*a_u3p/a_u3 + v_u12p;

d_u(i)=0 + v_u12p^2/(2*a_u2p) + (v_u34p^2 - v_u12p^2)/(2*a_u3p) - (v_u34p^2)/(2*a_u4p);
d_u_individual(i,:)=  [v_u12p^2/(2*a_u2p),  (v_u34p^2 - v_u12p^2)/(2*a_u3p),  - (v_u34p^2)/(2*a_u4p)];

tu(i,:)=[tu1 tu2 tu3 tu4];

end

d_u_b=-d_u;


%%%%%%%%% Utroque Cycle %%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v_c_average=d_c./tc(:,2) *1000;
v_u_average =d_u./tu(:,4)*1000;


v_c_average_forward = v_c_average ;
v_c_average_backward = - v_c_average ;

v_u_average_forward = v_u_average ;
v_u_average_backward = - v_u_average ;

d_desired ;

d_desired_new (1) = d_desired(1);
t_plus (1) = 0;
t_plus_t = 0;
y_minus (1) = 0; 
y_minus_t = 0; 

for i=1:8
   
            [Y , I]=min(abs(d_desired_new(i)/(1 + t_plus_t) - v_u_average_forward));
            v_u_average_forward(I);
            tu(I,4);
            floor((1 + t_plus (i))/tu(I,4)); %%
            profile(i) = I;
            % cycle(i)= floor(1 /tu(I,4)); %% Calculation of number of cycles
            cycle(i)= floor((1 + t_plus_t)/tu(I,4)); %% Calculation of number of cycles

            d_desired_new(i+1)= y(i+2) - d_u(I)*1000 * floor((1 + t_plus_t)/tu(I,4))- y_minus_t;
            % t_plus_correct (i)= ;

            y_minus (i+1) = d_u(I)*1000 * floor((1 + t_plus_t)/tu(I,4));
            y_minus_t = y_minus_t + y_minus (i+1);

            t_plus (i+1) = 1 - tu(I,4) * floor((1 + t_plus_t)/tu(I,4));
            t_plus_t = t_plus_t + t_plus (i+1);
            t_plus_t_view(i)=t_plus_t;


end

% find ()                  
time_sim_individual_cum_store = 0; 
%for p=1:(a-1)
for p=1:a
time_sim_individual(p)= tu(profile(p),4)*cycle(p);
time_sim_individual_cum (p) = time_sim_individual_cum_store + time_sim_individual(p);
time_sim_individual_cum_store = time_sim_individual_cum (p);
end



for i=9
     
            [Y , I] = min(abs(d_desired_new(i)/(1 + t_plus_t) - v_c_average_forward));
           
            tc(I,2);
            floor((1 + t_plus (i))/tu(I,4)); %%
            profile(i) = I;
            
            %cycle_c(i)= floor((1 + t_plus_t)/tu(I,4)); %% Calculation of number of cycles

            d_desired_new(i+1)= y(i+2) - d_c(I)*1000 - y_minus_t;
           

            y_minus (i+1) = d_c(I)*1000 ;
            y_minus_t = y_minus_t + y_minus (i+1);

            t_plus (i+1) = 1  - tc(I,2) ;
            t_plus_t = t_plus_t + t_plus (i+1);
            t_plus_t_view(i)=t_plus_t;

end

% find ()                  
%for p=1:(a-1)
for p=9
% cycle(p);
time_sim_individual(p)= tc(profile(p),2);
time_sim_individual_cum (p) = time_sim_individual_cum_store + time_sim_individual(p);
time_sim_individual_cum_store = time_sim_individual_cum (p);
end

y(18)=0;
for i=10:16
      %if d_desired_new(i)> 0
            [Y , I]=min(abs(d_desired_new(i)/(1 + t_plus_t) - v_u_average_backward));
            v_u_average_forward(I);
            tu(I,4);
            floor((1 + t_plus (i))/tu(I,4)); %%
            profile(i) = I;
            % cycle(i)= floor(1 /tu(I,4)); %% Calculation of number of cycles
            cycle(i)= floor((1 + t_plus_t)/tu(I,4)); %% Calculation of number of cycles

            d_desired_new(i+1)= y(i+2) - d_u_b(I)*1000 * floor((1 + t_plus_t)/tu(I,4)) -  y_minus_t ;
            % t_plus_correct (i)= ;

            y_minus (i+1) = d_u_b(I)*1000 * floor((1 + t_plus_t)/tu(I,4));
            y_minus_t = y_minus_t + y_minus (i+1);

            t_plus (i+1) = 1  - tu(I,4) * floor((1 + t_plus_t)/tu(I,4));
            t_plus_t = t_plus_t + t_plus (i+1);
            t_plus_t_view(i)=t_plus_t;

     % elseif d_desired_new(i)<0
          
      %end 


%t_min (i+1)
% t_plus_correct (i)=1  - floor((1 + t_plus)/tu(I,4));
end

% find ()                  

%for p=1:(a-1)
for p=10:16
% cycle(p);
time_sim_individual(p)= tu(profile(p),4)*cycle(p);
time_sim_individual_cum (p) = time_sim_individual_cum_store + time_sim_individual(p);
time_sim_individual_cum_store = time_sim_individual_cum (p);
end


% 
% 
time=28;
if time<= time_sim_individual_cum (1)
    t_add1=0;
for p=1:cycle(1)
   for q=0:001:49
            if q<=(t_add1 + tu(profile(1),1))
                    a_s =a_u_desired(profile(1),1) ;
            elseif (q>(t_add1 + tu(profile(1),1))) && (q<=(t_add1 + tu(profile(1),2))) 
                    a_s =a_u_desired(profile(1),2) ;
            elseif (q>(t_add1 + tu(profile(1),2))) && (q<=(t_add1 + tu(profile(1),3))) 
                    a_s =a_u_desired(profile(1),3) ;
            elseif (q>(t_add1 + tu(profile(1),3))) && (q<=(t_add1 + tu(profile(1),4)))
                    a_s = a_u_desired(profile(1),4) ;
%             else
%                     a_s = a_u_desired(profile(1),4) ;
            end
            
           t_add1 =tu(profile(1),4) ; 
    end
end 

end
%%%%%
plot(simoutLuGreeaM,'b')
hold on 
plot(simoutcoulombsa,'r')
ad = accelarationerror(:,1)+3
title('Capsubot accelearation with LuGree and Coulombs friction model')
xlabel('time'),ylabel('accelaration')
legend({'Capsubot accelaration with LuGree friction','Capsubot accelaration with Coulombs friction ,'},'Location','northeast')
plot(ad,'black')
hold off

%%%

plot(simoutLuGreeVc,'b')
hold on 
plot(simoutcoulombsVc,'r')
aa = velocityerror(:,1)+0.04
plot(aa,'black')
legend({'Capsubot Velocity with LuGree friction','Capsubot Velocity with Coulombs friction ','error'},'Location','northeast')
xlabel('time'),ylabel('Velocity')
title('Capsubot Velocity with LuGree and Coulombs friction model')
hold off

%%%

plot(simoutLuGreeD,'b')
hold on 
plot(simoutcoulombsD,'r')
title('Capsubot displacement with LuGree and Coulombs friction model')
xlabel('time'),ylabel('displacement')
legend({'Capsubot displacement with LuGree friction','Capsubot displacemment with Coulombs friction ,'},'Location','northeast')
hold off
