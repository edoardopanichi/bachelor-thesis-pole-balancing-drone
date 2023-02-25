%% Variable of the system
close all;
clear all;
clc;
% Numerical parameters of the system
m = 1.55;
g = 9.81;
Jx = 0.029145975;
Jy = 0.029418104;
Jz = 0.055519004;
L = 1;
% At the equilibrium (Linearization)
U_star = [-m*g 0 0 0];
X_star = [0 0 0 0 10 0 0 0 0 0 0 0 0 0 0 0];

%% Linear System & Analysis of Controllability and Observability

% Matrix A of the whole system
A = [0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 g 0 0 0 0 0 0 0 0;
     0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 -g 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
     0 0 0 0 0 0 0 -g 0 0 0 0 g/L 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
     0 0 0 0 0 0 g 0 0 0 0 0 0 0 g/L 0];
 
 % Matrix A consedering only the states of the pendulum
 Ap = [0 1 0 0;
       g/L 0 0 0;
       0 0 0 1;
       0 0 g/L 0];
   
 % Matrix B of the whole system
 B = [0 0 0 0;
      0 0 0 0;
      0 0 0 0;
      0 0 0 0;
      0 0 0 0;
      -1/m 0 0 0;
      0 0 0 0;
      0 0 0 0;
      0 0 0 0;
      0 1/Jx 0 0;
      0 0 1/Jy 0;
      0 0 0 1/Jz;
      0 0 0 0;
      0 0 0 0;
      0 0 0 0;
      0 0 0 0];
  
% Matrix C for the whole system, considering all the outputs (used for 
% Frequency shaped control
C = eye(16); % I am interested in the whole set of states

% Matrix C of the whole system, considering as exits the position of 
% the drone, RPY and the pole
%   C = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%        0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
%        0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
%        0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
%        0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
%        0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
%        0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
%        0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0];
 
% Matrix A consedering only the states of the pendulum
Cp = [1 0 0 0;
      0 0 1 0];
%Controllability Matrix of the whole system 
CO = ctrb(A,B);
Rc = rank(CO);

% Observability Matrix of the system considering the whole system
%   OB = obsv(A,C);
%   Ro = rank(OB);

%Observability Matrix of the pendulum only 
OBp = obsv(Ap,Cp);
Rop = rank(OBp);

%I define the polynomial that contains the solutions I'd have as eigenvalues
%of the matrix A+HC
P = [-1 -5 -12 -2];
P1 = [-10 -50 -120 -20];

%https://it.mathworks.com/help/control/getstart/pole-placement.html
H = transpose(place(Ap',-Cp',P));
H1 = transpose(place(Ap',-Cp',P1));

eig(Ap+H*Cp)
eig(Ap+H1*Cp)


%% Frequency Shaped Control
%ATTENTION:
%ATTENTION:
%ATTENTION:
%ATTENTION: This sketch is equal to the Drone_Control_CASY file but without
%the filters on X13 and X15. The comments here are not updated, then check
%the ones in "Drone_Control_CASY".

Fullsys = ss(A,B,C,0);
FullsysDist = c2d(Fullsys, 0.004, 'Tustin');
G = tf(Fullsys); %G is a matrix-like transfer function, as it is actually composed
%by 16x4 tf. As matter of fact, we have for each input (4 in total) a tf
%for each output (16 in total)

%Now I need to generate the filters for the inputs and for the combination
%of states (y), which in this case are just the states, as can be seen
%studying the matrix C that is an Identity matrix. 
%This means that for the 'y' I have to generate 16 filters, one per state.


%FILTER FOR X1, I.E. THE POSITION ALONG X
%The filter I will implement is an high-pass filter, in this way actually we 
%will filter-out the high frequencies, as with the LQR control we have to 
%minimize a cost function J where the component with an higher weight are
%the ones more penalized.

%My project is to implement a control that allows the drone to performe a
%circumference at constant height. For this reason the filter for the X and
%Y position must have a resonant peak at a frequency strictly related to
%the angular velocity with which we want to performe the circle:
%omega=2*pi*frequency. In this way the position along X and Y will be
%forced to have a circular behavior. 
%Since the FS LQR penalized the frequency with an higher weight, the resonant
%peak has to be negative. In this way all the frequency, except for the one
%choosen, will be penalized and then eliminated.

%Another consideration is keeping the gain of the whole frequency domain 
%lower than the gains of Roll, Pitch & Yaw, and respectives angular 
%velocities, as my first goal is to maintain the system in a linear region 
%of operation.


omega1 = 3; %I assume a radius of 0.5m with a tangential velocity of 1.5m/s
gain1 = 500;
Pq1 = tf(gain1, [1, 0, omega1^2]);
% Pq1 = tf(gain1*[1, 0, (omega1)^2], [1, 7, 10]); %The denominator is
% %randomly chosen in order to obtain a feasible transfer function (relative
% %degree <= 0), and being careful to get a stable filter.

bode(Pq1); grid on; legend;

[Aq1,Bq1,Cq1,Dq1] = tf2ss(gain1, [1, 0, omega1^2]); 

% [Aq1,Bq1,Cq1,Dq1] = tf2ss(gain1*[1, 0, (omega1)^2], [1, 7, 10]); 
%Matrices A,B,C,D for this filter


%FILTER FOR X2, I.E. THE VELOCITY ALONG X
%For the velocity along X we do not have any special request, then a
%constant gain will be applied.

gain2 = 20; %ten times higher (10dB) than the linear velocities 
Pq2 = tf(gain2, 1);

bode(Pq2); grid on; legend;

[Aq2,Bq2,Cq2,Dq2] = tf2ss(gain2, 1);%Matrices A,B,C,D for this filter


%FILTER FOR X3, I.E. THE POSITION ALONG Y
%Same filter of X1 for the same reason

omega3 = 3; %I assume a radius of 0.5m with a tangential velocity of 1.5m/s
gain3 = 500;
Pq3 = tf(gain3, [1, 0, omega3^2]);
% Pq3 = tf(gain3*[1, 0, (omega3)^2], [1, 7, 10]); %The denominator is
% % %randomly chosen in order to obtain a feasible transfer function (relative
% % %degree <= 0), and being careful to get a stable filter.

bode(Pq3); grid on; legend;

[Aq3,Bq3,Cq3,Dq3] = tf2ss(gain3, [1, 0, omega3^2]); 
% [Aq3,Bq3,Cq3,Dq3] = tf2ss(gain3*[1, 0, (omega3)^2], [1, 7, 10]); 
%Matrices A,B,C,D for this filter



%FILTER FOR X4, I.E. THE VELOCITY ALONG Y
%Same filter of X2 for the same reason

gain4 = 20; %ten times higher (10dB) than the linear velocities 
Pq4 = tf(gain4, 1);

bode(Pq4); grid on; legend;

[Aq4,Bq4,Cq4,Dq4] = tf2ss(gain4, 1);%Matrices A,B,C,D for this filter


%FILTER FOR X5, I.E. THE POSITION ALONG Z
%Here the filter we want is something that acts as an integrator. Indeed if
%we want to control our drone in order to let it make circle in air at a
%constant height we must have a zero error for the position along z. To
%achieve this result we must be sure that a constant disturb (such as a
%wrong data for the mass), will be removed, then an integrator will be
%implemented.
%Since the FS LQR penalized the frequency with an higher weight, if we
%implement an integrator, this will give an infinite weight to the
%component at omega=0, then the FS LQR will penalized this frequency the
%most.

gain5 = 20;
Pq5 = tf(gain5, [1, 0]);
bode(Pq5); grid on; legend;

[Aq5,Bq5,Cq5,Dq5] = tf2ss(gain5, [1, 0]); 
%Matrices A,B,C,D for this filter

%FILTER FOR X6, I.E. THE VELOCITY ALONG Z 
%same thing of X2

gain6 = 20; %ten times higher (10dB) than the linear velocities 
Pq6 = tf(gain6, 1);

bode(Pq6); grid on; legend;

[Aq6,Bq6,Cq6,Dq6] = tf2ss(gain6, 1);%Matrices A,B,C,D for this filter


%FILTER FOR X7, I.E. THE ROLL ANGLE
%What we need for RPY and their velocity is a filter that gives to them an
%high gain, then an high weight, in the whole frequency range. In this way,
%the cost function will try to minimize their variations. We want this to
%keep the model inside the linear region.
%Actually the linearity of RPY is the most important, then the gain for
%them will be slightly higher than the one on the angular velocities.

gain7 = 70; %ten times higher (10dB) than the linear velocities 
Pq7 = tf(gain7, 1);

bode(Pq7); grid on; legend;

[Aq7,Bq7,Cq7,Dq7] = tf2ss(gain7, 1);%Matrices A,B,C,D for this filter


%FILTER FOR X8, I.E. THE PITCH ANGLE
%Same of X7

gain8 = 70; %ten times higher (10dB) than the linear velocities 
Pq8 = tf(gain8, 1);

bode(Pq8); grid on; legend;

[Aq8,Bq8,Cq8,Dq8] = tf2ss(gain8, 1);%Matrices A,B,C,D for this filter


%FILTER FOR X9, I.E. THE YAW ANGLE
%Same of X7

gain9 = 70; %ten times higher (10dB) than the linear velocities 
Pq9 = tf(gain9, 1);

bode(Pq9); grid on; legend;

[Aq9,Bq9,Cq9,Dq9] = tf2ss(gain9, 1);%Matrices A,B,C,D for this filter


%FILTER FOR X10, I.E. THE ANGULAR VELOCITY ASSOCIATED TO ROLL
%Same of X7. The only difference is that on the velocity I choose a gain
%slightly lower than the one on RPY, as the linearity of RPY is more
%relevant to keep the drone in the linear model.

gain10 = 50; %ten times higher (10dB) than the linear velocities 
Pq10 = tf(gain10, 1);

bode(Pq10); grid on; legend;

[Aq10,Bq10,Cq10,Dq10] = tf2ss(gain10, 1);%Matrices A,B,C,D for this filter


%FILTER FOR X11, I.E. THE ANGULAR VELOCITY ASSOCIATED TO PITCH
%Same of X10

gain11 = 50; %ten times higher (10dB) than the linear velocities 
Pq11 = tf(gain11, 1);

bode(Pq11); grid on; legend;

[Aq11,Bq11,Cq11,Dq11] = tf2ss(gain11, 1);%Matrices A,B,C,D for this filter


%FILTER FOR X12, I.E. THE ANGULAR VELOCITY ASSOCIATED TO YAW
%Same of X10

gain12 = 50; %ten times higher (10dB) than the linear velocities 
Pq12 = tf(gain12, 1);

bode(Pq12); grid on; legend;

[Aq12,Bq12,Cq12,Dq12] = tf2ss(gain12, 1);%Matrices A,B,C,D for this filter


%FILTER FOR X13, I.E. THE POSITION a OF THE POLE
%As said for RPY, our goal is to keep the system in a configuration where
%the behavoir is well described by the linear model. For this reason, a
%similar filter to the ones of RPY can be built.

gain13 = 70; %ten times higher (10dB) than the linear velocities 
Pq13 = tf(gain13, 1);

bode(Pq13); grid on; legend;

[Aq13,Bq13,Cq13,Dq13] = tf2ss(gain13, 1);%Matrices A,B,C,D for this filter


%FILTER FOR X14, I.E. THE VELOCITY ASSOCIATED TO a 
%Same of X13. But here again the linearity on the position is more
%important than the linearity on the velocity, then a lower gain will be
%associated to 'Va' and 'Vb'.

gain14 = 50; %ten times higher (10dB) than the linear velocities 
Pq14 = tf(gain14, 1);

bode(Pq14); grid on; legend;

[Aq14,Bq14,Cq14,Dq14] = tf2ss(gain14, 1);%Matrices A,B,C,D for this filter


%FILTER FOR X15,  I.E. THE POSITION b OF THE POLE
%Same of X13

gain15 = 70; %ten times higher (10dB) than the linear velocities 
Pq15 = tf(gain15, 1);

bode(Pq15); grid on; legend;

[Aq15,Bq15,Cq15,Dq15] = tf2ss(gain15, 1);%Matrices A,B,C,D for this filter


%FILTER FOR X16, I.E. THE VELOCITY ASSOCIATED TO b
%Same of X13

gain16 = 50; %ten times higher (10dB) than the linear velocities 
Pq16 = tf(gain16, 1);

bode(Pq16); grid on; legend;

[Aq16,Bq16,Cq16,Dq16] = tf2ss(gain16, 1);%Matrices A,B,C,D for this filter



%Let's see now how to built the filters for the 4 inputs. The most
%resonable design is 4 filters to eliminate the high frequency content of
%the control action. Such a filter has to amplify high frequencies, in
%this way the cost function J will minimize that contribute.

%FILTER FOR U1, I.E. THE TRUST

tau1u1 = 0.5; %In this way the zero will be in -2
tau2u1 = tau1u1/20; %In this way the pole will be in -40

Pr1 = tf(0.2*[tau1u1, 1], [tau2u1, 1]); %I calculate the tf to obtain later the
%matrices A,B,C,D
bode(Pr1); grid on; legend;

[Ar1,Br1,Cr1,Dr1] = tf2ss(0.2*[tau1u1 1],[tau2u1 1]);%Matrices A,B,C,D for this filter


%FILTER FOR U2, I.E. THE TORQUE AROUND X

tau1u2 = 0.03; %In this way the zero will be in -2 (Not true anymore)
tau2u2 = tau1u2/20; %In this way the pole will be in -40 (Not true anymore)

Pr2 = tf([tau1u2, 1], [tau2u2, 1]); %I calculate the tf to obtain later the
%matrices A,B,C,D
bode(Pr2); grid on; legend;

[Ar2,Br2,Cr2,Dr2] = tf2ss([tau1u2 1],[tau2u2 1]);%Matrices A,B,C,D for this filter


%FILTER FOR U3, I.E. THE TORQUE AROUND Y

tau1u3 = 0.03; %In this way the zero will be in -2
tau2u3 = tau1u3/20; %In this way the pole will be in -40

Pr3 = tf([tau1u3, 1], [tau2u3, 1]); %I calculate the tf to obtain later the
%matrices A,B,C,D
bode(Pr3); grid on; legend;

[Ar3,Br3,Cr3,Dr3] = tf2ss([tau1u3 1],[tau2u3 1]);%Matrices A,B,C,D for this filter


%FILTER FOR U4, I.E. THE TORQUE AROUND Z

tau1u4 = 0.03; %In this way the zero will be in -2
tau2u4 = tau1u4/20; %In this way the pole will be in -40

Pr4 = tf([tau1u4, 1], [tau2u4, 1]); %I calculate the tf to obtain later the
%matrices A,B,C,D
bode(Pr4); grid on; legend;

[Ar4,Br4,Cr4,Dr4] = tf2ss([tau1u4 1],[tau2u4 1]);%Matrices A,B,C,D for this filter

%WRONG (da eliminare se non mi dovesse servire poi)
% [m1,n1]=size(Aq1); [m2,n2]=size(Aq2); [m3,n3]=size(Aq3); [m4,n4]=size(Aq4);
% [m5,n5]=size(Aq5); [m6,n6]=size(Aq6); [m7,n7]=size(Aq7); [m8,n8]=size(Aq8);
% [m9,n9]=size(Aq9); [m10,n10]=size(Aq10); [m11,n11]=size(Aq1); [m12,n12]=size(Aq12);
% [m13,n13]=size(Aq13); [m14,n14]=size(Aq14); [m15,n15]=size(Aq15); [m16,n16]=size(Aq16);

%EXTENDED SYSTEM 
%From the theory of the FS LQR I can understand which are the structures
%of the matrices Aa, Ba, Ca that define the extended system. The problem is
%that the known structure works when we have a single matrix Aq, Bq, Ar, Br
%whereas here I have multiple of those. Then, in order to apply the
%formulas of the theory I need to built a unique matrix Aq, Bq, Ar and Br. 
%The structure of those matrices can be deducted thinking about the result 
%we are trying to achive. What we need are diagonal matrices composed, on
%the diagonal, by the single matrices of the filters (e.g. Aq1, Aq2, Aq3
%etc..).

%Actually, as many matrix Aq_i do not exist, the matrix Aq is just 5x5.
%This are the states of the actual filter on the states of the system.
%Indeed, only on x1, x3, and x5 I have designed a real filter. For the
%other states I just used a constant gain that is pretty useless as those
%gains where designed for a different scenario that is not anymore used.
%Those gains are replaced by the matrix Q (see below).
Aq = blkdiag(Aq1,Aq2,Aq3,Aq4,Aq5,Aq6,Aq7,Aq8,Aq9,Aq10,Aq11,Aq12,Aq13,Aq14,Aq15,Aq16);
[mmm,nnn] = size(Aq);

%For the same reason explained for Aq, also the matrix Bq is just a 5x3
%matrix.
Bq = blkdiag(Bq1,Bq2,Bq3,Bq4,Bq5,Bq6,Bq7,Bq8,Bq9,Bq10,Bq11,Bq12,Bq13,Bq14,Bq15,Bq16);
[jjj,kkk] = size(Bq);

%Matrix used to select the state of the system to be filtered, i.e. x,y,z.
%This matrix is employed in the composition of the matrix Aa. In this way
%doing Bq*CC (so in matrix dimension 5x3*3x16) we get a 5x16 matrix that is
%suitable to interface the states of the 3 state filters (on x,y and z) with
%the 16 states of the system (drone+pole).
CC = [1, zeros(1,15);
      0, 0, 1, zeros(1,13);
      0, 0, 0, 0, 1, zeros(1,11)];
[eta,ni] = size(CC);

%For Cq and Dq the same approach used for Aq and Bq is employed, the only
%difference is that in the following case only the relevant matrices are
%inserted just for sake of semplicity (no real difference is present
%though).
Cq = [Cq1, 0, 0, 0;
      zeros(1, 5);
      0, 0, Cq3, 0;
      zeros(1, 5);
      0, 0, 0, 0,Cq5;
      zeros(11,5)];
[teta2,rho2] = size(Cq);
Dq = blkdiag(Dq1,Dq2,Dq3,Dq4,Dq5,Dq6,Dq7,Dq8,Dq9,Dq10,Dq11,Dq12,Dq13,Dq14,Dq15,Dq16);
[teta,rho] = size(Dq);

%WRONG (da eliminare se non mi dovesse servire poi)
%As many matrices Aq_i actually do not exist, if I put them as input of the
%matlab command 'blkdiag', this will neglect them. This will cause problem
%of dimensions when we try to built the matrices of the extend system (Aa,
%Ba and Ca). The solution is substitute those inexistent matrices with a
%zero, therefore a row of zero will be introduced in that position. 
% 
% Aqq = blkdiag(Aq1,0,Aq3,0,Aq5,0,0,0,0,0,0,0,0,0,0,0);
% [mm,nn] = size(Aqq);
% Bqq = blkdiag(Bq1,0,Bq3,0,Bq5,0,0,0,0,0,0,0,0,0,0,0);
% [jj,kk] = size(Bqq);
% Cqq = blkdiag(Cq1,0,Cq3,0,Cq5,0,0,0,0,0,0,0,0,0,0,0);
% [m,n] = size(Cqq);
% Dqq = blkdiag(Dq1,Dq2,Dq3,Dq4,Dq5,Dq6,Dq7,Dq8,Dq9,Dq10,Dq11,Dq12,Dq13,Dq14,Dq15,Dq16);
% [j,k] = size(Dqq);

%For this system Ca is nothing but a indentity matrix 38x38
% Ca = [C, zeros(16, nn+ss);
%       zeros(mm, 16), eye(mm), zeros(mm, ss);
%       zeros(tt, 16+nn), eye(tt)];
% FINE DEL WRONG

%For this 4 matrices the approach is quite straight-forward once understood
%the logic behind Aq, Bq, Cq and Dq.
Arr = blkdiag(Ar1,Ar2,Ar3,Ar4);
[tt,ss] = size(Arr);
Brr = blkdiag(Br1,Br2,Br3,Br4);
Crr = blkdiag(Cr1,Cr2,Cr3,Cr4);
[t,s] = size(Crr);
Drr = blkdiag(Dr1,Dr2,Dr3,Dr4);
[e,f] = size(Drr);

%Once the previous matrices are correctly defined, the generation of Aa,
%Ba and Ca is quite straight-forward using the standard formulas. The only
%task to be done is to correctly choose the zeros dimensions.  
Aa = [A, zeros(16, nnn+ss);
      Bq*CC, Aq, zeros(mmm, ss)
      zeros(tt, 16+nnn), Arr];
[ee,ff] = size(Aa);


Ba = [B;
      zeros(mmm, ss);
      Brr];
[gg,hh] = size(Ba);
 
%(To understand this part, check also the notes on GoodNotes pg. 76-77-78)
%I will built now other two matrix: Ca and Cb. Ca is used to performe
%ya=Ca*xa, and so is used to select 12 states (3 of the drone, 5
%related to the 3 filters built for x,y,z, finally 4 for the states of the 
%filters of the inputs).
%Cb is used to performe y=Cb*xa and is designed to collect all the states not
%included in ya.

%With the first 3 rows I select x, y and z, then with the following 5
%rows I select all the states Zq, finally with the last 4 rows I select Zr.
% Ca = [1, zeros(1,24);
%       0, 0, 1, zeros(1,22);
%       0, 0, 0, 0, 1, zeros(1,20);
%       zeros(5,16), eye(5), zeros(5,4);
%       zeros(4,16+5), eye(4)];
Ca = eye (25);
  
%With the first 3 rows I select Vx, Vy and Vz, then with the following 10
%rows I select all the remaining states of the system drone+pole.
% Cb = [0, 1, zeros(1,23);
%       0, 0, 0, 1, zeros(1,21);
%       0, 0, 0, 0, 0, 1, zeros(1,19);
%       zeros(10,6), eye(10), zeros(10,5+4)];

%Now we have all the tools to define the equivalent Q, R, N matrices (Qa, 
%Ra, Na) used in the cost function J of the extended system
%Note that with the approach we are using (the one described in the notes 
%on GoodNotes pg. 76-77-78, Qa includes the weights only for the states
%related to the filters, then x,y,z, 5 Zq and 4 Zr)
Qa = [(Dq')*Dq, (Dq')*Cq, zeros(rho, s);
      (Cq')*Dq, (Cq')*Cq, zeros(rho2, s);
      zeros(s, rho+rho2), (Crr')*Crr];
  
Na = [zeros(25-4, f); (Crr')*Drr]; 
Ra = Drr'*Drr;

%The matrices just found are useful to define the cost function of the
%extended system, but this cost fuction as a problem very similar to a
%cross-coupled LQR, then we can define a new A, Q and u (ABar, QBar, uBar)
%that allows us to redefine the same problem without the cross-coupling
%between x and u.

ABar = Aa-Ba*(inv(Ra))*(Na')*Ca;
QBar = Qa-Na*(inv(Ra))*(Na');

%Checking the controllability of the system 
digitsOld = digits(40); %This increases the number of significant digits 
%after the coma
COa = ctrb(vpa(ABar), vpa(Ba)); %vpa is used to performe calculations with 
%the number of digits previously selected
disp(rank(vpa(COa)));
%disp(rank(ctrb(ABar, Ba)));


%Even if thru the previous commands seems that the obtained system is not
%controllable, as the rank is not maximum, this is not true. Indeed Matlab,
%probably because the computations are really hard gets confused. To
%demonstrate that the rank is actually 25 I will "manually" compute the
%rank of a part of the controllability matrix. Doing so is possible to see
%that for all the components showed below the rank results to be 25. Then,
%just completing the controllability matrix, so just adding new colums, the
%rank starts decreasing, and this is absurd. Then we can say that 
%ctrb(ABar, Ba) as rank = 25.

% rank([Ba, Aa*Ba, Aa^2*Ba, Aa^3*Ba, Aa^4*Ba, Aa^5*Ba, Aa^6*Ba, Aa^7*Ba, Aa^8*Ba, Aa^9*Ba, Aa^10*Ba])
rank([Ba, ABar*Ba, ABar^2*Ba, ABar^3*Ba, ABar^4*Ba, ABar^5*Ba, ABar^6*Ba, ABar^7*Ba, ABar^8*Ba, ABar^9*Ba, ABar^10*Ba, ABar^11*Ba, ABar^12*Ba, ABar^13*Ba, ABar^14*Ba, ABar^15*Ba, ABar^16*Ba, ABar^17*Ba, ABar^18*Ba]);
rank([Ba, ABar*Ba, ABar^2*Ba, ABar^3*Ba, ABar^4*Ba, ABar^5*Ba, ABar^6*Ba, ABar^7*Ba]);
%CHECK A MANO DEL RANK
% H = sqrt(QBar);
H = [Dq, Cq, zeros(16,4)];
disp('Detectability FLQ:');
disp(rank(obsv(Aa, H*Ca)));

%{
H = [Dqq, Cqq, zeros(j, s)];

disp('Stabilizability FLQ:');
disp(rank(ctrb(ABar, Ba)));
% disp('Detectability FLQ:');
% disp(rank(obsv(ABar, H*Ca)));

 disp('Stabilizability FLQ:');
 disp(rank(ctrb(ABar, Ba)));

%To find uBar = -KBar*x_a, I must find KBar first (KBar = ((Ra)^(-1))*Ba'*S).
%But to find KBar I need to know S, that is the solution of the ARE. So
%thanks to the command 'are' I can determine the solution S:


[Kf,S,CLP] = lqr(Aa, Ba, (Ca')*Qa*Ca, Ra, Ca'*Na);

%}

% S = are(ABar, Ba*(inv(Ra))*(Ba'), Ca'*QBar*Ca);

%The vector q contains all the weights for the states of the filter that
%are not filtered, but the still need to be controlled be the cost function
%J. The logic behind the chosen values are the same explained for the
%'filters' (actually just a gain) previously uselessly designed for these
%states.

% q = [40 40 40 900 900 900 60 60 60 900 60 900 60];

%To be used in the cost function the vector q must be transform in a
%diagonal matrix.

% Q = diag (q);

%Qe is an extended matrix that takes into account both the weights of the
%filtered states (QBar) and of the un-filtered once (Q). The zeros are
%needed to make it squared.

% Qe = [QBar, zeros(12, 13);
%       zeros(13, 12), Q];
  
%Now that I have a matrix Qe that includes the weights for the whole
%extended system, I have also to define a matrix Ce that includes the
%caracteristics of the two matric C, I have already used)

% Ce = [Ca;
%       Cb];

%Now that the extended system has been fully described, and all the
%problems has been solved, the Algebraic Riccati Equation can be solved.
S = are(ABar, Ba*(inv(Ra))*(Ba'), Ca'*QBar*Ca);


%If we prefer to use the command 'lqr' this can be done. Be careful because
%we have to select as matrix N, the matrix zeros(25,4), because we have
%already solve the cross-coupled problem 'manually'. Indeed the command
%'lqr' foresee to solve by itself the cross-coupled problem.
%[Kf,S,CLP] = lqr(ABar, Ba, (Ce')*QBar*Ca, Ra, zeros(25,4));

KBar = inv(Ra)*(Ba')*S;
disp(eig(ABar-Ba*KBar));

%See the slide about the frequency shaped controll to understand why Ka is
%define like this.
Ka = KBar+inv(Ra)*Na'*Ca;

%If all the eigenvalues have a negative real part, then the system is
%stabilized
disp(eig(Aa-Ba*Ka));

%In the following lines we have an example of LQR control
q_lqr = [1, 1, 1, 1, 1, 1, 15, 15, 15, 10, 10, 10, 15, 10, 15, 10];
Q_lqr = diag (q_lqr);

r_lqr = [2, 2, 2, 2];
R_lqr = diag (r_lqr);
[K_lqr,S,CLP] = lqr(A, B, Q_lqr, R_lqr, zeros(16,4));


[K_dlqr,S_d,CLP_d] = dlqr(A, B, Q_lqr, R_lqr, zeros(16,4));

piru_I_II = [0 -1 0;
            1 0 0;
            0 0 1];

piru_B_I = [-1 0 0;
            0 1 0;
            0 0 -1];
        
        
piru_B_II = piru_B_I * piru_I_II;
 