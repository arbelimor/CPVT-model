%clear all

constants_CPVT

%initial conditions

Y0(1)=-64.5216286940; %Vm
Y0(2)=0.6246780312; %qa
Y0(3)=0.4537033169; %qi
Y0(4)=0.0016256324;	 %dT
Y0(5)=0.4264459666; %fT
Y0(6)=0.4043600437; %pa
Y0(7)=0.9250035423; %pi
Y0(8)=0.0127086259; %xs
Y0(9)=0.9968141226; %fL12
Y0(10)=0.0000045583; %dL12
Y0(11)=0.9809298233; %fL13
Y0(12)=0.0002036671; %dL13
Y0(13)=0.7649576191; %fCa
Y0(14)=0.0046263658; %r
Y0(15)=0.6107148187; %q
Y0(16)=0.4014088304; %m15
Y0(17)=0.2724817537; %h15
Y0(18)=0.0249208708; %j15
Y0(19)=0.1079085266; %m11
Y0(20)=0.4500098710; %h11
Y0(21)=0.0268486392; %j11
Y0(22)=0.0279984462; %y
Y0(23)=0.0000319121; %Cai
Y0(24)=0.1187281829; %CajSR
Y0(25)=1.5768287365; %CanSR
Y0(26)=0.0000560497; %Casub
Y0(27)=0.0063427103; %fTC
Y0(28)=0.1296677919; %fTMC
Y0(29)=0.7688656371; %fTMM
Y0(30)=0.0242054739; %fCMs
Y0(31)=0.0138533048; %fCMi
Y0(32)=0.1203184861; %fCQ
Y0(33)=0.7720290515; %R
Y0(34)=0.0000000760; %OO
Y0(35)=0.0000000213; %S
Y0(36)=0.2162168926; %RI
Y0(37)=0.0004; %w
Y0(38)=19.73; %cAMP
Y0(39)=0.23; %PLB
Y0(40)=0.06;   %A
Y0(41)=0.02;   %TT
Y0(42)=0.06;   %U
Y0(43)=1.75e-6;   %SL

%Model parameters

prompt = 'Insert total time (ms)';
t_end = input(prompt);

%ODE solver
%Mutations code: WT=0, R4496C=1, Casq2=2, HCN4=3

prompt = 'Insert mutations code (WT=0, R4496C=1, Casq2=2, HCN4=3)';
mutation = input(prompt);

prompt = 'Insert ISO (mM)';
con.ISO = input(prompt);

%Therapy
con.IBMX=0; % PDE blocker therapy
con.P_up=0.02; %was 0.02 - SERCA therapy
con.K_AC = 0.0735; %was 0.0735  [min-1] Non-Ca2+ AC activity 

if mutation==0
    con.CQ_tot = 10;
    con.k_om = 0.06; %[1/ms]
    con.n=2;
elseif mutation==1 && con.ISO==0 %R4496C, no ISO
    con.k_om = 0.045; %[1/ms]
    con.n=1.9;
elseif mutation==1 && con.ISO>0 %R4496C + ISO
    con.k_om = 0.045; %[1/ms]
    con.n=1.7;
elseif mutation==2 && con.ISO==0 %Casq2
    con.CQ_tot = 2.6;
    con.k_fCQ = 0.534;
elseif mutation==2 && con.ISO==10 %Casq2 + ISO
    con.CQ_tot = 2.6;
    con.k_fCQ = 0.3;
elseif mutation==2 && con.ISO==20 %Casq2 + ISO
    con.CQ_tot = 2.6;
    con.k_fCQ = 0.15;
elseif mutation==4
    %con.PKA_PLB=0.01;
    con.n_PLB = 1.5; % was 1

end
 
[t_ms,Y]=ode15s(@(t,Y)diff_eq_CPVT(t,Y,con,mutation),0:t_end,Y0);
I=CPVT_equations(Y,con,0);

%ms to sec
t=t_ms/1000;

% Naming currents

It0=I(:,1);
Isus=I(:,2);
Ib_Ca=I(:,3);
Ib_K=I(:,4);
Ib_Na=I(:,5);
ICaL=I(:,6);
ICaT=I(:,7);
If=I(:,8);
IK1=I(:,9);
IKr=I(:,10);
IKs=I(:,11);
INa=I(:,12);
INaK=I(:,13);
INaCa=I(:,14);
I_st=I(:,15);
I_KACh=I(:,16);
jSRCarel=I(:,17);
j_up=I(:,18);
j_Cadiff=I(:,19);
j_tr=I(:,20);
ATP=I(:,21);
PKA=I(:,22);

% Naming variables

Vm=Y(:,1); %[mV] Membrane potential
qa=Y(:,2); %Activation gating variable of I_st
qi=Y(:,3); %Inactivation gating variable of I_st
dT=Y(:,4); %Inactivation gate of I_CaT
fT=Y(:,5); %Activation gate of I_CaT
pa=Y(:,6); %Activation gating variable of I_Kr
pi=Y(:,7); %Inactivation gating variable of I_Kr
xs=Y(:,8); %Activation gating variable of I_Ks
fL12=Y(:,9); %Inactivation gate for I_(C aL,1.2)
dL12=Y(:,10); %Activation gate for I_(CaL,1.2)
fL13=Y(:,11); %Inactivation gate for I_(CaL,1.3)
dL13=Y(:,12); %Activation gate for I_(CaL,1.3)
fCa=Y(:,13); %Ca2+-dependent inactivation gating variable for I_(CaL,1.2)and I_(CaL,1.2)
r=Y(:,14); %Activation gating variable of I_toand I_sus
q=Y(:,15); %Inactivation gating variable of I_to
m15=Y(:,16); %Activation gating variable of Nav1.5
h15=Y(:,17); %Fast inactivation gating variable of Nav1.5
j15=Y(:,18); %Slow inactivation gating variable of Nav1.5
m11=Y(:,19); %Activation gating variable of Nav1.1
h11=Y(:,20); %Fast inactivation gating variable of Nav1.1
j11=Y(:,21); %Slow inactivation gating variable of Nav1.1
y=Y(:,22); %Activation gating variable of I_f
Cai=Y(:,23); %[mM] Intracellular Ca2+ concentration or Ca2+ concentration in the cytosol
CajSR=Y(:,24); %[mM] Ca2+ concentration in the JSR
CanSR=Y(:,25); %[mM] Ca2+ concentration in the NSR
Casub=Y(:,26); %[mM] Ca2+ concentration in the subspace
f_TC=Y(:,27); %Fractional occupancy of the troponin Ca2+ site by [Ca2+]i
f_TMC=Y(:,28); %Fractional occupancy of the troponin Mg2+ site by [Ca2+]i
f_TMM=Y(:,29); %Fractional occupancy of the troponin Mg2+ site by Mg2+
f_CMs=Y(:,30); %Fractional occupancy of calmodulin by [Ca2+]sub
f_CMi=Y(:,31); %Fractional occupancy of calmodulin by [Ca2+]i
f_CQ=Y(:,32); %Fractional occupancy of calsequestrin by [Ca2+]rel
R=Y(:,33); %Fraction of reactivated (closed) RyR channels
OO=Y(:,34); %Open fraction of RyR channels
S=Y(:,35); %Inactive fraction of RyR channels
RI=Y(:,36); %Fraction of RyR inactivated channels
w=Y(:,37); %IKACh ACh and voltage-dependent gating variable
cAMP=Y(:,38); %cAMP
PLB=Y(:,39); %PLB phosphorylation level
A=Y(:,40); %The density of regulatory units with bound Ca2+ and adjacent weak cross-bridges
TT=Y(:,41); %The density of regulatory units with bound Ca2+ and adjacent strong cross-bridge
U=Y(:,42); %The density of regulatory units without bound but with adjacent strong cross-bridge
SL=Y(:,43); %Sarcomere length