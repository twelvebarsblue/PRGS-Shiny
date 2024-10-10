$PROB

$PARAM @annotated @covariates
CRCL : 90 : CRCL
AGE : 20 : AGE
TIME_1: 12 : Time at first dose
TIME_2 : 24 : Time at second dose
TIME_3 : 36 : Time at third dose

II_1: 12 : Dose interval between T0 and T1
II_2: 24 : Dose interval between T1 and T2
II_3: 36 : Dose interval between T3 and T4

DOSE_1: 1000 : Initial dose (mg)
DOSE_2: 1000 : Second dose
DOSE_3 : 1000 : Third dose
DOSE_4 : 1000 : Fourth dose

UNTIL : 160 : Until when

$PARAM @annotated
ETA1 : 0 : CL (L/h)
ETA2 : 0 : VC (L)

$OMEGA 0.0099940932 0.08005596
$SIGMA 0.2679 0.002647

$CMT @annotated
CENTRAL : Central compartment (mg/L) [ADM, OBS]
AUC : Area under the curve

$TABLE
double CP = CENTRAL/V ;
double DV  = CP*(1+EPS(1)) + EPS(2) ;

$PLUGIN evtools

$PREAMBLE
reg.init(self);

$GLOBAL
double CRCL, AGE, TVCL, CL, V, K10, AUC, INTERVAL, DOSE_1, DOSE_2, DOSE_3, DOSE_4;
evt::regimen reg;

$THETA @annotated
2.828  : TVCL (initial estimate)
52.14 : TVV
0.00842 : CRCL
0.08143 : AGE

$PK
TVCL = THETA(1) * (1+THETA(3) * (CRCL-80) * std::pow((75/AGE), THETA(4))) ;
CL = TVCL*exp(ETA(1) + ETA1) ;
V = THETA(2) * exp(ETA(2) + ETA2)  ;
K10 = CL/V;


if(evt::near(TIME, TIME_1) && EVID > 0) {
  reg.ii(II_2);
  reg.amt(DOSE_2);
}

if(evt::near(TIME, TIME_2) && EVID > 0) {
  reg.ii(UNTIL);
  reg.amt(DOSE_3);
}

if(evt::near(TIME, TIME_3) && EVID > 0) {
  reg.amt(DOSE_4);
}


if(NEWIND <= 1) {
  reg.init(self);
  reg.amt(DOSE_1);
  reg.ii(II_1);
  reg.cmt(1);
  reg.until(UNTIL);
  reg.flagnext();
}

$ERROR
reg.execute();

$SET ss_cmt = "-AUC"

$ODE
dxdt_CENTRAL    = - K10*CENTRAL;
dxdt_AUC = CENTRAL/V;

$CAPTURE CP DV CL DOSE_1 DOSE_2 DOSE_3 DOSE_4  
