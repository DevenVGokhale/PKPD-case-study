[ plugin ] mrgx

[ param ]
TVCL   = 5    // Clearance (L/h)
TVV1   = 20   // Central volume (L)
TVQ    = 2    // Intercompartmental clearance (L/h)
TVV2   = 10   // Peripheral volume (L)

TVKG   = 0.05  // Tumor growth rate (1/day)
TVKS   = 0.3   // Drug killing rate (1/day)
TVLAMBDA = 0.01 // Resistance rate (1/day)

BSA    = 1.8   // Body surface area (m²), or other scaling covariate
IC50   = 1.0   // Concentration for 50% effect (mg/L)

ETA_CL = 0
ETA_V1 = 0
ETA_KG = 0
ETA_KS = 0

[ init ]
CENT   = 0
PERIPH = 0
TUM = 100     // Baseline tumor size (mm)
RES = 0.2       // Resistance effect (unitless)

[ main ]
double CL = TVCL * exp(ETA_CL);
double V1 = TVV1 * exp(ETA_V1);
double Q  = TVQ;
double V2 = TVV2;

double KG = TVKG * exp(ETA_KG);
double KS = TVKS * exp(ETA_KS);
double LAMBDA = TVLAMBDA;

double IC50_eff = IC50;

[ ode ]
double CP   = CENT / V1;
double EFF  = KS * CP / (IC50_eff + CP) * exp(-RES);

double dxdt_CENT   = - (CL / V1) * CENT - Q/V1 * CENT + Q/V2 * PERIPH;
double dxdt_PERIPH = Q/V1 * CENT - Q/V2 * PERIPH;
double dxdt_TUM    = KG * TUM - EFF * TUM;
double dxdt_RES    = LAMBDA;

[ capture ]
CP CL V1 KG KS
