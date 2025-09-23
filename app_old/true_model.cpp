///////////////////////////////////////////////////////////
// Model: Tislelizumab like PK + Cytokine PD Model
// Author: Deven Gokhale
// Date: 2025-08-2
// Description:
// - Three-compartment PK model for tislelizumab (IV dosing)
// - PK parameters include CL, VC, Q2, V2, Q3, V3
// - Inter-individual variability (IIV) handled via idata
// - Indirect response PD model for cytokine (e.g., IFN-γ)
//   -> Cytokine production is stimulated by CP via Emax model
//   -> No drug elimination via PD mechanism
//
// PK equations:
//   dCENT/dt     = input - elimination - distribution to P1/P2
//   dPERIPH1/dt  = distribution to/from CENT
//   dPERIPH2/dt  = distribution to/from CENT
//
// PD equation:
//   dCYT/dt = KIN * (1 + CP / (EC50 + CP)) - KOUT * CYT
//
// Units:
//   - Doses: mg
//   - Concentration: µg/mL
//   - KIN: pg/mL/day
//   - Time: days
///////////////////////////////////////////////////////////

[ plugin ] 
mrgx

[ param ]
// PK parameters
CL   = 0.153   // Clearance (L/day)
VC   = 3.05    // Central volume (L)
Q2   = 0.74    // Intercompartmental clearance 1 (L/day)
V2   = 1.27    // Peripheral volume 1 (L)
Q3   = 0.092   // Intercompartmental clearance 2 (L/day)
V3   = 2.10    // Peripheral volume 2 (L)

// PD parameters (cytokine model)
KIN   = 10     // Cytokine production rate (pg/mL/day)
KOUT  = 1      // Cytokine degradation rate (1/day)
EC50  = 5      // Drug concentration for 50% effect (µg/mL)

[ init ]
CENT     = 0    // Central compartment
PERIPH1  = 0    // Peripheral compartment 1
PERIPH2  = 0    // Peripheral compartment 2
CYTi      = 5    // Baseline cytokine concentration

[ main ]
double CP = CENT / VC;  // Central plasma concentration (µg/mL)

[ ode ]
// Declare state derivatives
// double dxdt_CENT, dxdt_PERIPH1, dxdt_PERIPH2, dxdt_CYTi;

// PK model (3-compartment)
dxdt_CENT     = - (CL / VC) * CENT 
                - (Q2 / VC) * CENT + (Q2 / V2) * PERIPH1 
                - (Q3 / VC) * CENT + (Q3 / V3) * PERIPH2;

dxdt_PERIPH1  = (Q2 / VC) * CENT - (Q2 / V2) * PERIPH1;
dxdt_PERIPH2  = (Q3 / VC) * CENT - (Q3 / V3) * PERIPH2;

// PD model (indirect response: stimulated production)
dxdt_CYTi = KIN * (1 + CP / (EC50 + CP)) - KOUT * CYTi;

[ capture ]
CP CYT=CYTi CL VC Q2 V2 Q3 V3 KIN KOUT EC50