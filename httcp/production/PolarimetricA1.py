import os
from columnflow.util import maybe_import

np     = maybe_import("numpy")
ak     = maybe_import("awkward")
coffea = maybe_import("coffea")

from httcp.production.TComplex import TComplex
from IPython import embed


class PolarimetricA1:
    def __init__(self,
                 p4_tau    : ak.Array,
                 p4_os_pi  : ak.Array,
                 p4_ss1_pi : ak.Array,
                 p4_ss2_pi : ak.Array,
                 taucharge : ak.Array) -> None:
        """
        Calculate Polarimetric vector for tau to a1 decay
        All the vectors in the arguments must be wrt the rest frame
        """

        self.p4_tau         =  p4_tau       # tau
        self.p4_os_pi       =  p4_os_pi     # os-pion
        self.p4_ss1_pi      =  p4_ss1_pi    # ss1-pion
        self.p4_ss2_pi      =  p4_ss2_pi    # ss2-pion
        
        self.mpi            =  0.13957018   # GeV
        self.mpi0           =  0.1349766    # GeV
        self.mtau           =  1.776        # GeV
        self.coscab         =  0.975
        self.mrho           =  0.773        # GeV
        self.mrhoprime      =  1.370        # GeV
        self.ma1            =  1.251        # GeV
        self.mpiprime       =  1.300        # GeV
        self.Gamma0rho      =  0.145        # GeV
        self.Gamma0rhoprime =  0.510        # GeV
        self.Gamma0a1       =  0.599        # GeV
        self.Gamma0piprime  =  0.3          # GeV
        self.fpi            =  0.093        # GeV
        self.fpiprime       =  0.08         # GeV
        self.gpiprimerhopi  =  5.8          # GeV
        self.grhopipi       =  6.08         # GeV
        self.beta           = -0.145
        self.COEF1          =  2.0*np.sqrt(2.)/3.0
        self.COEF2          = -2.0*np.sqrt(2.)/3.0
        # C AJW 2/98: Add in the D-wave and I=0 3pi substructure:
        self.COEF3          =  2.0*np.sqrt(2.)/3.0
        self.SIGN           = -taucharge
        self.doSystematic   =  False
        self.systType       =  "UP"


    def PVC(self) -> ak.Array:
        P  = self.p4_tau
        q1 = self.p4_ss1_pi
        q2 = self.p4_ss2_pi
        q3 = self.p4_os_pi

        a1 = q1+q2+q3

        N = P.subtract(a1)

        s1 = (q2+q3).mass2
        s2 = (q1+q3).mass2
        s3 = (q1+q2).mass2

        # Three Lorentzvector: Why??? : No idea!!!
        getvec = lambda a,b,c: b - c - a*(a.dot((b - c))*(1/a.mass2))

        vec1 = getvec(a1, q2, q3)
        vec2 = getvec(a1, q3, q1)
        vec3 = getvec(a1, q1, q2)
        
        F1 = TComplex(self.COEF1)*self.F3PI(1, a1.mass2, s1, s2)
        F2 = TComplex(self.COEF2)*self.F3PI(2, a1.mass2, s2, s1)
        F3 = TComplex(self.COEF3)*self.F3PI(3, a1.mass2, s3, s1)

        HADCUR = []

        HADCUR.append(TComplex(vec1.energy)*F1 + TComplex(vec2.energy)*F2 + TComplex(vec3.energy)*F3)
        HADCUR.append(TComplex(vec1.px)*F1 + TComplex(vec2.px)*F2 + TComplex(vec3.px)*F3)
        HADCUR.append(TComplex(vec1.py)*F1 + TComplex(vec2.py)*F2 + TComplex(vec3.py)*F3)
        HADCUR.append(TComplex(vec1.pz)*F1 + TComplex(vec2.pz)*F2 + TComplex(vec3.pz)*F3)

        HADCURC = [val.Conjugate() for val in HADCUR]

        CLV = self.CLVEC(HADCUR, HADCURC, N)
        CLA = self.CLAXI(HADCUR, HADCURC, N)

        pclv    = P.dot(CLV)
        pcla    = P.dot(CLA)
        omega   = pclv - pcla
        A       = (P.mass)**2
        CLAmCLV = CLA.subtract(CLV)
        
        out = ((P.mass)*(P.mass)*(CLA-CLV) - P*(P.dot(CLA) - P.dot(CLV)))*(1/omega/P.mass)

        return out

    
    def F3PI(self,
             IFORM: float,
             QQ: ak.Array,
             SA: ak.Array,
             SB: ak.Array):
        """
            Calculate the F3PIFactor.
        """
        MRO = 0.7743
        GRO = 0.1491
        MRP = 1.370
        GRP = 0.386
        MF2 = 1.275
        GF2 = 0.185
        MF0 = 1.186
        GF0 = 0.350
        MSG = 0.860
        GSG = 0.880
        MPIZ = self.mpi0
        MPIC = self.mpi

        M1 = 0
        M2 = 0
        M3 = 0

        IDK = 1  # It is 3pi

        if IDK == 1:
            M1 = MPIZ
            M2 = MPIZ
            M3 = MPIC
        elif IDK == 2:
            M1 = MPIC
            M2 = MPIC
            M3 = MPIC

        M1SQ = M1*M1
        M2SQ = M2*M2
        M3SQ = M3*M3
        
        
        # parameter varioation for
        # systematics from https://arxiv.org/pdf/hep-ex/9902022.pdf
        db2, dph2 = 0.094, 0.253
        db3, dph3 = 0.094, 0.104
        db4, dph4 = 0.296, 0.170
        db5, dph5 = 0.167, 0.104
        db6, dph6 = 0.284, 0.036
        db7, dph7 = 0.148, 0.063

        scale = 0.0
        if self.doSystematic:
            if self.systType == "UP":
                scale = 1
            elif self.systType == "DOWN":
                scale = -1
                
        # Breit-Wigner functions with isotropic decay angular distribution
        # Real part must be equal to one, stupid polar implemenation in root
        BT1 = TComplex(1., 0.)
        BT2 = TComplex(0.12  + scale*db2, 0.) * TComplex(1, (0.99   +  scale*dph2)*np.pi, True)
        BT3 = TComplex(0.37  + scale*db3, 0.) * TComplex(1, (-0.15  +  scale*dph3)*np.pi, True)
        BT4 = TComplex(0.87  + scale*db4, 0.) * TComplex(1, (0.53   +  scale*dph4)*np.pi, True)
        BT5 = TComplex(0.71  + scale*db5, 0.) * TComplex(1, (0.56   +  scale*dph5)*np.pi, True)
        BT6 = TComplex(2.10  + scale*db6, 0.) * TComplex(1, (0.23   +  scale*dph6)*np.pi, True)
        BT7 = TComplex(0.77  + scale*db7, 0.) * TComplex(1, (-0.54  +  scale*dph7)*np.pi, True)

        F3PIFactor = None
        
        if IDK == 2:
            if IFORM == 1 or IFORM == 2:
                S1 = SA
                S2 = SB 
                S3 = QQ - SA - SB + M1SQ + M2SQ + M3SQ

                F134 = -(1 / 3.) * ((S3 - M3SQ) - (S1 - M1SQ))
                F15A = -(1 / 2.) * ((S2 - M2SQ) - (S3 - M3SQ))
                F15B = -(1 / 18.) * (QQ - M2SQ + S2) * (2 * M1SQ + 2 * M3SQ - S2) / S2
                F167 = -(2 / 3.)
    
                # Breit Wigners for all the contributions:
                FRO1 = self.BWIGML(S1, MRO, GRO, M2, M3, 1)
                FRP1 = self.BWIGML(S1, MRP, GRP, M2, M3, 1)
                FRO2 = self.BWIGML(S2, MRO, GRO, M3, M1, 1)
                FRP2 = self.BWIGML(S2, MRP, GRP, M3, M1, 1)
                FF21 = self.BWIGML(S1, MF2, GF2, M2, M3, 2)
                FF22 = self.BWIGML(S2, MF2, GF2, M3, M1, 2)
                FSG2 = self.BWIGML(S2, MSG, GSG, M3, M1, 0)
                FF02 = self.BWIGML(S2, MF0, GF0, M3, M1, 0)

                F3PIFactor = BT1*FRO1 \
                           + BT2*FRP1 \
                           + BT3*TComplex(F134, 0.)*FRO2 \
                           + BT4*TComplex(F134, 0.)*FRP2 \
                           - BT5*TComplex(F15A, 0.)*FF21 \
                           - BT5*TComplex(F15B, 0.)*FF22 \
                           - BT6*TComplex(F167, 0.)*FSG2 \
                           - BT7*TComplex(F167, 0.)*FF02

            elif IFORM == 3:
                S3 = SA 
                S1 = SB
                S2 = QQ - SA - SB + M1SQ + M2SQ + M3SQ

                F34A =  (1 / 3) * ((S2 - M2SQ) - (S3 - M3SQ))
                F34B =  (1 / 3) * ((S3 - M3SQ) - (S1 - M1SQ))
                F35A = -(1 / 18) * (QQ - M1SQ + S1) * (2 * M2SQ + 2 * M3SQ - S1) / S1
                F35B =  (1 / 18) * (QQ - M2SQ + S2) * (2 * M3SQ + 2 * M1SQ - S2) / S2
                F36A = -(2 / 3)
                F36B =  (2 / 3)

                FRO1 = self.BWIGML(S1, MRO, GRO, M2, M3, 1)
                FRP1 = self.BWIGML(S1, MRP, GRP, M2, M3, 1)
                FRO2 = self.BWIGML(S2, MRO, GRO, M3, M1, 1)
                FRP2 = self.BWIGML(S2, MRP, GRP, M3, M1, 1)
                FF21 = self.BWIGML(S1, MF2, GF2, M2, M3, 2)
                FF22 = self.BWIGML(S2, MF2, GF2, M3, M1, 2)
                FSG1 = self.BWIGML(S1, MSG, GSG, M2, M3, 0)
                FSG2 = self.BWIGML(S2, MSG, GSG, M3, M1, 0)
                FF01 = self.BWIGML(S1, MF0, GF0, M2, M3, 0)
                FF02 = self.BWIGML(S2, MF0, GF0, M3, M1, 0)

                F3PIFactor = BT3*(TComplex(F34A, 0.)*FRO1 + TComplex(F34B, 0.)*FRO2) \
                           + BT4*(TComplex(F34A, 0.)*FRP1 + TComplex(F34B, 0.)*FRP2) \
                           - BT5*(TComplex(F35A, 0.)*FF21 + TComplex(F35B, 0.)*FF22) \
                           - BT6*(TComplex(F36A, 0.)*FSG1 + TComplex(F36B, 0.)*FSG2) \
                           - BT7*(TComplex(F36A, 0.)*FF01 + TComplex(F36B, 0.)*FF02)

        if IDK == 1:
            if IFORM == 1 or IFORM == 2:
                S1 = SA 
                S2 = SB
                S3 = QQ - SA - SB + M1SQ + M2SQ + M3SQ

                # C it is 2pi0pi-
                # C Lorentz invariants for all the contributions:
                F134 = -(1 / 3.)  * ((S3 - M3SQ) - (S1 - M1SQ))                      # array
                F150 =  (1 / 18.) * (QQ - M3SQ + S3) * (2*M1SQ + 2*M2SQ - S3) / S3   # array
                F167 =  (2 / 3.)                                                     # scalar

                # FR**: all are TComplex
                FRO1 = self.BWIGML(S1, MRO, GRO, M2, M3, 1)
                FRP1 = self.BWIGML(S1, MRP, GRP, M2, M3, 1)
                FRO2 = self.BWIGML(S2, MRO, GRO, M3, M1, 1)
                FRP2 = self.BWIGML(S2, MRP, GRP, M3, M1, 1)
                FF23 = self.BWIGML(S3, MF2, GF2, M1, M2, 2)
                FSG3 = self.BWIGML(S3, MSG, GSG, M1, M2, 0)
                FF03 = self.BWIGML(S3, MF0, GF0, M1, M2, 0)

                
                F3PIFactor = BT1*FRO1 \
                           + BT2*FRP1 \
                           + BT3*TComplex(F134)*FRO2 \
                           + BT4*TComplex(F134)*FRP2 \
                           + BT5*TComplex(F150)*FF23 \
                           + BT6*TComplex(F167)*FSG3 \
                           + BT7*TComplex(F167)*FF03
                
            elif IFORM == 3:
                S3 = SA
                S1 = SB 
                S2 = QQ - SA - SB + M1SQ + M2SQ + M3SQ

                F34A = (1 / 3.) * ((S2 - M2SQ) - (S3 - M3SQ)) # array
                F34B = (1 / 3.) * ((S3 - M3SQ) - (S1 - M1SQ)) # array
                F35 = -(1 / 2.) * ((S1 - M1SQ) - (S2 - M2SQ)) # array

                FRO1 = self.BWIGML(S1, MRO, GRO, M2, M3, 1)
                FRP1 = self.BWIGML(S1, MRP, GRP, M2, M3, 1)
                FRO2 = self.BWIGML(S2, MRO, GRO, M3, M1, 1)
                FRP2 = self.BWIGML(S2, MRP, GRP, M3, M1, 1)
                FF23 = self.BWIGML(S3, MF2, GF2, M1, M2, 2)                    

                F3PIFactor = BT3*(TComplex(F34A)*FRO1 + TComplex(F34B)*FRO2) \
                           + BT4*(TComplex(F34A)*FRP1 + TComplex(F34B)*FRP2) \
                           + BT5*TComplex(F35)*FF23

        FORMA1 = self.FA1A1P(QQ) # TComplex
        out = F3PIFactor*FORMA1

        return  out



    # ------- L-wave BreightWigner for rho
    # Breit-Wigner function with isotropic decay angular distribution
    def BWIGML(self,
               S: ak.Array,
               M: float,
               G: float,
               m1: float,
               m2: float,
               L: int):
        MP = (m1 + m2)**2 # scalar
        MM = (m1 - m2)**2 # scalar
        MSQ = M**2        # scalar
        W = np.sqrt(S)    # array

        # Prepare a WGS array
        # use ak.where 
        wgs = self.GetWGS(S, MP, MM, MSQ, L, G, W, M)
        dummy = ak.zeros_like(wgs)
        mask = (W > m1+m2)
        WGS = ak.where(mask, wgs, dummy)

        num = TComplex(MSQ, 0)
        den = TComplex((MSQ - S), -WGS)
        
        out = num/den

        return out

    
    def GetWGS(self, S: ak.Array, MP: float, MM: float, MSQ: float, L: int,
               G: float, W: ak.Array, M: float):

        QS = np.sqrt(np.abs((S - MP) * (S - MM))) / W
        QM = np.sqrt(np.abs((MSQ - MP) * (MSQ - MM))) / M
        IPOW = 2 * L + 1
        WGS = G * (MSQ / W) * np.power((QS / QM), IPOW)

        return WGS

    
    def FA1A1P(self, XMSQ: ak.Array) -> TComplex:
        XM1 = 1.275000
        XG1 = 0.700
        XM2 = 1.461000
        XG2 = 0.250
        BET = TComplex(0.0, 0.0)

        GG1 = XM1*XG1/(1.3281*0.806)
        GG2 = XM2*XG2/(1.3281*0.806)
        XM1SQ = XM1*XM1
        XM2SQ = XM2*XM2

        GF = self.WGA1(XMSQ) # array
        FG1 = GG1*GF
        FG2 = GG2*GF
          
        F1 = TComplex(-XM1SQ)/TComplex(XMSQ - XM1SQ, FG1)
        F2 = TComplex(-XM2SQ)/TComplex(XMSQ - XM2SQ, FG2)
        FA1A1P = F1 + (BET*F2)

        return FA1A1P


    def WGA1(self, QQ: ak.Array):
        # C mass-dependent M*Gamma of a1 through its decays to
        # C.[(rho-pi S-wave) + (rho-pi D-wave) +
        # C.(f2 pi D-wave) + (f0pi S-wave)]
        # C.AND simple K*K S-wave
        MKST = 0.894
        MK   = 0.496
        MK1SQ = (MKST+MK)**2
        MK2SQ = (MKST-MK)**2
        # C coupling constants squared:
        C3PI = 0.2384*0.2384
        CKST = 4.7621*4.7621*C3PI
        # C Parameterization of numerical integral of total width of a1 to 3pi.
        # C From M. Schmidtler, CBX-97-64-Update.

        S = QQ
        WG3PIC = self.WGA1C(S)
        WG3PIN = self.WGA1N(S)
            
        # C Contribution to M*Gamma(m(3pi)^2) from S-wave K*K, if above threshold
        #GKST = 0.0
        mask = S > MK1SQ
        gskt = np.sqrt((S-MK1SQ)*(S-MK2SQ))/(2.0*S)
        dummy = ak.zeros_like(gskt)
        GKST = ak.where(mask, gskt, dummy)

        out = C3PI*(WG3PIC+WG3PIN) + (CKST*GKST)

        return out



    def WGA1C(self, S: ak.Array):
        STH = 0.1753 
        Q0  = 5.80900
        Q1  = -3.00980 
        Q2  = 4.57920 
        P0  = -13.91400 
        P1  = 27.67900 
        P2  = -13.39300 
        P3  = 3.19240 
        P4  = -0.10487

        mask1 = S < STH
        mask2 = (S > STH) & (S < 0.823)

        dummy = ak.zeros_like(S)
        ifmask2 = Q0 * ((S - STH)*(S - STH)*(S - STH)) * (1.0 + Q1 * (S - STH) + Q2 * (S - STH)*(S - STH))
        elsemask2 = P0 + P1*S + P2*S*S + P3*S*S*S + P4*S*S*S*S
        
        G1_IM = ak.where(
            mask1,
            dummy,
            ak.where(
                mask2, ifmask2, elsemask2
            )
        )

        return G1_IM

    
    def WGA1N(self, S: ak.Array):
        Q0 = 6.28450
        Q1 = -2.95950
        Q2 = 4.33550
        P0 = -15.41100
        P1 = 32.08800
        P2 = -17.66600
        P3 = 4.93550
        P4 = -0.37498
        STH = 0.1676

        mask1 = S < STH
        mask2 = (S > STH) & (S < 0.823)

        dummy = ak.zeros_like(S)
        ifmask2 = Q0 * ((S - STH)*(S - STH)*(S - STH)) * (1.0 + Q1 * (S - STH) + Q2 * (S - STH)*(S - STH))
        elsemask2 = P0 + P1*S + P2*S*S + P3*S*S*S + P4*S*S*S*S
        
        G1_IM = ak.where(
            mask1,
            dummy,
            ak.where(
                mask2, ifmask2, elsemask2
            )
        )

        return G1_IM

        
    def CLVEC(self, H: list, HC: list, N: ak.Array) -> ak.Array:
        HN  = H[0]*N.energy - H[1]*N.px - H[2]*N.py - H[3]*N.pz        # TComplex
        HCN = HC[0]*N.energy - HC[1]*N.px - HC[2]*N.py - HC[3]*N.pz    # TComplex
        HH  = (H[0]*HC[0] - H[1]*HC[1] - H[2]*HC[2] - H[3]*HC[3]).Re() # ak.Array

        PIVEC0 = 2*( 2*(HN*HC[0]).Re() - HH*N.energy )
        PIVEC1 = 2*( 2*(HN*HC[1]).Re() - HH*N.px     )
        PIVEC2 = 2*( 2*(HN*HC[2]).Re() - HH*N.py     )
        PIVEC3 = 2*( 2*(HN*HC[3]).Re() - HH*N.pz     )

        out = ak.zip(
            {
                "x": PIVEC1,
                "y": PIVEC2,
                "z": PIVEC3,
                "t": PIVEC0,
            },
            with_name="LorentzVector",
            behavior=coffea.nanoevents.methods.vector.behavior,
        )

        return out
     
    
    def CLAXI(self, H: list, HC: list, N: ak.Array) -> ak.Array:
        a1 = HC[1]
        a2 = HC[2]
        a3 = HC[3]
        a4 = HC[0]

        b1 = H[1]
        b2 = H[2]
        b3 = H[3]
        b4 = H[0]

        c1 = N.px
        c2 = N.py
        c3 = N.pz
        c4 = N.energy   
        
        d34 = (a3*b4 - a4*b3).Im()
        d24 = (a2*b4 - a4*b2).Im()
        d23 = (a2*b3 - a3*b2).Im()
        d14 = (a1*b4 - a4*b1).Im()
        d13 = (a1*b3 - a3*b1).Im()
        d12 = (a1*b2 - a2*b1).Im()

        PIAX0 = -self.SIGN*2*(-c1*d23 + c2*d13 - c3*d12)
        PIAX1 = self.SIGN*2*(c2*d34 - c3*d24 + c4*d23)
        PIAX2 = self.SIGN*2*(-c1*d34 + c3*d14 - c4*d13)
        PIAX3 = self.SIGN*2*(c1*d24 - c2*d14 + c4*d12)
        
        out = ak.zip(
            {
                "x": PIAX1,
                "y": PIAX2,
                "z": PIAX3,
                "t": PIAX0,
            },
            with_name="LorentzVector",
            behavior=coffea.nanoevents.methods.vector.behavior,
        )

        return out
