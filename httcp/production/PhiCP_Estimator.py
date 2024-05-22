import os
from columnflow.util import maybe_import

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")

from httcp.production.PolarimetricA1 import PolarimetricA1
from IPython import embed


def GetPhiCP(
        p4hcandinfodict: dict, 
        method_leg1: str, 
        method_leg2: str,
        mode_leg1: str, 
        mode_leg2: str,
        **kwars
) -> ak.Array:
    #print(f"{method_leg1}-{method_leg2}")
    #print(f"{mode_leg1}-{mode_leg2}")
    """
    Inputs:
      p4hcandinfodict : full hcand dict
      method_leg1     : "DP/PV/IP"
      method_leg1     : "DP/PV/IP"
      mode_leg1       : "e/mu/pi/rho/a1"
      mode_leg1       : "pi/rho/a1"
    Output:
      Final PhiCP array for any methods / configurations
    Steps:
      Prepares the input vectors for PhiCP
      Compute the PhiCP
    """
    phicp_input_vec_dict = PrepareVecsForPhiCP(p4hcandinfodict,
                                               method_leg1, method_leg2, 
                                               mode_leg1, mode_leg2)
    phicp = ak.values_astype(ComputeAcopAngle(phicp_input_vec_dict), np.float32)
    phicp = ak.enforce_type(phicp, "var * float32")
    
    return phicp


def PrepareVecsForPhiCP(
        p4hcandinfodict: dict, 
        method_leg1: str, 
        method_leg2: str, 
        mode_leg1: str, 
        mode_leg2: str
) -> dict:
    """
    Inputs:
      p4hcandinfodict : full hcand dict
      method_leg1     : "DP/PV/IP"
      method_leg1     : "DP/PV/IP"
      mode_leg1       : "e/mu/pi/rho/a1"
      mode_leg1       : "pi/rho/a1"
    Output:
      dict{
        P1 : tau/pi/... for hcand1
        R1 : k/pi0/...  for hcand1  [to get phicp]
        H1 : h/pi0/...  for hcand1  [to calculate angle]
        C1 : charge of hcand1
        P2 : tau/pi/... for hcand2
        R2 : h/pi0/...  for hcand2  [to get phicp]
        H2 : k/pi0/...  for hcand2  [to calculate angle]
        C2 : charge of hcand2
        Y1 : Phase-shift for hacnd1 [only for DP]
        Y2 : Phase-shift for hacnd2 [only for DP]
      }
    Steps:
      Prepare input vectors for PhiCP calculation
    """
    p4h1    = p4hcandinfodict["p4h1"]
    p4h1pi  = p4hcandinfodict["p4h1pi"]
    p4h1pi0 = p4hcandinfodict["p4h1pi0"]
    p4h2    = p4hcandinfodict["p4h2"]
    p4h2pi  = p4hcandinfodict["p4h2pi"]
    p4h2pi0 = p4hcandinfodict["p4h2pi0"]
    
    _P1, _R1, _C1 = _prepareVecs(p4h1, p4h1pi, p4h1pi0, method_leg1, mode_leg1)
    _P2, _R2, _C2 = _prepareVecs(p4h2, p4h2pi, p4h2pi0, method_leg2, mode_leg2)

    # Get Boost
    boostv = _getBoost(_P1, _R1,
                       _P2, _R2,
                       method_leg1, 
                       method_leg2,
                       mode_leg1, 
                       mode_leg2)

    P1, R1, H1, Y1 = _reStructureVecs(boostv, _P1, _R1, p4h1pi0, method_leg1, mode_leg1)
    P2, R2, H2, Y2 = _reStructureVecs(boostv, _P2, _R2, p4h2pi0, method_leg2, mode_leg2)

    return {"P1" : P1, "R1" : R1, "H1": H1, "C1": _C1,
            "P2" : P2, "R2" : R2, "H2": H2, "C2": _C2,
            "Y1" : Y1, "Y2" : Y2}



def ComputeAcopAngle(vecsdict):
    """
    Geometrical estimation of PhiCP
    For DPDP: 
       P1/P2  --> vecPiPlus
       R1/R2  --> vecPiZeroPlustransv
       H1/H2  --> R1/R2
       Y1     --> phase-shift for leg 1
       Y2     --> phase-shift for leg 2
    """
    #embed()
    assert len(vecsdict) == 10, "input dict to ComputeAcopAngle does not have proper structure"
    P1 = vecsdict["P1"]
    R1 = vecsdict["R1"]
    H1 = vecsdict["H1"]
    C1 = vecsdict["C1"]
    P2 = vecsdict["P2"]
    R2 = vecsdict["R2"]
    H2 = vecsdict["H2"]
    C2 = vecsdict["C2"]
    Y1 = vecsdict["Y1"]
    Y2 = vecsdict["Y2"]
    Y  = None
    if Y1 is not None:
        if Y2 is not None:
            Y = Y1*Y2
        else:
            Y = Y1
    else:
        if Y2 is not None:
            Y = Y2
            
    acop = np.arccos(R2.dot(R1))
    sign = _getSign(P1,H1,P2,H2,C1,C2) #P1.dot(H2.cross(H1))
    acop = ak.where(sign < 0.0, 2*np.pi - acop, acop)
    if Y is not None:
        acop  = ak.where(Y < 0.0, acop + np.pi, acop)
        acop  = ak.where((Y < 0.0) & (acop > 2*np.pi), 
                         acop - 2*np.pi, 
                         acop)
        
    return acop


# PRIVATE
def _getSign(P1,H1,P2,H2,C1,C2):
    #embed()
    Pm = ak.where(C1 < 0, P1, P2) # check tau-
    Hm = ak.where(C1 < 0, H1, H2) # Sort according to tau charge
    Hp = ak.where(C1 < 0, H2, H1) # Same
    angle = Pm.dot(Hp.cross(Hm))
    return angle
    

def _getBoost(
        P1   : ak.Array, 
        R1   : ak.Array, 
        P2   : ak.Array, 
        R2   : ak.Array, 
        leg1_method: str, 
        leg2_method: str, 
        leg1_mode: str, 
        leg2_mode: str,
        **kwargs
) -> ak.Array:

    frame = None
    if leg1_method == "PV" and leg2_method == "PV":
        # P1/P2 -- tau1/2, R1/R2 -- pion_tau1/2
        if leg1_mode == "pi" and leg2_mode == "pi":
            frame = R1 + R2
        else:
            frame = P1 + P2
    else:
        frame = P1 + P2
    """
    elif leg1_method == "DP" and leg2_method == "DP":
        #v1p -- pi+
        #v2p -- pi0 from tau+
        #v1m -- pi-
        #v2m -- pi0 from tau-
        #frame -- pi+ + pi-
        frame = v1p+ v1m
    elif leg1_method == "IP" and leg2_method == "IP":
        # v1p -- pion+
        # v2p -- eta+
        # v1m -- pion-
        # v2m -- eta-
        frame = v1p+ v1m
    elif leg1_method == "PV" and leg2_method == "IP":
        # v1p -- tau+
        # v2p -- pion+
        # v1m -- pion-
        # v2m -- eta-
        frame = v1p + v1m
    elif leg1_method == "IP" and leg2_method == "PV":
        # v1p -- pion+
        # v2p -- eta+
        # v1m -- tau-
        # v2m -- pion-
        frame = v1p + v1m
    elif leg1_method == "DP" and leg2_method == "IP":
        # v1p -- pion+
        # v2p -- pi0+
        # v1m -- pion-
        # v2m -- eta-
        frame = v1p  + v1m
    elif leg1_method == "IP" and leg2_method == "DP":
        # v1p -- pion+
        # v2p -- eta+
        # v1m -- pion-
        # v2m -- pi0-
        frame = v1p + v1m
    """
    return frame.boostvec


def _prepareVecs(
        hcand: ak.Array, 
        hcand_pi: ak.Array, 
        hcand_pi0: ak.Array, 
        leg_method: str, 
        leg_mode: str
) -> tuple[ak.Array, ak.Array, ak.Array]:
    """
    A private function to be used in PrepareVecsForPhiCP
    This mainly returns the required vectors for different modes and methods
    Returned vectors will further be modified in the reStructure step
    Method : PV
      Mode : Pi, Rho, a1 [only possible modes]
        WARNING: Frame will be different for Pi
                 and it will use R instead of P
        P --> tau (hcand)
        R --> pi  (pions)
        It will be converted to polarimertric vector in next function
    Method : DP
      Mode : Rho, a1 [only possible modes]
        P --> pi  (diff for a1)
        R --> pi0 (diff for a1)
    Method : IP
      Mode : e, mu, pi [only needed]
        P -->
        R -->
      
    """
    P = None
    R = None

    if leg_method == "DP":
        if leg_mode == "rho":
            P = hcand_pi
            R = hcand_pi0
        elif leg_mode == "a1":
            P = _get_pi_a1_DP(hcand_pi)
            R = hcand_pi[:,:1]
        else:
            raise RuntimeError(f"Wrong mode : {leg_mode}")


    elif leg_method == "PV":
        P = hcand
        R = hcand_pi

    elif leg_method == "IP":
        pass

    else:
        raise RuntimeError(f"Wrong {leg_method}")

    #embed()
    hcand_charge = hcand.charge
    return P, R, hcand_charge


def _reStructureVecs(
        boostv: ak.Array, 
        V1: ak.Array, 
        V2: ak.Array, 
        V3: ak.Array, 
        leg_method: str, 
        leg_mode: str
) -> tuple[ak.Array, ak.Array, ak.Array, ak.Array]:
    """
    PostProcess the inputs to phicp with boost
    Output : 
    Method : PV
      Mode : Pi, Rho, a1 [only possible modes]
        WARNING: Frame will be different for Pi
        P --> tau (hcand)
        R --> k   (tau x polarimetric vector)
        H --> h   (pv)
    Method : DP
      Mode : Rho, a1 [only possible modes]
        P --> pi  (diff for a1)
        R --> pi0 (diff for a1)
        H --> pi0 (same as R)
    Method : IP
      Mode : e, mu, pi [only needed]
        P -->
        R -->
        H -->
    """
    P = None
    R = None
    H = None # to calculate angle
    Y = None
    
    
    if leg_method == "DP":
        """
        Input
          V1   : Pi   [from _prepareVecs]
          V2   : Pi0  [from _prepareVecs]
          V3   : hcand_pi0 [not required]
        Output
          P : Pion at zero momentum frame
          R : Pi0 (Transverse) at zero momentum frame
          H : R [same, as it will be used to compute the sign]
          Y : Phase shift
        """
        pi_ZMF  = V1.boost(boostv.negative())
        #embed()

        pi0_ZMF = V2.boost(boostv.negative())
        pi_ZMF_unit  = pi_ZMF.pvec.unit
        pi0_ZMF_unit = pi0_ZMF.pvec.unit
        pi0_ZMF_unit_T = (pi0_ZMF_unit - pi_ZMF_unit*(pi_ZMF_unit.dot(pi0_ZMF_unit))).unit
    
        P = pi_ZMF_unit
        R = pi0_ZMF_unit_T
        H = R
        Y = (pi_ZMF.energy - pi0_ZMF.energy)/(pi_ZMF.energy + pi0_ZMF.energy)
    
    elif leg_method == "PV":
        """
        Input
          V1   : tau [i.e. hcand]
          V2   : Pi  [hcand decay]
          V3   : Pi0 [hcand decay]
        """

        P, R, H = _pv(boostv, V1, V2, V3, leg_mode)

    elif leg_method == "IP":
        pass


    else:
        raise RuntimeError(f"Wrong method: {leg_method}")

    return P, R, H, Y


def _get_pi_a1_DP(p4_pi):
    Minv1 = (p4_pi[:,:1] + p4_pi[:,1:2]).mass
    Minv2 = (p4_pi[:,:1] + p4_pi[:,2:3]).mass
    Pi = ak.where(np.abs(0.77526-Minv1) < np.abs(0.77526-Minv2), p4_pi[:,1:2], p4_pi[:,2:3])
    return Pi


def _pv(boostv: ak.Array, 
        tau : ak.Array, 
        pi  : ak.Array, 
        pi0 : ak.Array, 
        leg_mode: str) -> ak.Array:
    P  = tau.boost(boostv.negative())
    if leg_mode == "pi":
        pv = pi.boost(boostv.negative()).pvec
    elif leg_mode == "rho":
        pi  = pi.boost(boostv.negative())
        pi0 = pi0.boost(boostv.negative())
        q   = pi.subtract(pi0)
        N   = P.subtract(pi.add(pi0))
        pv  = (((2*(q.dot(N))*q.pvec).subtract(q.mass2*N.pvec)))
    elif leg_mode == "a1":
        os_pi_HRF   = pi[:, 0:1].boost(boostv.negative())
        ss1_pi_HRF  = pi[:, 1:2].boost(boostv.negative()) 
        ss2_pi_HRF  = pi[:, 2:3].boost(boostv.negative()) 
        a1pol       = PolarimetricA1(P,
                                     os_pi_HRF,
                                     ss1_pi_HRF,
                                     ss2_pi_HRF,
                                     tau.charge)
        pv          = -a1pol.PVC().pvec
    else:
        raise RuntimeError(f"Wrong mode: {leg_mode}")
    
    P = P.pvec.unit
    H = pv.unit
    R = (H.cross(P)).unit

    return P, R, H
