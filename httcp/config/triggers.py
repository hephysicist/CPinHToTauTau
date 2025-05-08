# coding: utf-8

"""
Definition of triggers
"""

import order as od

from httcp.config.trigger_util import Trigger, TriggerLeg
    
def add_triggers_run3(config: od.Config) -> None:
    """
    Adds all triggers to a *config*. For the conversion from filter names to trigger bits, see
    https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py.
    """
    config.x.triggers = od.UniqueObjectIndex(Trigger,[
        
        #
        # single muon
        #
        # https://github.com/cms-sw/cmssw/blob/203834e3ae301f2564423dd1cc84bebf660519b9/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L118
        Trigger(
            name="HLT_IsoMu24",
            id=132, #13 is for muon pdg_id, 1 because it's first muon trigger
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    trigger_bits=[2,4],
                    min_pt=26.0,
                    min_eta=2.4,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        # #
        # # single electron
        # #    
        # # https://github.com/cms-sw/cmssw/blob/203834e3ae301f2564423dd1cc84bebf660519b9/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L63
        # Trigger(
        # name="HLT_Ele30_WPTight_Gsf",
        # id=111,
        # legs=[
        #     TriggerLeg(
        #         pdg_id=11,
        #         min_pt=31.0,
        #         min_eta=2.4,
        #         # filter names:
        #         # hltEle32WPTightGsfTrackIsoFilter
        #         trigger_bits=[2],
        #     ),
        # ],
        # tags={"single_trigger", "single_e", "channel_e_tau"},
        # ),
        # mu tauh
        #
        # Muon : https://github.com/cms-sw/cmssw/blob/203834e3ae301f2564423dd1cc84bebf660519b9/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L121
        # Tau  : https://github.com/cms-sw/cmssw/blob/203834e3ae301f2564423dd1cc84bebf660519b9/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L149
        Trigger(
            name="HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1",
            id=13,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=21.0,
                    min_eta=2.1,
                    # filter names:
                    # hltL3crIsoBigORMu18erTauXXer2p1L1f0L2f10QL3f20QL3trkIsoFiltered  # TODO Twiki sugests 2
                    # hltHpsOverlapFilterIsoMu20LooseMuTauWPDeepTauPFTau27L1Seeded (OverlapFilter PFTau) # TODO Twiki sugests 4 + 64  # noqa
                    trigger_bits= [2,4],#[2,3,7], #4,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=32.0,
                    min_eta=2.1,
                    # filter names:
                    # (DeepTau + HPS) # TODO Twiki sugests 8 + 32 + 512 + 262144
                    # hltHpsOverlapFilterIsoMu20LooseMuTauWPDeepTauPFTau27L1Seeded
                    trigger_bits= [4,10],#[4,14,29],
                ),
            ],
            tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        ),
     ])

def add_triggers_run3_2022_preEE(config: od.Config) -> None:
    """
    Adds all triggers to a *config*. For the conversion from filter names to trigger bits, see
    https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py.
    """
    config.x.triggers = od.UniqueObjectIndex(Trigger,[
        
        #
        # single muon
        #
        # https://github.com/cms-sw/cmssw/blob/203834e3ae301f2564423dd1cc84bebf660519b9/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L118
        Trigger(
            name="HLT_IsoMu24",
            id=132, #13 is for muon pdg_id, 1 because it's first muon trigger
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    trigger_bits=[2,4],
                    min_pt=25.0,
                    min_eta=2.4,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        
        # single electron
        #    
        # https://github.com/cms-sw/cmssw/blob/203834e3ae301f2564423dd1cc84bebf660519b9/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L63
        Trigger(
        name="HLT_Ele30_WPTight_Gsf",
        id=111,
        legs=[
            TriggerLeg(
                pdg_id=11,
                min_pt=31.0,
                min_eta=2.4,
                # filter names:
                # hltEle32WPTightGsfTrackIsoFilter
                trigger_bits=[2],
            ),
        ],
        tags={"single_trigger", "single_e", "channel_e_tau"},
        ),

        #
        # e tauh
        #
        # Electron leg: https://github.com/cms-sw/cmssw/blob/203834e3ae301f2564423dd1cc84bebf660519b9/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L69
        # Tau leg: https://github.com/cms-sw/cmssw/blob/203834e3ae301f2564423dd1cc84bebf660519b9/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L148
        # Trigger(
        #     name="HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1",
        #     id=15,
        #     legs=[
        #         TriggerLeg(
        #             pdg_id=11,
        #             min_pt=27.0,
        #             min_eta=2.1,
        #             # filter names:
        #             # hltEle24erWPTightGsfTrackIsoFilterForTau
        #             # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
        #             trigger_bits= [2], #[2,4,8],
        #         ),
        #         TriggerLeg(
        #             pdg_id=15,
        #             min_pt=35.0,
        #             min_eta=2.1,
        #             # filter names:
        #             # hltSelectedPFTau30LooseChargedIsolationL1HLTMatched
        #             # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
        #             trigger_bits= [4,9],#[4,13,29],
        #         ),
        #     ],
        #     tags={"cross_trigger", "cross_e_tau", "channel_e_tau"},
        # ),
        
        # # mu tauh
        # #
        # # Muon : https://github.com/cms-sw/cmssw/blob/203834e3ae301f2564423dd1cc84bebf660519b9/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L121
        # # Tau  : https://github.com/cms-sw/cmssw/blob/203834e3ae301f2564423dd1cc84bebf660519b9/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L149
        # Trigger(
        #     name="HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1",
        #     id=13,
        #     legs=[
        #         TriggerLeg(
        #             pdg_id=13,
        #             min_pt=21.0,
        #             min_eta=2.1,
        #             # filter names:
        #             # hltL3crIsoBigORMu18erTauXXer2p1L1f0L2f10QL3f20QL3trkIsoFiltered  # TODO Twiki sugests 2
        #             # hltHpsOverlapFilterIsoMu20LooseMuTauWPDeepTauPFTau27L1Seeded (OverlapFilter PFTau) # TODO Twiki sugests 4 + 64  # noqa
        #             trigger_bits= [2,4],#[2,3,7], #4,
        #         ),
        #         TriggerLeg(
        #             pdg_id=15,
        #             min_pt=32.0,
        #             min_eta=2.1,
        #             # filter names:
        #             # (DeepTau + HPS) # TODO Twiki sugests 8 + 32 + 512 + 262144
        #             # hltHpsOverlapFilterIsoMu20LooseMuTauWPDeepTauPFTau27L1Seeded
        #             trigger_bits= [4,10],#[4,14,29],
        #         ),
        #     ],
        #     tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        # ),
        # # #
        # # # tauh tauh
        # # #
        # # # https://github.com/cms-sw/cmssw/blob/203834e3ae301f2564423dd1cc84bebf660519b9/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L147 
        # Trigger(
        #     name="HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1",
        #     id=5,
        #     legs=[
        #         TriggerLeg(
        #             pdg_id=15,
        #             min_pt=40.0,
        #             min_eta=2.1,
        #             # filter names:
        #             # hltHpsDoublePFTau35MediumDitauWPDeepTauDz02 (Deeptau + HPS)
        #             trigger_bits= [4,8], # [4,12,29],
        #         ),
        #         TriggerLeg(
        #             pdg_id=15,
        #             min_pt=40.0,
        #             min_eta=2.1,
        #             # filter names:
        #             # hltHpsDoublePFTau35MediumDitauWPDeepTauDz02 (Deeptau + HPS)
        #             trigger_bits= [4,8], #[4,12,29], 
        #         ),
        #     ],
        #     tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        # ),
     ])

