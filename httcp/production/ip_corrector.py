import functools
import law
from columnflow.production import Producer, producer
from columnflow.util import maybe_import, safe_div
from law.util import InsertableDict
from columnflow.columnar_util import set_ak_column, EMPTY_FLOAT, flat_np_view, optional_column as optional
from columnflow.production.util import attach_coffea_behavior

np = maybe_import("numpy")
ak = maybe_import("awkward")


@producer(
    uses={
        'event', 'hcand_*', optional('gen_lep*'),
    },
    produces={
        'hcand_*'
    },
)
def ip_correction(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    ch_str     = self.config_inst.channels.names()[0]
    ch_objects  = self.config_inst.x.ch_objects[ch_str]
    if ch_str == 'mutau':
        ip_correction = self.ip_corr_mm
    elif ch_str == 'etau':
        ip_correction = self.ip_corr_ee
    else:
        raise NotImplementedError(f"IP correction is not implemented for {ch_str} channel")
    
    tau_part_flav = {
                        "prompt_e"  : 1,
                        "prompt_mu" : 2,
                        "tau->e"    : 3,
                        "tau->mu"   : 4,
                        "tau_had"   : 5
                    }
    
    hcand = events[f'hcand_{ch_str}']
    output_hcand = {}
    
    for lep_str in ['lep0','lep1']:
        lep = hcand[lep_str]
        if not self.dataset_inst.is_mc: #For data duplicate values of IP variable to IP_qm
            for the_comp in ['x','y','z']:
                #print(f'Do not apply IP correction: Filling IP_{the_comp}_qm with  the values from IP_{the_comp}')
                lep[f'IP{the_comp}_qm'] = lep[f'IP{the_comp}']
        else:
            gen_lep = events.gen_lep[lep_str]
            abseta = flat_np_view(abs(lep.eta), axis=1)
            max_eta = 2.4 if ch_str=='mutau' else 2.1
            mask = (abseta <= max_eta) 
            genmatch = flat_np_view(lep.genPartFlav, axis=1)
            if ch_objects[lep_str] == 'Tau':
                if ch_str == 'mutau': 
                    mask = mask & ((genmatch == tau_part_flav["prompt_mu"]) | (genmatch == tau_part_flav["tau->mu"]) | (genmatch == tau_part_flav["tau_had"])) 
                elif ch_str == 'etau': 
                    mask = mask & ((genmatch == tau_part_flav["prompt_e"]) | (genmatch == tau_part_flav["tau->e"]) | (genmatch == tau_part_flav["tau_had"])) 
                else:
                    pass
            elif ch_objects[lep_str] == 'Muon':
                mask = mask & ((genmatch == tau_part_flav["prompt_mu"]) | genmatch == 15)
            for the_comp in ['x','y','z']:
                #print(f'performing IP correction for {the_comp} component...')
                gen_ip_comp = flat_np_view(gen_lep[f'gen_ip_{the_comp}'], axis=1)
                ip_comp     = flat_np_view(lep[f'IP{the_comp}'], axis=1)
                #DIRTY HACK PLEASE REPLACE IT ASAP
                if lep_str == 'lep1':
                    gen2reco_shift = ip_comp 
                else:
                    gen2reco_shift = (ip_comp - gen_ip_comp)
                mask = mask & (abs(gen2reco_shift) < 0.03) #Input for the correction should have IP component magnitude < 0.03
                corrected_ip = ip_comp.copy()
                corrected_ip[mask] = ip_correction.evaluate(gen2reco_shift[mask],
                                                            the_comp, 
                                                            abseta[mask])
                #DIRTY HACK PLEASE REPLACE IT ASAP
                if lep_str == 'lep1':
                    pass
                else:
                    corrected_ip[mask] = corrected_ip[mask] + gen_ip_comp[mask]
               
                lep[f'IP{the_comp}_qm'] = ak.unflatten(corrected_ip, ak.num(lep.pt, axis=1)) 
        output_hcand[lep_str] = lep 
    events = set_ak_column(events, f'hcand_{ch_str}', ak.zip(output_hcand))     
    return events

@ip_correction.requires
def ip_correction_requires(self: Producer, task: law.Task, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(task)

@ip_correction.setup
def ip_correction_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
    task: law.task
) -> None:
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    
    correction_set = correctionlib.CorrectionSet.from_file(bundle.files.ip_corr.path)
    self.ip_corr_ee = correction_set['ip_correction_ee']
    self.ip_corr_mm = correction_set['ip_correction_mm']
    
    
    
    


@producer(
    uses={
        'event', 'hcand_*', optional('gen_lep*'),
    },
    produces={
        'hcand_*'
    },
)
def ip_sig_correction(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    ch_str     = self.config_inst.channels.names()[0]
    ch_objects  = self.config_inst.x.ch_objects[ch_str]
    if ch_str == 'mutau':
        ip_correction = self.ip_corr_mm
    elif ch_str == 'etau':
        ip_correction = self.ip_corr_ee
    else:
        raise NotImplementedError(f"IP correction is not implemented for {ch_str} channel")
    
    tau_part_flav = {
                        "prompt_e"  : 1,
                        "prompt_mu" : 2,
                        "tau->e"    : 3,
                        "tau->mu"   : 4,
                        "tau_had"   : 5
                    }
    
    hcand = events[f'hcand_{ch_str}']
    output_hcand = {}
    
    for lep_str in ['lep0','lep1']:
        lep = hcand[lep_str]
        if not self.dataset_inst.is_mc: #For data duplicate values of IP variable to IP_qm
            for the_comp in ['x','y','z']:
                #print(f'Do not apply IP correction: Filling IP_{the_comp}_qm with  the values from IP_{the_comp}')
                lep[f'IP{the_comp}_qm'] = lep[f'IP{the_comp}']
        else:
            gen_lep = events.gen_lep[lep_str]
            abseta = flat_np_view(abs(lep.eta), axis=1)
            max_eta = 2.4 if ch_str=='mutau' else 2.1
            mask = (abseta <= max_eta) 
            genmatch = flat_np_view(lep.genPartFlav, axis=1)
            if ch_objects[lep_str] == 'Tau':
                if ch_str == 'mutau': 
                    mask = mask & ((genmatch == tau_part_flav["prompt_mu"]) | (genmatch == tau_part_flav["tau->mu"]) | (genmatch == tau_part_flav["tau_had"])) 
                elif ch_str == 'etau': 
                    mask = mask & ((genmatch == tau_part_flav["prompt_e"]) | (genmatch == tau_part_flav["tau->e"]) | (genmatch == tau_part_flav["tau_had"])) 
                else:
                    pass
            elif ch_objects[lep_str] == 'Muon':
                mask = mask & ((genmatch == tau_part_flav["prompt_mu"]) | genmatch == 15)
            for the_comp in ['x','y','z']:
                #print(f'performing IP correction for {the_comp} component...')
                gen_ip_comp = flat_np_view(gen_lep[f'gen_ip_{the_comp}'], axis=1)
                ip_comp     = flat_np_view(lep[f'IP{the_comp}'], axis=1)
                #DIRTY HACK PLEASE REPLACE IT ASAP
                if lep_str == 'lep1':
                    gen2reco_shift = ip_comp 
                else:
                    gen2reco_shift = (ip_comp - gen_ip_comp)
                mask = mask & (abs(gen2reco_shift) < 0.03) #Input for the correction should have IP component magnitude < 0.03
                corrected_ip = ip_comp.copy()
                corrected_ip[mask] = ip_correction.evaluate(gen2reco_shift[mask],
                                                            the_comp, 
                                                            abseta[mask])
                #DIRTY HACK PLEASE REPLACE IT ASAP
                if lep_str == 'lep1':
                    pass
                else:
                    corrected_ip[mask] = corrected_ip[mask] + gen_ip_comp[mask]
               
                lep[f'IP{the_comp}_qm'] = ak.unflatten(corrected_ip, ak.num(lep.pt, axis=1)) 
        output_hcand[lep_str] = lep 
    events = set_ak_column(events, f'hcand_{ch_str}', ak.zip(output_hcand))     
    return events

@ip_sig_correction.requires
def ip_correction_requires(self: Producer, task: law.Task, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(task)

@ip_sig_correction.setup
def ip_correction_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
    task: law.task
) -> None:
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    
    correction_set = correctionlib.CorrectionSet.from_file(bundle.files.ip_corr.path)
    self.ip_corr_ee = correction_set['ip_correction_ee']
    self.ip_corr_mm = correction_set['ip_correction_mm']