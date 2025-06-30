from HiggsAnalysis.CombinedLimit.PhysicsModel import PhysicsModel
from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder
import numpy as np

class CPMixtureDecays(PhysicsModel):
    def __init__(self):
        PhysicsModel.__init__(self)
        self.do_kappas = False

    def setPhysicsOptions(self, physOptions):
        for po in physOptions:
            if po.startswith("do_kappas"):
                self.do_kappas = True

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        # --- POI and other parameters ----

        poiNames = []

        params = {
                'pi': np.pi
        }
 
        if not self.do_kappas: 
          self.modelBuilder.doVar('alpha[0,-90,90]')
          poiNames.append('alpha')
          self.modelBuilder.doVar('mutautau[1]')
          self.modelBuilder.doVar('muV[1,0,10]')
          self.modelBuilder.doVar('muggH[1,0,10]')

          self.modelBuilder.doSet('POI', ','.join(poiNames))

          self.modelBuilder.factory_('expr::a1("sqrt(@0)*cos(({pi}/2)*(@1/90))", mutautau, alpha)'.format(**params))
          self.modelBuilder.factory_('expr::a3("sqrt(@0)*sin(({pi}/2)*(@1/90))", mutautau, alpha)'.format(**params))

          self.modelBuilder.factory_('expr::sm_scaling("@0*@0 - @0*@1", a1, a3)')
          self.modelBuilder.factory_('expr::ps_scaling("@1*@1 - @0*@1", a1, a3)')
          self.modelBuilder.factory_('expr::mm_scaling("2*@0*@1", a1, a3)')

        else: 
          self.modelBuilder.doVar('kappaH[1,-2,2]')
          self.modelBuilder.doVar('kappaA[0,-2,2]')
          poiNames.append('kappaH')
          poiNames.append('kappaA')
          self.modelBuilder.doSet('POI', ','.join(poiNames))
          self.modelBuilder.doVar('muV[1]')
          self.modelBuilder.doVar('muggH[1]')
          self.modelBuilder.doVar('mutautau[1]')

          self.modelBuilder.factory_('expr::a1("@0", kappaH)'.format(**params))
          self.modelBuilder.factory_('expr::a3("@0", kappaA)'.format(**params))

          self.modelBuilder.factory_('expr::SM_width_correction("1/(1 + 6.272E-02*(@0*@0+@1*@1 -1))", kappaA, kappaH)')

          self.modelBuilder.factory_('expr::sm_scaling("(@0*@0 - @0*@1)*@2", a1, a3, SM_width_correction)')
          self.modelBuilder.factory_('expr::ps_scaling("(@1*@1 - @0*@1)*@2", a1, a3, SM_width_correction)')
          self.modelBuilder.factory_('expr::mm_scaling("2*@0*@1*@2", a1, a3, SM_width_correction)')

        self.modelBuilder.factory_('expr::muV_mutautau("@0*@1", muV, mutautau)')
        for x in ['sm_scaling', 'ps_scaling', 'mm_scaling']:
            self.modelBuilder.factory_('expr::muV_%s("@0*@1", muV, %s)' % (x,x))
            self.modelBuilder.factory_('expr::muggH_%s("@0*@1", muggH, %s)' % (x,x))

    def getYieldScale(self, bin_, process):

        scalings = []

        if 'ggf' in process: 
            scalings.append('muggH')
        if 'qqH' in process or 'ZH' in process or 'WH' in process: 
            scalings.append('muV')

        if ('ggf' in process or 'qqH' in process or 'WH' in process or 'ZH' in process) and 'hww' not in process:
            if "sm" in process:
                scalings.append('sm_scaling')
            elif "cpo" in process:
                scalings.append('ps_scaling')
            elif "mm" in process:
                scalings.append('mm_scaling')
        if scalings:
            scaling = '_'.join(scalings)
            
            print('Scaling %s/%s as %s' % (bin_, process,scaling))
            return scaling
        else:
            return 1

CPMixtureDecays = CPMixtureDecays()



