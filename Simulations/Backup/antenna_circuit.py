"""
Arnold's Laws of Documentation:
    (1) If it should exist, it doesn't.
    (2) If it does exist, it's out of date.
    (3) Only documentation for useless programs transcends the first two laws.

Created on 27 Jun 2022

@author: frederic
"""
__updated__ = "2022-08-26 17:13:09"

import os

# set the number of threads for the OpenBLAS libs
#   (needs to happen before importing numpy)

os.environ['OPENBLAS_NUM_THREADS']='{:d}'.format(2)

import numpy as np
import matplotlib.pyplot as pl
from scipy.constants import speed_of_light as c0


from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject

pF_, MHz_ = 1e-12, 1e6 # unit conversion factors 
OhmSign = u'\N{GREEK CAPITAL LETTER OMEGA}'

os.makedirs('figures', exist_ok=True)

#===============================================================================
#
# T O M A S
#
class TOMAS():
    #===========================================================================
    #
    # _ _ i n i t _ _ 
    #
    def __init__(self, **kwargs):
        """
        build the tomas antenna circuit:
        
          sc [A1] t   t  [A2] cap ant  vcw vcw  cap cap [Ca] sc
        |>---------+ +-----------+--------+--------+-----||----<|
                    +             [A2Ca1]   [A2Ca]
                    | t
                    |
                    | [T2Cs1]
                    |
                    | vcw
                    +
                    | vcw
                    |
                    | [T2Cs]
                    |
                    | cap
                    +
                    | cap
                    |
                   === [Cs]
                    |
                    | toCp
                    +
                    | Cs
                    |
                    | [Cs2Cp]
                    |
                    | Cp   Cp [CptoV0] V0 V0 [V0toV1] V1 V1 [V1toV2] V2 V2 [V2toV3] V3
                    +     +--------------+--------------+--------------+--------------+
                    | Cp
                    |
                   === [Cp]
                    |
                    | sc
                  |>+
                
        """
        
        self.Zbase = kwargs.pop('Zbase', 50.)
        self.CT = rfCircuit(Zbase = self.Zbase)
        
        # antenna strap
        
        self.tsStrap = kwargs.pop('tsStrap', None)
        
        self.A2Ca1_L = kwargs.pop('A2Ca1_L', 0.1715)
        self.A2Ca1_Z = kwargs.pop('A2Ca1_Z', 50.)
        self.T2Cs1_L = kwargs.pop('T2Cs1_L', 0.121)
        self.T2Cs1_Z = kwargs.pop('T2Cs1_Z', 50.)
        
        
        
        if self.tsStrap:
            print('using touchstone file for straps')
            self.strap = rfObject(touchstone=self.tsStrap, ports=['cap','t'])
            
            if self.strap.fs[0] == 0.:
                self.strap.fs[0] = 1e-3 # avoid divisions by 0 at 0 Hz in deembed
            # print(self.strap.fs)
            self.strap.deembed({'cap':(self.A2Ca1_L, self.A2Ca1_Z), 
                                't'  :(self.T2Cs1_L, self.T2Cs1_Z)})
                    
            ignored = [kw for kw in ['La1', 'La2', 'Z0a', 'v0a', 'Ra'] if kw in kwargs]
            for kw in ignored:
                kwargs.pop(kw,'ignored')
            if ignored:
                print(f'kwargs {", ".join(ignored)} ignored because touchstone '
                      'file (tsStrap) was provided.')

        else:
            print(f'using TL-model of straps')
            self.La1 = kwargs.pop('La1', 0.19 + 0.045) # feeders to strap lengths
            self.La2 = kwargs.pop('La2', 0.19 + 0.045) # feeders to strap lengths
            self.Z0a = kwargs.pop('Z0a', 45.0)
            self.v0a = kwargs.pop('v0a', 1.0)
            self.Ra = kwargs.pop('Ra', 2.0)
            
            self.Lant = self.Z0a / (self.v0a * c0)
            self.Cant = 1 / (self.Z0a * self.v0a * c0)
        
            self.strap =rfCircuit(Zbase=self.Zbase)
            self.strap.addblock(
                'A1', 
                rfTRL(L=self.La1, rTL=self.Ra, LTL=self.Lant, CTL=self.Cant, 
                      ports= ['sc','t'])
            )
            self.strap.addblock(
                'A2', 
                rfTRL(L=self.La2, rTL=self.Ra, LTL=self.Lant, CTL=self.Cant, 
                      ports= ['t','cap'])
            )
            self.strap.terminate('A1.sc', Z=0.)
                        
            self.strap.connect('A1.t', 'A2.t', 't')
            

            self.strap.connect('A2.cap', 'cap')
            
        # print(self.strap.__str__(full=1))
        
        # circuit
        
        self.T2Cs_L = kwargs.pop('T2Cs_L', 0.16)
        self.T2Cs_Z = kwargs.pop('T2Cs_Z', 50.)

        self.Cs2Cp_L = kwargs.pop('Cs2Cp_L', 0.2)
        self.Cs2Cp_Z = kwargs.pop('Cs2Cp_Z', 50.)
        
        self.A2Ca_L = kwargs.pop('A2Ca_L', 0.29)
        self.A2Ca_Z = kwargs.pop('A2Ca_Z', 50.)
        
        self.CptoV0_L = kwargs.pop('CptoV0_L', 0.235)
        self.CptoV0_Z = kwargs.pop('CptoV0_L', 50.)

        self.V0toV1_L = kwargs.pop('V0toV1_L', 0.66)
        self.V0toV1_Z = kwargs.pop('V0toV1_L', 50.)

        self.V1toV2_L = kwargs.pop('V1toV2_L', 0.795)
        self.V1toV2_Z = kwargs.pop('V1toV2_L', 50.)
         
        self.V2toV3_L = kwargs.pop('V2toV3_L', 0.66)
        self.V2toV3_Z = kwargs.pop('V2toV3_L', 50.)
        
        self.LsCaps_H = kwargs.pop('LsCaps_H', 20e-9)
        
        if kwargs:
            print(f'unused kwargs: {", ".join(kwargs)}')

        # build the circuit
        
        self.CT.addblock('strap', self.strap, 
                         # ports=['t','cap']
                         )
        # print(self.CT.__str__(full=1))
        
        self.CT.addblock('A2Ca1',
            rfTRL(L=self.A2Ca1_L, Z0TL=self.A2Ca1_Z,
            ports= ['ant','vcw'])
        )
        
        self.CT.connect('strap.cap','A2Ca1.ant')
        
        self.CT.addblock('A2Ca', 
            rfTRL(L=self.A2Ca_L, Z0TL=self.A2Ca_Z,
            ports= ['vcw','cap'])
        )
        self.CT.connect('A2Ca1.vcw', 'A2Ca.vcw')
        
        self.CT.addblock('Ca', rfRLC(Cs=1e-12, Ls=self.LsCaps_H, ports= ['cap','sc']))
        
        self.CT.connect('Ca.cap', 'A2Ca.cap')
        self.CT.terminate('Ca.sc',Z=0.)
        
        self.CT.addblock('T2Cs1',
            rfTRL(L=self.T2Cs1_L, Z0TL=self.T2Cs1_Z,
            ports= ['t','vcw'])
        )
        
        self.CT.connect('strap.t', 'T2Cs1.t')
        
        self.CT.addblock('T2Cs',
            rfTRL(L=self.T2Cs_L, Z0TL=self.T2Cs_Z,
            ports= ['vcw','cap'])
        )
        
        self.CT.connect('T2Cs1.vcw', 'T2Cs.vcw')
        
        self.CT.addblock('Cs', rfRLC(Cs=1e-12, Ls=self.LsCaps_H, ports= ['cap','toCp']))
        self.CT.connect('T2Cs.cap','Cs.cap')
        
        self.CT.addblock('Cs2Cp',
            rfTRL(L=self.Cs2Cp_L, Z0TL=self.Cs2Cp_Z,
            ports= ['Cs','Cp'])
        )

        self.CT.connect('Cs.toCp', 'Cs2Cp.Cs')
        
        self.CT.addblock('Cp', rfRLC(Cs=1e-12, Ls=self.LsCaps_H, ports= ['Cp','sc']))
        self.CT.terminate('Cp.sc', Z=0)
        
        self.CT.addblock('CptoV0', 
            rfTRL(L=self.CptoV0_L, Z0TL=self.CptoV0_Z,
            ports= ['Cp','V0'])
        )
    
        self.CT.connect('Cs2Cp.Cp', 'Cp.Cp', 'CptoV0.Cp')

        self.CT.addblock('V0toV1', 
            rfTRL(L=self.V0toV1_L, Z0TL=self.V0toV1_Z,
            ports= ['V0','V1'])
        )
        self.CT.addblock('V1toV2', 
            rfTRL(L=self.V1toV2_L, Z0TL=self.V1toV2_Z,
            ports= ['V1','V2'])
        )
    
        self.CT.connect('CptoV0.V0','V0toV1.V0')
        self.CT.connect('V0toV1.V1','V1toV2.V1')

        self.CT.addblock('V2toV3', 
            rfTRL(L=self.V2toV3_L, Z0TL=self.V2toV3_Z,
            ports= ['V2','V3'])
        )

        self.CT.connect('V1toV2.V2','V2toV3.V2')
        
    #===========================================================================
    #
    # G e t C a p s 
    #
    def GetCaps(self, pF=1):
        return dict(
            [(ck, self.CT.blocks[ck]['object'].Cs/pF) for ck in ['Ca', 'Cs', 'Cp']])
    
    #===========================================================================
    #
    # S e t C a p s 
    #
    def SetCaps(self, **kwargs):
        for tC in ['Ca', 'Cs', 'Cp']:
            tCval = kwargs.pop(tC, None)
            if tCval is not None:
                self.CT.addblock(tC, rfRLC(Cs=tCval))
                
        if kwargs:
            raise ValueError(f'SetCaps: Unrecognized kwargs :{",".join(kwargs)}')
        
    #===========================================================================
    #
    # g e t V I s
    #
    def getVIs(self, fs, Agen=1.0):
        
        tmap = {
            'Vt':   ('strap.t', 0),
            'It':   ('strap.t', 1),
            'V1':   ('V0toV1.V1',0),
            'V2':   ('V1toV2.V2',0),
            'V3':   ('V2toV3.V3',0),
            'Ia':   ('A2Ca.cap',1),
            'Is':   ('Cs.cap',1),
            'Ip':   ('Cp.Cp',1),
            'VpS':  ('Cs2Cp.Cp', 0),
            'IpS':  ('Cs2Cp.Cp', 1),
            'Iant': ('strap.cap',1), 
            'Vfwd': ('V2toV3.V3',2),
            'Vrfl': ('V2toV3.V3',3),
        }
        
        if self.tsStrap == None:
            # we can only have those for a TL model of the straps
            tmap['IAL'] = ('strap.A1.t',1)
            tmap['IAU'] = ('strap.A2.t',1)
            
        tSol = {}
        for kw in tmap:
            tSol[kw] = []

        for fs0 in fs if hasattr(fs,'__iter__') else [fs]:
            Sol = self.CT.Solution(fs0, {'V2toV3.V3':Agen})
            
            # for kw, val in Sol.items():
            #     print(f'{kw:<15s}: {val}')
            
            for kw, (node, typ) in tmap.items():
                tSol[kw].append(Sol[node][typ])
                
        return tSol
    
    #===========================================================================
    #
    # e r r o r s _ C s C p 
    #
    def errors_CsCp(self, f, dgdCs=1, dbdCp=-1, dCmax=10E-12):
        gm = self.CT.getS(f)[0]
        tSol = self.CT.Solution(f, {'V2toV3.V3':np.sqrt(2*50*5*1E3)}, nodes=['CptoV0.Cp'])
        yCp = -tSol['CptoV0.Cp'][1]/tSol['CptoV0.Cp'][0]*self.Zbase
        dCs = max(-dCmax, min((yCp.real-1)*dgdCs, dCmax))
        dCp = max(-dCmax, min(yCp.imag*dbdCp, dCmax))
        return gm, dCs, dCp
        
    #===========================================================================
    #
    # m a t c h 
    #
    def match(self, f, gmmax=0.01, maxstep=1000, dgdCs=0.002, dbdCp=-0.002, 
              dCmax=10E-12, Cmin=8e-12, Cmax=1000e-12):
        gms, Css, Cps = [], [], []
        for kstep in range(maxstep):
            Cks = self.GetCaps()
            Cs, Cp = Cks['Cs'], Cks['Cp']
            gm, dCs, dCp = self.errors_CsCp(f, dgdCs, dbdCp, dCmax)
            gms.append(gm)
            Css.append(Cs)
            Cps.append(Cp)
            nCs = max(Cmin, min(Cs+dCs, Cmax))
            nCp = max(Cmin, min(Cp+dCp, Cmax))
            self.SetCaps(Cs=nCs, Cp=nCp)
            if np.abs(gm) < gmmax:
                break
        
        return gms, Css, Cps

    #===========================================================================
    #
    # Realistic Match
    #
    def RealisticMatch(self, f, gmmax=0.01, maxstep=1000, dgdCs=0.002, dbdCp=-0.002, 
              dCmax=10E-12, Cmin=8e-12, Cmax=1000e-12):
        gms, Css, Cps = [], [], []
        for kstep in range(maxstep):
            Cks = self.GetCaps()
            Cs, Cp = Cks['Cs'], Cks['Cp']
            gm, dCs, dCp = self.errors_CsCp(f, dgdCs, dbdCp, dCmax)
            gms.append(gm)
            Css.append(Cs)
            Cps.append(Cp)
            nCs = max(Cmin, min(Cs+dCs, Cmax))
            nCp = max(Cmin, min(Cp+dCp, Cmax))
            self.SetCaps(Cs=nCs, Cp=nCp)
            if np.abs(gm) < gmmax:
                break
        
        return gms, Css, Cps

    
#===============================================================================
#
# g e t S 
#
def getS(par, ant, f0MHz, Cmin=7., Cmax=1000., verbose=False):
       
    ant.SetCaps(**dict([(kC, val*pF_) for kC, val in zip(['Cs', 'Cp', 'Ca'], par)]))
    f0 = np.squeeze(np.abs(ant.CT.getS(f0MHz * MHz_)))
    if len(par) == 1:
        w = 2e6*np.pi*f0MHz
        # in this case Ca is fixed and Cg is determined by trying to cancel the
        # remaining reactamce at the Cg node
        tSol = ant.getVIs(f0MHz*MHz_)
        YP = tSol['IpS'][0]/tSol['VpS'][0]
        yCp = np.imag(YP)
        #Cp = min(max(-yCp/(1 - w*ant.LsCaps_H*yCp)/w, Cmin*pF_), Cmax*pF_)
        Cp = min(max(-yCp/w, Cmin*pF_), Cmax*pF_)
        ant.SetCaps(Cp=Cp)
        f0 = np.squeeze(np.abs(ant.CT.getS(f0MHz * MHz_)))
        tSols = ant.CT.Solution(f0MHz*MHz_, {'V2toV3.V3':np.sqrt(2*50*5.0*1E3)}, nodes=['CptoV0.Cp'])
        YCp = -tSols['CptoV0.Cp'][1]/tSols['CptoV0.Cp'][0]
        if verbose:
            print(ant.CT.__str__(full=2))
            print(YP*ant.Zbase, Cp/pF_, YCp*ant.Zbase)
    # print(f'Ca={Ca:.3f}pF, Cs={Cs:.3f}pF, Cp={Cp:.3f}pF, GdB={20*np.log10(f0):.3f}')
    return f0

#===============================================================================
#
# s c a n _ C a
#
def scan_Ca(ant, Cas):
    
    if ant.tsStrap:
        model = f'{os.path.basename(ant.tsStrap)}'
    else:
        model = f'TL-model Ra={ant.Ra:0.2f} Ohm'
    model += f' @ f={fMHz0:0.1f}MHz'

    res = dict(
        [(kw, []) for kw in [
            'Ca', 'Ia', 'Cs', 'Is', 'Cp', 'Ip', 'Iant', 'It', 'Vt', 'Pt', 'gamma',
            'IAL', 'IAU'
            ]
        ])
    
    print(f'case:{CASE} @ {fMHz0:.1f} MHz\n')
    #      12345678 12345678 12345678 12345678 12345678 12345678 12345678 12345678 12345678 12345678 12345678
    print('   Ca/pF     Ia/A    Cs/pF     Is/A    Cp/pF     Ip/A   Iant/A     It/A     Vt/V    Pt/kW  |gamma|')
    
    for Ca in Cas: 

        ant.SetCaps(Ca=Ca*pF_)
        
        # print(TomasAnt.GetCaps(pF))
        nfev = 0
        if True:
            # perform Cs scan
            Csmin, gmmin = -1, +np.Inf
            
            for Cs in np.linspace(8,1000,992//32 +1):
                gm = getS([Cs], ant, fMHz0)
                nfev += 1
                if gm < gmmin:
                    gmmin, Csmin = gm, Cs

        else:
            Csmin = 500
            
        ant.SetCaps(Cs=Csmin*pF_)
            
        tmin = minimize(getS, [ant.GetCaps(pF_)['Cs']], 
                        bounds=[(8,1000)],
                        args=(ant, fMHz0), tol=1e-9)
        # print(tmin)
        Csmin = tmin.x[0]
        nfev += tmin.nfev
            
        fmin = getS([Csmin], ant, fMHz0, verbose=False) # this will set Cs and adjust Cp to the optimal
        
        # print(TomasAnt.GetCaps(pF))
        tSol = ant.getVIs(fMHz0*MHz_, Agen=Agen)
        # pprint(tSol)
    
        for ck in ['a','s','p']:
            res[f'C{ck}'].append(TomasAnt.CT.blocks[f"C{ck}"]["object"].Cs/pF_)
            res[f'I{ck}'].append(np.abs(tSol[f"I{ck}"][0]))
        for ck in ['Iant', 'It','Vt']:
            res[ck].append(np.abs(tSol[ck][0]))
        res['Pt'].append(0.5E-3 * np.real(np.conj(tSol["It"][0])*tSol["Vt"][0]))
        # res['gamma'].append(tmin.fun)
        res['gamma'].append(np.abs(tSol[f"Vrfl"][0]/tSol[f"Vfwd"][0]))
        try:
            res['IAL'].append(tSol['IAL'][0])
            res['IAU'].append(tSol['IAU'][0])
        except:
            res['IAL'].append(None)
            res['IAU'].append(None)

        
        print(
            f'{ant.CT.blocks["Ca"]["object"].Cs/pF_:8.2f} {np.abs(tSol["Ia"][0]):8.2f} '
            f'{Csmin:8.2f} {np.abs(tSol["Is"][0]):8.2f} '
            f'{ant.CT.blocks["Cp"]["object"].Cs/pF_:8.2f} {np.abs(tSol["Ip"][0]):8.2f} '
            f'{np.abs(tSol["Iant"][0]):8.1f} '
            f'{np.abs(tSol["It"][0]):8.1f} '
            f'{np.abs(tSol["Vt"][0]):8.1f} '
            f'{0.5E-3 * np.real(np.conj(tSol["It"][0])*tSol["Vt"][0]):8.3f} '
            f'{fmin:8.3f} ({nfev})')

    kmin, best = -1, +np.Inf
    for k, gamma in enumerate(res['gamma']):
        if gamma <= 0.001 and res['It'][k] < best:
            kmin, best = k, res['It'][k]
            Ca, Cs, Cp = res['Ca'][k], res['Cs'][k], res['Cp'][k]
            
    if kmin == -1:
        print('no matched solutions found for the scanned Ca\'s!')
        
    else:
        print(f'kmin= {kmin}, Ca= {Ca:7.1f}, Cs= {Cs:7.1f}, Cp= {Cp:7.1f}')

        ant.SetCaps(Cs=Cs*pF_, Ca=Ca*pF_, Cp=Cp*pF_)
        rhos = ant.CT.getS(fMHzs*MHz_)
                
        tfig, axs = pl.subplots(3, 2, squeeze=True, sharex=False, 
                                figsize=(12,12))
        
        pl.suptitle(f'{model}', fontsize=16)
        
        pl.sca(axs[0,0])
        pl.plot(res['Ca'], res['Ia'], 'r-', label='Ia')
        pl.plot(res['Ca'], res['Is'], 'b-', label='Is')
        pl.plot(res['Ca'], res['Ip'], 'g-', label='Ip')
        pl.plot(res['Ca'], res['Iant'], 'm-', label='Iant')
        pl.plot(res['Ca'], res['It'], 'c-', label='It')
        pl.axvline(Ca, color='k', ls='--', lw=2)
        
        pl.title('Currents')
        pl.xlabel('C$_a$ [pF]')
        pl.ylabel('Ia, Is, Ig [A]')
        pl.legend(loc='best')
        pl.grid(True)
        pl.axvline(Ca, color='k', ls='--', lw=2)
        
        ## T voltage and current ratio
        
        pl.sca(axs[1,0])
        pl.plot(res['Ca'], res['Vt'],'r-', label='Vt')
        pl.title('T-voltage')
        pl.xlabel('C$_a$ [pF]')
        pl.ylabel('Vt [V]')
        pl.legend(loc='best')
        pl.grid(True)
        # ax2.ylabel('Current ratio Iant/It')
        pl.axvline(Ca, color='k', ls='--', lw=2)
        
        ##
        ## Gamma
        ##
        pl.sca(axs[2,0])
        pl.plot(res['Ca'], res['gamma'],'r-', label='|$\Gamma$|')
        
        pl.title('|$\Gamma$|')
        pl.xlabel('C$_a$ [pF]')
        pl.ylabel('|$\Gamma$|')
        pl.legend(loc='best')
        pl.grid(True)
        pl.axvline(Ca, color='k', ls='--', lw=2)
        
        ##
        ## Capacitors
        ##
        pl.sca(axs[0,1])
        pl.plot(res['Ca'], res['Cp'],'g-', label='Cp')
        pl.plot(res['Ca'], res['Cs'],'b-', label='Cs')
        
        pl.title('Capacitor values')
        pl.xlabel('C$_a$ [pF]')
        pl.ylabel('Cp, Cs [pF]')
        pl.legend(loc='best')
        pl.grid(True)
        pl.axvline(Ca, color='k', ls='--', lw=2)
    
        ##
        ## mismatch vs frequency
        ##
        pl.sca(axs[2,1])
        pl.plot(fMHzs, np.squeeze(np.abs(rhos)), 
                label=f'Ca={Ca:.1f}pF, Cs={Cs:.1f}pF,Cp={Cp:.1f}pF'
        )
        pl.xlabel('frequency [MHz]')
        pl.ylabel('|$\Gamma_{V_3}$|')
        # pl.ylim(top=1,bottom=-60)
        pl.title('mismatch vs. frequency')
        pl.legend(loc='best')
        pl.tight_layout()
        pl.grid(True)
        pl.axvline(fMHz0, color='k', ls='--', lw=2)
        pl.ylim(bottom=0.)
        
        ##
        ## dGdCs
        ##
        pl.sca(axs[1,1])
        if False:
            dCs, dgs = 5e-1, []
            for Cak, Csk, Cpk in zip(res['Ca'], res['Cs'], res['Cp']):
                
                ant.SetCaps(Ca=Cak*pF_, Cs=(Csk+dCs)*pF_, Cp=Cpk*pF_)
                tSol = ant.CT.Solution(fMHz0*MHz_, {'V2toV3.V3':np.sqrt(2*50*5*1E3)})
                dg = np.real(-tSol['CptoV0.Cp'][1]/tSol['CptoV0.Cp'][0] * TomasAnt.Zbase)
                
                ant.SetCaps(Cs=(Csk-dCs)*pF_)
                tSol = TomasAnt.CT.Solution(fMHz0*MHz_, {'V2toV3.V3':np.sqrt(2*50*5*1E3)})
                dg -= np.real(-tSol['CptoV0.Cp'][1]/tSol['CptoV0.Cp'][0] * TomasAnt.Zbase)
                
                dgs.append(dg/2/dCs)
                
            pl.plot(res['Ca'], dgs, label=f'f={fMHz0:.1f}MHz')
            
            pl.title('dG/dCs')
            pl.xlabel('C$_a$ [pF]')
            pl.ylabel('dG/dCs[Mho/pF]')
            pl.legend(loc='best')
            pl.grid(True)
            pl.axvline(Ca, color='k', ls='--', lw=2)
        
        else:
            ##
            ## current ratios
            ##
            pl.plot(res['Ca'], (np.abs(res['It'])/np.abs(res['Iant'])),'b', 
                    label='|I$_{t}$ / I$_{ant}$|')
            
            if all([kw in tSol for kw in ['IAL','IAU']]):
                pl.plot(res['Ca'], (np.abs(res['IAL'])/np.abs(res['IAU'])), 'r',
                        label = '|I$_{t,L}$ / I$_{t,U}$|')
                
                pl.plot(res['Ca'], 
                        np.angle(np.array(res['IAU'])/np.array(res['IAL']),deg=1)/180, 
                        'g',
                        label = '$\phi$(I$_{t,U}$ / I$_{t,L}$)')
                        
                pl.plot()
            pl.xlabel('C$_a$ [pF]')
            pl.ylabel('|I$_t$/I$_{ant}$|, |I$_{t,L}$/I$_{t,U}$|, $\phi/\pi$')
            pl.title('Current ratios')
            pl.grid(True)
            pl.legend(loc='best')
            pl.axvline(Ca, color='k', ls='--', lw=2)
            pl.ylim(bottom=0., top=3.)

        pl.tight_layout()
                
        pl.savefig(f'figures/{model}.png')

    return res, kmin

#===============================================================================
#
# _ _ m a i n _ _
#
if __name__ == '__main__':
    
    from scipy.optimize import minimize
    import csv
    
    fMHzs = np.linspace(5,50,261)
    fMHz0s = np.linspace(25,50,25)
    PkW = 5.00

    f = open('MatchedSystems.csv','w')
    writer = csv.writer(f)
    CASE = 'Antenna/Tomas-Ref_geo-R=200-Diel_eps=0500.s2p'

    TomasAnt = TOMAS(
            tsStrap=CASE if CASE != 'TL model' else None,
            Ra = 3,                     # Ohm/m
            LsCaps_H = 30e-9            # H
        )

    for fMHz0 in fMHz0s:
        Vmatched = np.sqrt(2 * 50. * PkW *1E3)
        print(f'Vmatched = {Vmatched:0.1f} V')
            # CASE = 'TL model'
        
        if TomasAnt.tsStrap:
            model = f'{os.path.basename(TomasAnt.tsStrap)}'
        else:
            model = f'TL-model Ra={TomasAnt.Ra:0.2f} Ohm'
        model += f' @ f={fMHz0:0.1f}MHz'

        TomasAnt.SetCaps(Ca=20*pF_, Cs=750*pF_, Cp=750*pF_)
        # print(TomasAnt.CT.__str__(full=1))
        Agen = np.sqrt(2*TomasAnt.Zbase*PkW*1E3)
        
        ############################################################################
        ##                                                                        ##
        ##                     S C A N _ C A                                      ##
        ##                                                                        ##
        ############################################################################
        
        SCAN_CA = True
        
        if SCAN_CA:
            Cas = np.linspace(20,200,180+1)
            res, kmin = scan_Ca(TomasAnt, Cas)
            Ca = res['Ca'][kmin]
            
        ## END SCAN_CA == True #####################################################
                
        if not SCAN_CA or kmin == -1: # SCAN_CA == False
            
            Ca = 50. # pF
            
        TomasAnt.SetCaps(Ca=Ca*pF_)
            
        print(f'Computing dgdCs for Ca= {TomasAnt.GetCaps(pF_)["Ca"]:0.1f}pF '
              f'at {fMHz0:0.1f} MHz  ...')
        
        Css, yCps = np.linspace(7, 1000, 994), []
        for Cs in Css:
            TomasAnt.SetCaps(Cs=Cs*pF_)
            tSol = TomasAnt.getVIs(fMHz0*MHz_)
            yCp = tSol['IpS'][0]/tSol['VpS'][0]*TomasAnt.Zbase
            yCps.append(yCp)
        yCps = np.array(yCps)
        
        # find the best Cs that minimizes abs(real(yCp)-1))
        errgs = np.abs(yCps.real - 1)
        errgmin = np.min(errgs)
        Cs = Css[np.where(errgs == errgmin)][0]
        TomasAnt.SetCaps(Cs=Cs*pF_)
        fmin = getS([Cs], TomasAnt, fMHz0) # also sets Cp
        Cp = TomasAnt.GetCaps(pF_)['Cp']
        print(f'gamma= {fmin:0.6f} for Ca= {Ca:.1f}pF, Cs= {Cs:0.2f}pF and '
              f'Cp = {TomasAnt.GetCaps(pF_)["Cp"]:0.2f}pF at {fMHz0:0.1f}MHz')
        
        dgdCs = np.real(yCps[2:]-yCps[:-2])/(Css[2]-Css[0])
        
        # TomasAnt.SetCaps(Ca=Ca*pF, Cs=Cs*pF,Cp=Cp*pF)
        rhos = TomasAnt.CT.getS(fs=fMHzs*MHz_)
        
        print('simulating a match ...')
        
        TomasAnt.SetCaps(Cs=1000*pF_, Cp=700*pF_, Ca=Ca*pF_)
        gms, Css, Cps = TomasAnt.match(fMHz0*MHz_, gmmax=1e-2, 
                                       dgdCs=2E-11, dbdCp=-10E-11, dCmax=50e-12)
        writer.writerow([fMHz0,Cs,Cp,Ca])
        print([fMHz0,Css[-1],Cps[-1],Ca])
        print(np.array(gms).shape,np.array(Css).shape, np.array(Cps).shape)
        print(f'gamma= {np.squeeze(np.abs(gms))[-1]:0.6f} for Ca= {Ca:.1f}pF, '
              f'Cs= {Css[-1]/pF_:0.2f}pF and '
              f'Cp= {Cps[-1]/pF_:0.2f}pF at {fMHz0:0.1f}MHz')

f.close()
