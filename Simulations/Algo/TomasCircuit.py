from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject, rfCoupler

def ArthurMod(CaVal=150*1e-12,CpVal=47.94*1e-12,CsVal=133*1e-12):
    Zbase = 50
    ct = rfCircuit(Zbase = Zbase)

    # antenna strap

    A2Ca1_L = 0.95 #Antenna to pre-matching capacitor length
    A2Ca1_Z = 50
    T2Cs1_L = 0.087
    T2Cs1_Z = 50

    tsStrap = 'Ansys/TOMAS_plasma_box_HFSSDesign9.s2p'

    strap = rfObject(touchstone=tsStrap, ports=['cap','t'])

    if strap.fs[0] == 0.:
        strap.fs[0] = 1e-3 # avoid divisions by 0 at 0 Hz in deembed

    strap.deembed({'cap':(A2Ca1_L, A2Ca1_Z), 
                        't'  :(T2Cs1_L, T2Cs1_Z)})
            
    # circuit
    T2Cs_L = 0.16
    T2Cs_Z = 50
    Cs2Cp_L = 0.2
    Cs2Cp_Z = 50
    A2Ca_L = 0.29
    A2Ca_Z = 50.
    CptoV0_Z = 50
    V0toV1_Z = 50
    V1toV2_Z = 50
    V2toV3_Z = 50
    V3toIn_Z = 50

    CptoV0_L = 0.235
    V0toV1_L = 0.66
    V1toV2_L = 0.795
    V2toV3_L = 0.66
    V3toIn_L = 0.235

    LsCaps_H = 30e-9

    # build the circuit
    ct.addblock('strap', strap, 
                     # ports=['t','cap']
                     )
    ct.addblock('A2Ca1',
        rfTRL(L=A2Ca1_L, Z0TL=A2Ca1_Z,
        ports= ['ant','vcw'])
    )

    ct.connect('A2Ca1.ant','strap.cap')

    ct.addblock('A2Ca', 
        rfTRL(L=A2Ca_L, Z0TL=A2Ca_Z,
        ports= ['vcw','cap'])
    )
    ct.connect('A2Ca1.vcw', 'A2Ca.vcw')
    ct.addblock('Ca', rfRLC(Cs=CaVal, Ls=LsCaps_H, ports= ['cap','sc']))
    ct.connect('Ca.cap', 'A2Ca.cap')
    ct.terminate('Ca.sc',Z=0.)
    ct.addblock('T2Cs1',
        rfTRL(L=T2Cs1_L, Z0TL=T2Cs1_Z,
        ports= ['t','vcw'])
    )
    ct.connect('strap.t', 'T2Cs1.t')
    ct.addblock('T2Cs',
        rfTRL(L=T2Cs_L, Z0TL=T2Cs_Z,
        ports= ['vcw','cap'])
    )
    ct.connect('T2Cs1.vcw', 'T2Cs.vcw')
    ct.addblock('Cs', rfRLC(Cs=CsVal, Ls=LsCaps_H, ports= ['cap','toCp']))
    ct.connect('T2Cs.cap','Cs.cap')
    ct.addblock('Cs2Cp',
        rfTRL(L=Cs2Cp_L, Z0TL=Cs2Cp_Z,
        ports= ['Cs','Cp'])
    )
    ct.connect('Cs.toCp', 'Cs2Cp.Cs')
    ct.addblock('Cp', rfRLC(Cs=CpVal, Ls=LsCaps_H, ports= ['Cp','sc']))
    ct.terminate('Cp.sc', Z=0)


    # Coaxial cable

    ct.addblock('IntoV3', 
        rfTRL(L=V3toIn_L, Z0TL=V3toIn_Z,
        ports= ['In','V3']),
        relpos=0
    )
    ct.addblock('V3toV2', 
        rfTRL(L=V2toV3_L, Z0TL=V2toV3_Z,
        ports= ['V3','V2']),
        relpos=V3toIn_L
    )
    ct.addblock('V2toV1', 
        rfTRL(L=V1toV2_L, Z0TL=V0toV1_Z,
        ports= ['V2','V1']),
        relpos=V3toIn_L + V2toV3_L
    )
    ct.addblock('V1toV0', 
        rfTRL(L=V0toV1_L, Z0TL=V1toV2_Z,
        ports= ['V1','V0']),
        relpos=V3toIn_L + V2toV3_L + V0toV1_L
    )
    ct.addblock('V0toCp', 
        rfTRL(L=CptoV0_L, Z0TL=CptoV0_Z,
        ports= ['V0','Cp']),
        relpos=V3toIn_L + V2toV3_L + V0toV1_L + V1toV2_L
    )

    ct.connect('Cs2Cp.Cp', 'Cp.Cp', 'V0toCp.Cp')

    
    ct.connect('V0toCp.V0','V1toV0.V0')
    ct.connect('V1toV0.V1','V2toV1.V1')
    ct.connect('V2toV1.V2','V3toV2.V2')
    ct.connect('V3toV2.V3','IntoV3.V3')
    ct.connect('IntoV3.In','Source')

    return ct


def AntennaAndCaps(CaVal=150*1e-12,CpVal=47.94*1e-12,CsVal=133*1e-12):
    Zbase = 50
    ct = rfCircuit(Zbase = Zbase)

    # antenna strap

    A2Ca1_L = 0.95 #Antenna to pre-matching capacitor length
    A2Ca1_Z = 50
    T2Cs1_L = 0.087
    T2Cs1_Z = 50

    tsStrap = 'Ansys/TOMAS_plasma_box_HFSSDesign9.s2p'

    strap = rfObject(touchstone=tsStrap, ports=['cap','t'])

    if strap.fs[0] == 0.:
        strap.fs[0] = 1e-3 # avoid divisions by 0 at 0 Hz in deembed

    strap.deembed({'cap':(A2Ca1_L, A2Ca1_Z), 
                        't'  :(T2Cs1_L, T2Cs1_Z)})
            
    # circuit
    T2Cs_L = 0.16
    T2Cs_Z = 50
    Cs2Cp_L = 0.2
    Cs2Cp_Z = 50
    A2Ca_L = 0.29
    A2Ca_Z = 50.
    CptoV0_Z = 50
    V0toV1_Z = 50
    V1toV2_Z = 50
    V2toV3_Z = 50
    V3toIn_Z = 50

    CptoV0_L = 0.235
    V0toV1_L = 0.66
    V1toV2_L = 0.795
    V2toV3_L = 0.66
    V3toIn_L = 0.235

    LsCaps_H = 30e-9

    # build the circuit
    ct.addblock('strap', strap, 
                     # ports=['t','cap']
                     )
    ct.addblock('A2Ca1',
        rfTRL(L=A2Ca1_L, Z0TL=A2Ca1_Z,
        ports= ['ant','vcw'])
    )

    ct.connect('A2Ca1.ant','strap.cap')

    ct.addblock('A2Ca', 
        rfTRL(L=A2Ca_L, Z0TL=A2Ca_Z,
        ports= ['vcw','cap'])
    )
    ct.connect('A2Ca1.vcw', 'A2Ca.vcw')
    ct.addblock('Ca', rfRLC(Cs=CaVal, Ls=LsCaps_H, ports= ['cap','sc']))
    ct.connect('Ca.cap', 'A2Ca.cap')
    ct.terminate('Ca.sc',Z=0.)
    ct.addblock('T2Cs1',
        rfTRL(L=T2Cs1_L, Z0TL=T2Cs1_Z,
        ports= ['t','vcw'])
    )
    ct.connect('strap.t', 'T2Cs1.t')
    ct.addblock('T2Cs',
        rfTRL(L=T2Cs_L, Z0TL=T2Cs_Z,
        ports= ['vcw','cap'])
    )
    ct.connect('T2Cs1.vcw', 'T2Cs.vcw')
    ct.addblock('Cs', rfRLC(Cs=CsVal, Ls=LsCaps_H, ports= ['cap','toCp']))
    ct.connect('T2Cs.cap','Cs.cap')
    ct.addblock('Cs2Cp',
        rfTRL(L=Cs2Cp_L, Z0TL=Cs2Cp_Z,
        ports= ['Cs','Cp'])
    )
    ct.connect('Cs.toCp', 'Cs2Cp.Cs')
    ct.addblock('Cp', rfRLC(Cs=CpVal, Ls=LsCaps_H, ports= ['Cp','sc']))
    ct.terminate('Cp.sc', Z=0)


    ct.connect('Cs2Cp.Cp', 'Cp.Cp')

    return ct

def FredericCircuit(CaVal=150*1e-12,CpVal=47.94*1e-12,CsVal=133*1e-12,WhichStrap=0):
    Zbase = 50
    ct = rfCircuit(Zbase = Zbase)

    # antenna strap

    A2Ca1_L = 0.1715
    A2Ca1_Z = 50
    T2Cs1_L = 0.121
    T2Cs1_Z = 50

    Straps = ['Antenna/Tomas-Ref_geo-R=200-Diel_eps=0100.s2p','Antenna/Tomas-Ref_geo-R=200-Diel_eps=0200.s2p','Antenna/Tomas-Ref_geo-R=200-Diel_eps=0500.s2p','Antenna/Tomas-Ref_geo-R=200-Diel_eps=0842.s2p','Antenna/Tomas-Ref_geo-R=200-Diel_eps=1000.s2p','Antenna/Tomas-Ref_geo-R=200-Diel_eps=1500.s2p','Antenna/Tomas-Ref_geo-R=200-Diel_eps=2000.s2p','Antenna/Tomas-Ref_geo-R=200-Diel_eps=2500.s2p','Antenna/Tomas-Ref_geo-R=200-Diel_eps=3000.s2p']
    tsStrap = Straps[WhichStrap]


    strap = rfObject(touchstone=tsStrap, ports=['cap','t'])

    if strap.fs[0] == 0.:
        strap.fs[0] = 1e-3 # avoid divisions by 0 at 0 Hz in deembed

    strap.deembed({'cap':(A2Ca1_L, A2Ca1_Z), 
                        't'  :(T2Cs1_L, T2Cs1_Z)})
            
    # circuit
    T2Cs_L = 0.16
    T2Cs_Z = 50
    Cs2Cp_L = 0.2
    Cs2Cp_Z = 50
    A2Ca_L = 0.29
    A2Ca_Z = 50.
    CptoV0_Z = 50
    V0toV1_Z = 50
    V1toV2_Z = 50
    V2toV3_Z = 50
    V3toIn_Z = 50

    CptoV0_L = 0.235
    V0toV1_L = 0.66
    V1toV2_L = 0.795
    V2toV3_L = 0.66
    V3toIn_L = 0.235

    LsCaps_H = 30e-9

    # build the circuit
    ct.addblock('strap', strap, 
                     # ports=['t','cap']
                     )
    ct.addblock('A2Ca1',
        rfTRL(L=A2Ca1_L, Z0TL=A2Ca1_Z,
        ports= ['ant','vcw'])
    )
    ct.addblock('PlasmaLike',
        rfRLC(Rs=0.0,
        ports= ['in','out'])
    )
    ct.connect('PlasmaLike.in','A2Ca1.ant')
    ct.connect('strap.cap','PlasmaLike.out')
    ct.addblock('A2Ca', 
        rfTRL(L=A2Ca_L, Z0TL=A2Ca_Z,
        ports= ['vcw','cap'])
    )
    ct.connect('A2Ca1.vcw', 'A2Ca.vcw')
    ct.addblock('Ca', rfRLC(Cs=CaVal, Ls=LsCaps_H, ports= ['cap','sc']))
    ct.connect('Ca.cap', 'A2Ca.cap')
    ct.terminate('Ca.sc',Z=0.)
    ct.addblock('T2Cs1',
        rfTRL(L=T2Cs1_L, Z0TL=T2Cs1_Z,
        ports= ['t','vcw'])
    )
    ct.connect('strap.t', 'T2Cs1.t')
    ct.addblock('T2Cs',
        rfTRL(L=T2Cs_L, Z0TL=T2Cs_Z,
        ports= ['vcw','cap'])
    )
    ct.connect('T2Cs1.vcw', 'T2Cs.vcw')
    ct.addblock('Cs', rfRLC(Cs=CsVal, Ls=LsCaps_H, ports= ['cap','toCp']))
    ct.connect('T2Cs.cap','Cs.cap')
    ct.addblock('Cs2Cp',
        rfTRL(L=Cs2Cp_L, Z0TL=Cs2Cp_Z,
        ports= ['Cs','Cp'])
    )
    ct.connect('Cs.toCp', 'Cs2Cp.Cs')
    ct.addblock('Cp', rfRLC(Cs=CpVal, Ls=LsCaps_H, ports= ['Cp','sc']))
    ct.terminate('Cp.sc', Z=0)


    # Coaxial cable
    ct.addblock('CptoV0', 
        rfTRL(L=CptoV0_L, Z0TL=CptoV0_Z,
        ports= ['Cp','V0'])
    )

    ct.connect('Cs2Cp.Cp', 'Cp.Cp', 'CptoV0.Cp')

    ct.addblock('V0toV1', 
        rfTRL(L=V0toV1_L, Z0TL=V0toV1_Z,
        ports= ['V0','V1'])
    )
    ct.addblock('V1toV2', 
        rfTRL(L=V1toV2_L, Z0TL=V1toV2_Z,
        ports= ['V1','V2'])
    )
    ct.connect('CptoV0.V0','V0toV1.V0')
    ct.connect('V0toV1.V1','V1toV2.V1')

    ct.addblock('V2toV3', 
        rfTRL(L=V2toV3_L, Z0TL=V2toV3_Z,
        ports= ['V2','V3'])
    )
    ct.addblock('V3toIn', 
        rfTRL(L=V3toIn_L, Z0TL=V3toIn_Z,
        ports= ['V3','In'])
    )
    ct.connect('V1toV2.V2','V2toV3.V2')
    ct.connect('V2toV3.V3','V3toIn.V3')
    ct.connect('V3toIn.In','Source')

    return ct


def StrapCircuit2(CaVal=150*1e-12,CpVal=47.94*1e-12,CsVal=133*1e-12,WhichStrap=3):

    Zbase = 50.
    ct = rfCircuit(Zbase = Zbase)
    
    # antenna strap
    
    A2Ca1_L = 0.1715
    A2Ca1_Z = 50.
    T2Cs1_L = 0.121
    T2Cs1_Z = 50.

    Straps = ['Antenna/Tomas-Ref_geo-R=200-Diel_eps=0100.s2p','Antenna/Tomas-Ref_geo-R=200-Diel_eps=0200.s2p','Antenna/Tomas-Ref_geo-R=200-Diel_eps=0500.s2p','Antenna/Tomas-Ref_geo-R=200-Diel_eps=0842.s2p','Antenna/Tomas-Ref_geo-R=200-Diel_eps=1000.s2p','Antenna/Tomas-Ref_geo-R=200-Diel_eps=1500.s2p','Antenna/Tomas-Ref_geo-R=200-Diel_eps=2000.s2p','Antenna/Tomas-Ref_geo-R=200-Diel_eps=2500.s2p','Antenna/Tomas-Ref_geo-R=200-Diel_eps=3000.s2p']
    tsStrap = Straps[WhichStrap]

    strap = rfObject(touchstone=tsStrap, ports=['cap','t'])
    if strap.fs[0] == 0.:
        strap.fs[0] = 1e-3 # avoid divisions by 0 at 0 Hz in deembed

    strap.deembed({'cap':(A2Ca1_L, A2Ca1_Z), 
                                't'  :(T2Cs1_L, T2Cs1_Z)})
    T2Cs_L = 0.16
    T2Cs_Z = 50.

    Cs2Cp_L = 0.2
    Cs2Cp_Z = 50.

    A2Ca_L = 0.29
    A2Ca_Z = 50.

    CptoV0_L = 0.235
    CptoV0_Z = 50.
    V0toV1_L = 0.66
    V0toV1_Z = 50.
    V1toV2_L = 0.795
    V1toV2_Z = 50.
    V2toV3_L = 0.66
    V2toV3_Z = 50.
    V3toin_L = 0.235
    V3toin_Z = 50.
   
    LsCaps_H = 30e-9 #maybe 30?

    # build the circuit
        
    ct.addblock('strap', strap, 
                         # ports=['t','cap']
                         )
        # print(self.ct.__str__(full=1))
        
    ct.addblock('A2Ca1',
            rfTRL(L=A2Ca1_L, Z0TL=A2Ca1_Z,
            ports= ['ant','vcw'])
        )
        
    ct.connect('strap.cap','A2Ca1.ant')
        
    ct.addblock('A2Ca', 
            rfTRL(L=A2Ca_L, Z0TL=A2Ca_Z,
            ports= ['vcw','cap'])
        )
    ct.connect('A2Ca1.vcw', 'A2Ca.vcw')
        
    ct.addblock('Ca', rfRLC(Cs=CaVal, Ls=LsCaps_H, ports= ['cap','sc']))
        
    ct.connect('Ca.cap', 'A2Ca.cap')
    ct.terminate('Ca.sc',Z=0.)
        
    ct.addblock('T2Cs1',
            rfTRL(L=T2Cs1_L, Z0TL=T2Cs1_Z,
            ports= ['t','vcw'])
        )
        
    ct.connect('strap.t', 'T2Cs1.t')
        
    ct.addblock('T2Cs',
            rfTRL(L=T2Cs_L, Z0TL=T2Cs_Z,
            ports= ['vcw','cap'])
        )
        
    ct.connect('T2Cs1.vcw', 'T2Cs.vcw')
        
    ct.addblock('Cs', rfRLC(Cs=CsVal, Ls=LsCaps_H, ports= ['cap','toCp']))
    ct.connect('T2Cs.cap','Cs.cap')
        
    ct.addblock('Cs2Cp',
            rfTRL(L=Cs2Cp_L, Z0TL=Cs2Cp_Z,
            ports= ['Cs','Cp'])
        )

    ct.connect('Cs.toCp', 'Cs2Cp.Cs')
        
    ct.addblock('Cp', rfRLC(Cs=CpVal, Ls=LsCaps_H, ports= ['Cp','sc']))
    ct.terminate('Cp.sc', Z=0)
        
    ct.addblock('CptoV0',rfTRL(L=CptoV0_L, Z0TL=CptoV0_Z,ports= ['Cp','V0']))

    ct.addblock('V0toV1', 
            rfTRL(L=V0toV1_L, Z0TL=V0toV1_Z,
            ports= ['V0','V1'])
        )

    ct.addblock('V1toV2', 
            rfTRL(L=V1toV2_L, Z0TL=V1toV2_Z,
            ports= ['V1','V2'])
        )

    ct.addblock('V2toV3', 
            rfTRL(L=V2toV3_L, Z0TL=V2toV3_Z,
            ports= ['V2','V3'])
        )
    ct.addblock('V3toin', 
            rfTRL(L=V3toin_L, Z0TL=V3toin_Z,
            ports= ['V3','in'])
        )

    ct.connect('Cs2Cp.Cp', 'Cp.Cp', 'CptoV0.Cp')
    ct.connect('CptoV0.V0','V0toV1.V0')

    ct.connect('V0toV1.V1','V1toV2.V1')

    
    ct.connect('V1toV2.V2','V2toV3.V2')
    ct.connect('V2toV3.V3','V3toin.V3')
    ct.connect('Source','V3toin.in')

    return ct

def ArthurCircuit(CaVal=150*1e-12,CpVal=47.94*1e-12,CsVal=133*1e-12):

    # A simplified circuit

    TRLStoV3 = rfTRL(L=0.235, OD=0.041, ID=0.017, dx=500)
    TRLV3toV2 = rfTRL(L=0.66, OD=0.041, ID=0.017, dx=500)
    TRLV2toV1 = rfTRL(L=0.795, OD=0.041, ID=0.017, dx=500)
    TRLV1toV0 = rfTRL(L=0.66, OD=0.041, ID=0.017, dx=500)
    TRLV0toM = rfTRL(L=0.235, OD=0.041, ID=0.017, dx=500)
    """
    #load S matrix of antenna
    Antenna = rfObject(touchstone='Ansys/TOMAS_plasma_box_HFSSDesign9_test.s2p')
    RLCMatchBox = rfRLC(Cs=CsVal,Cp=CpVal)
    RLCPreMatchBox = rfRLC(Cs=CaVal)
    ct.addblock('Antenna',Antenna,ports=['G','H'])
    """
    #load S matrix of antenna and deembed:
    A2Ca1_L = 0.95 #Antenna to pre-matching capacitor length
    A2Ca1_Z = 50
    T2Cs1_L = 0.087
    T2Cs1_Z = 50
    tsStrap = 'Ansys/TOMAS_plasma_box_HFSSDesign9_test.s2p'
    Antenna = rfObject(touchstone=tsStrap, ports=['cap','t'])
    if strap.fs[0] == 0.:
        strap.fs[0] = 1e-3 # avoid divisions by 0 at 0 Hz in deembed
    Antenna.deembed({'cap':(A2Ca1_L, A2Ca1_Z), 't'  :(T2Cs1_L, T2Cs1_Z)})
    ct.addblock('Antenna', Antenna)
    ct.addblock('A2Ca1',
        rfTRL(L=A2Ca1_L, Z0TL=A2Ca1_Z,
        ports= ['ant','vcw'])
    )

    ct = rfCircuit(Id="TOMAS_ICRH")

    ct.addblock('StoV3',TRLStoV3,ports=['S','V3'],relpos=0)
    ct.addblock('V3toV2',TRLV3toV2,ports=['V3','V2'],relpos=TRLStoV3.L)
    ct.addblock('V2toV1',TRLV2toV1,ports=['V2','V1'],relpos=TRLStoV3.L + TRLV3toV2.L)
    ct.addblock('V1toV0',TRLV1toV0,ports=['V1','V0'],relpos=TRLStoV3.L + TRLV3toV2.L + TRLV2toV1.L)
    ct.addblock('V0toM',TRLV0toM,ports=['V0','M'],relpos=TRLStoV3.L + TRLV3toV2.L + TRLV2toV1.L + TRLV1toV0.L)

    EndOfLine =TRLStoV3.L + TRLV3toV2.L + TRLV2toV1.L + TRLV1toV0.L + TRLV0toM.L

    ct.addblock('Matching',RLCMatchBox,ports=['F','G'])
    ct.addblock('PreMatch',RLCPreMatchBox,ports=['H','I'])

    ct.connect('StoV3.S','Source')
    ct.connect('V3toV2.V3','StoV3.V3')
    ct.connect('V3toV2.V2','V2toV1.V2')
    ct.connect('V2toV1.V1','V1toV0.V1')
    ct.connect('V1toV0.V0','V0toM.V0')
    ct.connect('V0toM.M','Matching.F')
    ct.connect('Matching.G','Antenna.G')
    ct.connect('Antenna.H','PreMatch.H')
    ct.terminate('PreMatch.I',Z=0) #Grounding

    return ct
