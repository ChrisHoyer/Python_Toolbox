#############################################################################
#   Importiert s2p mithilfe von scikit-rf
#      
#   - Network2Sparam: s2p -> S Paramter Dictionary
#   - SParam2HParam: S Paramter Dictionary -> H Paramter Dictionary
#   - SParam2GParam: S Paramter Dictionary -> G Paramter Dictionary
#
#   Autor: C. Hoyer 
#   Stand: 14.02.19
#############################################################################

import skrf as rf
import numpy as np

#############################################################################
###         Extract S Param from Networkfile (.s2p)
#############################################################################
def Network2Sparam(File, S3P=False, S1P=False, **kwargs):
    
    # generate Network
    Network = rf.Network(File)
    
     # Extract S Param
    S11 = Network.s11.s.reshape(-2)
    
    if not S1P:
        S12 = Network.s12.s.reshape(-2)
        S21 = Network.s21.s.reshape(-2)
        S22 = Network.s22.s.reshape(-2)

    freq = Network.frequency.f
    
    if S3P:
        S33 = Network.s33.s.reshape(-2)        
         
    #return types
    SNetwork = {}
    SNetwork['S11'] = S11
    
    if not S1P:
        SNetwork['S12'] = S12
        SNetwork['S21'] = S21
        SNetwork['S22'] = S22
             
    #return types in dB
    SNetwork['S11_db'] = 20*np.log10(S11)
    
    if not S1P:
        SNetwork['S12_dB'] = 20*np.log10(S12)
        SNetwork['S21_dB'] = 20*np.log10(S21)
        SNetwork['S22_dB'] = 20*np.log10(S22)
    
    if S3P:
        SNetwork['S33'] = S33
        SNetwork['S33_dB'] = 20*np.log10(S33)
    
    #return freq
    SNetwork['freq'] = freq
    
    
    return SNetwork
    
#############################################################################
###         Extract H Param from S Param
#############################################################################
def Network2HParam(File, convert=False, **kwargs):
    
        # generate Network
    Network = rf.Network(File)
    
    # Extract S Param
    S11 = Network.s11.s.reshape(-2)
    S12 = Network.s12.s.reshape(-2)
    S21 = Network.s21.s.reshape(-2)
    S22 = Network.s22.s.reshape(-2)
    freq = Network.frequency.f
    
    # convert to dB
    if convert:
        S11 = 20*np.log10(S11)
        S12 = 20*np.log10(S12)
        S21 = 20*np.log10(S21)
        S22 = 20*np.log10(S22)
    
    # Equations from HHHS Building Blocks - Active Devices Page 29
    deltaS = S11*S22-S12*S21
    denominator = -deltaS - S11 + S22 + 1
    
    # H Param
    H11 = (deltaS + S11 + S22 + 1)/denominator
    H12 = (2*S12)/denominator
    H21 = (-2*S21)/denominator
    H22 = (deltaS - S11 + S22 + 1)/denominator
    
    #return types
    HNetwork = {}
    HNetwork['H11'] = H11
    HNetwork['H12'] = H12
    HNetwork['H21'] = H21
    HNetwork['H22'] = H22
    
    #return types in dB
    HNetwork['H11_db'] = 20*np.log10(H11)
    HNetwork['H12_db'] = 20*np.log10(H12)
    HNetwork['H21_db'] = 20*np.log10(H21)
    HNetwork['H22_db'] = 20*np.log10(H22)
    
    #return freq
    HNetwork['freq'] = freq
    
    
    return HNetwork

#############################################################################
###         Extract G Param from S Param
#############################################################################  
def Network2GParam (File, convert=False, **kwargs):
    
        # generate Network
    Network = rf.Network(File)
    
    # Extract S Param
    S11 = Network.s11.s.reshape(-2)
    S12 = Network.s12.s.reshape(-2)
    S21 = Network.s21.s.reshape(-2)
    S22 = Network.s22.s.reshape(-2)
    freq = Network.frequency.f
    
    # convert to dB
    if convert:
        S11 = 20*np.log10(S11)
        S12 = 20*np.log10(S12)
        S21 = 20*np.log10(S21)
        S22 = 20*np.log10(S22)
    
    # Equations from HHHS Building Blocks - Active Devices Page 29
    deltaS = S11*S22-S12*S21
    denominator = -deltaS - S11 + S22 + 1
    
    # H Param
    H11 = (deltaS + S11 + S22 + 1)/denominator
    H12 = (2*S12)/denominator
    H21 = (-2*S21)/denominator
    H22 = (deltaS - S11 + S22 + 1)/denominator
    
    # From H to G Param
    detH = []
    for i in range(len(H11)):
        HMatrix = np.array([[H11[i], H12[i]], [H21[i], H22[i]]])
        detH.append(np.linalg.det(HMatrix))
        
    # set dataformat
    detH = np.array(detH)
    
    # G Matrix, Equations from Wikipedia "Two Port Network"
    G11 = 1/detH * H22
    G12 = 1/detH * -H12
    G21 = 1/detH * -H21
    G22 = 1/detH * H11
    
    #return types
    GNetwork = {}
    GNetwork['G11'] = G11
    GNetwork['G12'] = G12
    GNetwork['G21'] = G21
    GNetwork['G22'] = G22
    
    #return types in dB
    GNetwork['G11_db'] = 20*np.log10(G11)
    GNetwork['G12_db'] = 20*np.log10(G12)
    GNetwork['G21_db'] = 20*np.log10(G21)
    GNetwork['G22_db'] = 20*np.log10(G22)
    
    #return freq
    GNetwork['freq'] = freq
    
    
    return GNetwork
    
#############################################################################
###         Extract Z Param from Networkfile (.s2p)
#############################################################################
def Network2Zparam(File, Z1P=False, **kwargs):
    
        # generate Network
    Network = rf.Network(File)
    
     # Extract S Param
    Z11 = Network.z[:,0,0]

    Z0 = Network.z0
    freq = Network.frequency.f
        
    if not Z1P:
        Z12 = Network.z[:,0,1]
        Z21 = Network.z[:,1,0]
        Z22 = Network.z[:,1,1]

    freq = Network.frequency.f
        
    #return types
    ZNetwork = {}
    ZNetwork['Z11'] = Z11
    
    if not Z1P:
        ZNetwork['Z12'] = Z12
        ZNetwork['Z21'] = Z21
        ZNetwork['Z22'] = Z22
                  
    # umgebung
    ZNetwork['Z0'] = Z0   
    #return freq
    ZNetwork['freq'] = freq
    
    
    return ZNetwork
 


