# Preliminary Design of a Centrifugal Compressor Stage


import math
from math import atan,degrees,pi,radians,sin,tan,sqrt,cos,sinh
import matplotlib.pyplot as plt

 
#================================LIST===================================#


stage_list   = []   # Stage list index   [-]


p_01_list     = []   # Inlet impeller total pressure list      [Pa]
T_01_list     = []   # Inlet impeller total temperature list   [K]
h_01_list     = []   # Inlet impeller total enthalpy list      [J/kg]
rho_01_list   = []   # Inlet impeller total density list       [kg/m^3]
p_1_list      = []   # Inlet impeller pressure list            [Pa]
T_1_list      = []   # Inlet impeller temperature list         [K]
h_1_list      = []   # Inlet impeller enthalpy list            [J/kg]
rho_1_list    = []   # Inlet impeller density list             [kg/m^3]


p_02_list     = []   # Exit impeller total pressure list      [Pa]
T_02_list     = []   # Exit impeller total temperature list   [K]
h_02_list     = []   # Exit impeller total enthalpy list      [J/kg]
rho_02_list   = []   # Exit impeller total density list       [kg/m^3]
p_2_list      = []   # Exit impeller pressure list            [Pa]
T_2_list      = []   # Exit impeller temperature list         [K]
h_2_list      = []   # Exit impeller enthalpy list            [J/kg]
rho_2_list    = []   # Exit impeller density list             [kg/m^3]


T_04_list     = []   # Exit diffuser total temperature list   [K]
h_04_list     = []   # Exit diffuser total enthalpy list      [J/kg]
p_4_list      = []   # Exit diffuser pressure list            [Pa]
T_4_list      = []   # Exit diffuser temperature list         [K]
h_4_list      = []   # Exit diffuser enthalpy list            [J/kg]
rho_4_list    = []   # Exit diffuser density list             [kg/m^3]


p_04_list     = []   # Exit stage total pressure list      [Pa]
T_04_list     = []   # Exit stage total temperature list   [K]


c_1_list     = []   # Inlet impeller absolute velocity         [m/s]
w_1_h_list   = []   # Inlet impeller hub relative velocity     [m/s]
w_1_t_list   = []   # Inlet impeller tip relative velocity     [m/s]
u_1_h_list   = []   # Inlet impeller hub peripheral velocity   [m/s]
u_1_t_list   = []   # Inlet impeller tip peripheral velocity   [m/s]


c_2_list     = []   # Exit impeller absolute velocity     [m/s]
w_2_list     = []   # Exit impeller relative velocity     [m/s]
u_2_list     = []   # Exit impeller peripheral velocity   [m/s]


c_4_list     = []   # Exit diffuser absolute velocity   [m/s]  


alpha_1_list    = []   # Inlet impeller absolute flow angle (relative to radial)       [degree]
beta_1_h_list   = []   # Inlet impeller hub relative flow angle (relative to radial)   [degree]
beta_1_t_list   = []   # Inlet impeller tip relative flow angle (relative to radial)   [degree]


alpha_2_list   = []   # Exit impeller absolute flow angle (relative to radial)   [degree]
beta_2_list    = []   # exit impeller relative flow angle (relarive to radial)   [degree]


alpha_4_list   = []   # Exit diffuser absolute flow angle (relative to radial)   [degree]
beta_4_list    = []   # Exit diffuser relative flow angle (relative to radial)   [degree]


Tau_list   = []   # Work coefficient        [-] 
Pot_list   = []   # Power                   [kW]
Eta_list   = []   # Isentropic efficiency   [-]

#=======================================================================#
# Input
#=======================================================================#

i              = 1    # Stage index           [-]
number_stage   = 8    # Numer of all stages   [-]


p_01       = 1.03*10**5                                                  # Inlet impeller total pressure      [Pa]
T_01       = 293.00                                                      # Inlet impeller total temperature   [K]


while i <= number_stage:


    gas        = 'AIR'                                                       # AIR|CO2|CH4|H2|N2
    if i == 1:
            pass
    else:
        p_01       = p_04                                                    # Inlet impeller total pressure                              [Pa]
        T_01       = T_04                                                    # Inlet impeller total temperature                           [K]
    m_f        = 25                                                          # Mass flow rate                                             [kg/s]
    m_w_air    = 0.028965                                                    # Molar mass of air                                          [kg/mol]
    c_p        = 1006.2                                                      # Specific heat capacity at constant pressure                [J/kg K]
    c_v        = 719.15                                                      # Specific heat capacity at constant volume                  [J/kg K]
    R          = 8.314                                                       # Universal gas constant                                     [J/kg K]
    R_specific = R / m_w_air                                                 # Specific gas constant of air                               [-]
    k          = c_p/c_v                                                     # Isentropic expansion factor                                [-]
    rho_01     = p_01/ (R_specific* T_01)                                    # Inlet impeller total density                               [kg/m^3]
    Eta_is_n   = 0.90                                                        # Value of first attempt for isentropic efficiency           [-]
    beta_tt    = 80**(1/number_stage)                                        # Total pressure ratio of one stage                          [-]
    N_s        = 0.12                                                        # Specific speed                                             [-]
    Mu         = 1.789*10**(-5)                                              # Dynamic viscocity                                          [Pa.s]
    h_01       = c_p*T_01                                                    # Inlet impeller total enthalpy                              [J/kg]
    phi_2      = 0.22                                                        # Exit flow coefficient                                      [-]
    beta_2_inf = 50
    t          = 0.004                                                       # Impeller thickness                                         [m]


#####################FIRST STAGE#########################################
#=======================================================================#
# Impeller Section
#=======================================================================#
    Ext_loss   = 0   # External loss   [J/kg]
    
    
    Error_eta   = 1   # Isentropic efficiency error   [-]
    while Error_eta > 0.001:
        
        
        Eta_is     = Eta_is_n                                                    # Isentropic efficiency
        delta_h_0  = c_p* T_01 * ((beta_tt**((k - 1) / k)) - 1) / Eta_is_n       # Impeller total enthalpy drop                               [J/kg]
        delta_T_0  = delta_h_0/c_p                                               # Impeller total temperature drop                            [K]
        T_02       = T_01+delta_T_0                                              # Exit impeller total temperature                            [K]
        P_t        = m_f*delta_h_0/1000                                          # Total Power                                                [kW]
        Z          = int((90.0 - beta_2_inf) / 3.0)+2                            # Number of Blades                                           [-]
        Q_1        = m_f/rho_01                                                  # Volumetric flow rate at inlet impeller                     [m^3/s]
        if i == 1:
            N          = N_s * ( (delta_h_0 * Eta_is) ** 0.75 ) / (Q_1 ** 0.5)       # Rotational speed                                           [round/s]
            RPM        = N*60                                                        # Rotational Speed                                           [round/min]
            Omega      = 2*pi*N                                                      # Rotational Speed                                           [rad/s]
            slip_f     = 1 - (cos(radians(beta_2_inf)) ** 0.5) / (Z ** 0.7)          # Slip factor                                                [-]
            tau        = slip_f-phi_2*tan(radians(beta_2_inf))                       # Work coefficient                                           [-]
            # u_2        = (delta_h_0-Ext_loss/tau)**0.5                               # Exit impeller peripheral speed                             [m/s]
            u_2        = (delta_h_0/tau)**0.5                                        # Exit impeller peripheral speed                             [m/s]
            D_2        = (2*u_2/Omega)                                               # Exit impeller diameter                                     [m]
            r_2        = D_2/2                                                       # Exit impeller radius                                       [m]
            cr2        = phi_2*u_2                                                   # Exit impeller meridional component of absolute velocity    [m/s]                                                                                                               [m/s]
            cth2       = tau*u_2                                                     # Exit impeller tangential component of absolute velocity    [m/s]
        else:
            slip_f     = 1 - (cos(radians(beta_2_inf)) ** 0.5) / (Z ** 0.7)          # Slip factor                                                [-]
            # tau        = (delta_h_0-Ext_loss)/(u_2**2)                               # Work coefficient                                           [-]            
            # cr2        = (slip_f-tau)*u_2/tan(radians(beta_2_inf))                   # Exit impeller meridional component of absolute velocity    [m/s]                
            tau        = slip_f-phi_2*tan(radians(beta_2_inf))                       # Work coefficient                                           [-]
            cr2        = phi_2*u_2                                                   # Exit impeller meridional component of absolute velocity    [m/s]                
            cth2       = tau*u_2                                                     # Exit impeller tangential component of absolute velocity    [m/s]
            wth2       = u_2-cth2                                                    # Exit impeller tangential component of relative velocity   [m/s]
        wth2       = u_2-cth2                                                    # Exit impeller tangential component of relative velocity   [m/s]
        w2         = (cr2**2+wth2**2)**0.5                                       # Exit impeller relative velocity                            [m/s]
        c_2        = (cr2**2+cth2**2)**0.5                                       # Exit impeller absolute velocity                            [m/s]
        beta2      = atan(wth2/cr2)                                              # Exit impeller relative flow angle (relative to radial)     [rad]
        alpha2     = atan(cth2/cr2)                                              # Exit impeller absolute flow angle (relative to radial)     [rad]
        a_02       = (k*R_specific*T_02)**0.5                                    # Exit impeller total speed of sound                         [m/s]
        a_01       = (k*R_specific*T_01)**0.5                                    # Inlet impeller total speed of sound                        [m/s] 
        M_u        = u_2/a_01                                                    # Peripheral Mach number                                     [-]
        
        
        if i == 1:
            Si_a     = 235.00/5                                  # Allowable stress                                                                    [MPa]
            Tau_s    = Si_a/((3)**(1/2))                         # Allowable torsional stress                                                          [MPa]
            D_1htt   = ((16*P_t)/(math.pi*Omega*Tau_s))**(1/3)   # Minimum inlet impeller hub diameter based on torque trasmission (t = trasmission)   [mm]
            
    
# Shaft diameter based on torsional deformation  

            G        = 220/(2*(1+0.30))                              # Shear module                                                                            [GPa]
            D_1htd   = ((32*P_t)/(math.pi*Omega*G*(1/0.25)))*(1/4)   # Minimum inlet impeller hub diameter based on torsionale deformation (d = deformation)   [mm]
      
    
# Shaft diameter based on preliminary criterion   

            D_1hpc   = 0.2*D_2   # Mimimum inlet impeller hub diameter based on preliminary criterion (pc = preliminary criterion)   [m]   
        
        
            D1h    = max(D_1htt*0.001,D_1htd*0.001,D_1hpc)   # Inlet impeller hub diameter   [m]                             
            r_1_h  = D1h/2                                   # Hub Radius at the Inlet       [m]
            
            
        else:
            pass
            
        
        if i == 1:
                M_1_t      = 4                                                          # Mach Number of the Inlet at the Tip                        [-]
                r_1_t_n    = r_1_h+0.15                                                 # Value of first attempt for impeller tip radius             [m]
                u_1_h      = Omega*r_1_h                                                # Peripheral Speed at the hub                                [m/s]
    
        
                while r_1_t_n > r_1_h:
                        rho_1n     = p_01/(T_01*R_specific)                              # Value of first attempt fo inlet impeller density           [kg/m^3]
                        Error_d    = 1                                                   # Inlet impeller density error                               [-]
                        r_1_t      = r_1_t_n                                             # Inlet impeller tip radius                                  [m] 
                        while Error_d > 0.001:
                                rho_1     = rho_1n                                       # Inlet impeller density                                     [kg/m^3]
                                A1        = pi*((r_1_t**2)-(r_1_h**2))                   # Inlet impeller area                                        [m^2]   
                                C_1_m     = m_f / (rho_1 * A1)                           # Inlet impeller meridional component of absolute velocity   [m/s] 
                                C_1       = C_1_m                                        # Inlet impeller absolute velocity                           [m/s]
                                h_1       = h_01-(0.5*C_1**2)                            # Inlet impeller enthalpy                                    [J/kg]
                                T_1       = h_1/c_p                                      # Inlet impeller temperature                                 [K]
                                p_1       = p_01-0.5*rho_1*C_1**2                        # Inlet impeller pressure                                    [Pa]
                                rho_1n    = p_1/(T_1*R_specific)                         # New value of inlet impeller density                        [kg/m^3]
                                Error_d   = (abs(rho_1n-rho_1)/rho_1)                    # Inlet impeller density error                               [-] 
                        u_1_t    = Omega*r_1_t                                           # Inlet impeller tip peripheral speed                        [m/s]
                        w_1_t    = sqrt(C_1**2+u_1_t**2)                                 # Inlet impeller tip relative velocity                       [m/s]
                        a_1_t    = sqrt(k * R_specific * T_01)                           # Inlet impeller tip speed of sound                          [m/s]
                        M_1_t_n  = w_1_t / a_1_t                                         # Inlet impeller tip relative Mach number                    [-]
                        if M_1_t_n < M_1_t:
                                r_1_t_n   = r_1_t+0.001                                  # New value of tip radius                                    [m]
                                M_1_t     = M_1_t_n                                      # Old value of inlet impeller tip relative Mach number       [-]
                        else:
                                break
                            
    
        else:
                rho_1n     = p_01/(T_01*R_specific)                              # Value of first attempt fo inlet impeller density           [kg/m^3]
                Error_d    = 1                                                   # Inlet impeller density error                               [-]
                while Error_d > 0.001:
                        rho_1     = rho_1n                                       # Inlet impeller density                                     [kg/m^3]
                        A1        = pi*((r_1_t**2)-(r_1_h**2))                   # Inlet impeller area                                        [m^2]   
                        C_1_m     = m_f / (rho_1 * A1)                           # Inlet impeller meridional component of absolute velocity   [m/s] 
                        C_1       = C_1_m                                        # Inlet impeller absolute velocity                           [m/s]
                        h_1       = h_01-(0.5*C_1**2)                            # Inlet impeller enthalpy                                    [J/kg]
                        T_1       = h_1/c_p                                      # Inlet impeller temperature                                 [K]
                        p_1       = p_01-0.5*rho_1*C_1**2                        # Inlet impeller pressure                                    [Pa]
                        rho_1n    = p_1/(T_1*R_specific)                         # New value of inlet impeller density                        [kg/m^3]
                        Error_d   = (abs(rho_1n-rho_1)/rho_1)                    # Inlet impeller density error                               [-] 
                u_1_t    = Omega*r_1_t                                           # Inlet impeller tip peripheral speed                        [m/s]
                w_1_t    = sqrt(C_1**2+u_1_t**2)                                 # Inlet impeller tip relative velocity                       [m/s]
             

        w_1_h      = sqrt(u_1_h**2+C_1_m**2)   # Inlet impeller hub relative velocity     [m/s]
        beta_1_h   = math.acos(C_1_m/w_1_h)    # Inlet impeller hub relative flow angle   [rad.]
        beta_1_t   = math.acos(C_1_m/w_1_t)    # Inlet impeller tip relative flow angle   [rad.]
        
        
        D_1_t   = 2*r_1_t   # Inlet impeller tip diameter   [m]      
        
        
        Eta_p         = 0.95
        n             = Eta_p*k/(k*(Eta_p-1)+1)
        beta_tt_imp   = ((delta_h_0+Ext_loss)/(c_p*T_01)+1)**(n/(n-1))   # Impeller total to total pressure ratio   [-]   
        p_02          = p_01*beta_tt_imp                                 # Exit impeller total pressure             [Pa]
        T_2           = T_02-c_2**2/(2*c_p)                              # Exit impeller temperature                [K]
        p_2           = p_02*(T_2/T_02)**(k/(k-1))                       # Exit impeller pressure                   [Pa]
        rho_2         = p_2/(R_specific*T_2)                             # Exit impeller density                    [kg/m^3]
        
        
        rho_02   = p_02/(R_specific*T_02)                     # Exit impeller total density              [kg/m^3]
        Q_2      = m_f/rho_2                                  # Exit impeller volumetric rlow rate       [m^3/s]
        f        = 1-Z*t/(pi*D_2*cos(radians(beta_2_inf)))    # Impeller incumberence coefficient        [-]
        b2       = Q_2/(f*pi*D_2*cr2)                         # Impeller width                           [m]
 
    
#=======================================================================#
# Diffuser Section
#=======================================================================#
    
    
        b_4       = b2          # Exit diffuser passage width                       [mm]
        cp        = 0.65        # Pressure recovery coefficient                     [-]


# Physical constants                               
        if i == 1:
                r_4n      = r_2+0.001   # Value of first attempt for exit diffuser radius   [m]
                while r_4n > r_2:
                        r_4   = r_4n   # Exit diffuser radius   [m]
                        
                        
                        Error_d   = 1.00   # Exit diffuser density error                        [-]                
                        rho_4n   = rho_2   # Value of first attemot for exit diffuser density   [kg/m^3]
                        while Error_d > 0.001:
                                rho_4   = rho_4n
                                
                                
                                A_4     = pi * 2 * r_4 * b_4   # Flow area at diffuser exit   [m^2]
                                
                                
                                c_4r    = m_f/(rho_4*A_4)          # Exit diffuser radial component of absolute velocity       [m/s]                                       
                                c_4th   = (r_2/r_4)*cth2           # Exit diffuser tangential component of absolute velocoty   [m/s] 
                                c_4     = sqrt(c_4th**2+c_4r**2)   # Exit diffuser absolute velocity                           [m/s]                                                                              
                                alpha4  = atan(c_4th/c_4r)         # Exit diffuser absolute flow angle (relative to radial)    [rad.]
                                h_02    = h_01+delta_h_0           # Exit impeller total enthalpy                              [J/kg]
                                h_04    = h_02                     # Exit diffuser total enthalpy                              [J/kg]
                                h_4     = h_04-0.5*c_4**2          # Exit diffuser enthalpy                                    [J/kg]
                                
                                
                                T_04    = h_04/c_p                 # Exit diffuser total temperature   [K]
                                T_4     = T_04-c_4**2/(2*c_p)      # Exit diffuser temperature         [Pa]                          
                                p_4     = (p_02-p_2)*cp+p_2        # Exit diffuser pressure            [Pa]
                                p_04    = p_4+0.5*rho_4*c_4**2     # Exit diffuser total pressure      [Pa]

                                
                                rho_4n     = p_4/(R_specific*T_4)        # New value of exit diffuser density   [kg/m^3]                                       
                                Error_d    = abs(rho_4n-rho_4)/rho_4     # Exit diffuser density error          [-] 
                                
                        
                        if c_4**2 <= 0.2*c_2**2 or 0.3*c_2**2 <= c_4**2:    
                                r_4n   = r_4+0.001   # New value of exit diffuser raidus   [m]
                        else:
                                break


        else:
                Error_d   = 1.00   # Exit diffuser density error                        [-]                
                rho_4n   = rho_2   # Value of first attemot for exit diffuser density   [kg/m^3]
                while Error_d > 0.001:
                        rho_4   = rho_4n
                        
                        
                        A_4     = pi * 2 * r_4 * b_4   # Flow area at diffuser exit   [m^2]
                        
                        
                        c_4r    = m_f/(rho_4*A_4)          # Exit diffuser radial component of absolute velocity       [m/s]                                       
                        c_4th   = (r_2/r_4)*cth2           # Exit diffuser tangential component of absolute velocoty   [m/s] 
                        c_4     = sqrt(c_4th**2+c_4r**2)   # Exit diffuser absolute velocity                           [m/s]                                                                              
                        alpha4  = atan(c_4th/c_4r)         # Exit diffuser absolute flow angle (relative to radial)    [rad.]
                        h_02    = h_01+delta_h_0           # Exit impeller total enthalpy                              [J/kg]
                        h_04    = h_02                     # Exit diffuser total enthalpy                              [J/kg]
                        h_4     = h_04-0.5*c_4**2          # Exit diffuser enthalpy                                    [J/kg]
                        
                        
                        T_04    = h_04/c_p                 # Exit diffuser total temperature   [K]
                        T_4     = T_04-c_4**2/(2*c_p)      # Exit diffuser temperature         [Pa]                          
                        p_4     = (p_02-p_2)*cp+p_2        # Exit diffuser pressure            [Pa]
                        p_04    = p_4+0.5*rho_4*c_4**2     # Exit diffuser total pressure      [Pa]

                        
                        rho_4n     = p_4/(R_specific*T_4)        # New value of exit diffuser density   [kg/m^3]                                       
                        Error_d    = abs(rho_4n-rho_4)/rho_4     # Exit diffuser density error          [-] 
                        
                
#=======================================================================#
# Impeller Losses
#=======================================================================#
    
    
        D_h   = (2 * (r_4 - r_2))      # Hydraulic diameter of annular diffuser passage   [mm]
        Re    = (rho_2 * c_2*D_h)/Mu   # Diffuser exit Reynolds number                    [-]
        Cf    = 0.079 / (Re ** 0.25)   # Skin friction coefficient                        [-]
    
    
# Impeller Blade Loading 
        Zspl        = 0.5*Z                                                                                                    # Number of splitter blades                                                [-]
        Neff        = Z+0.75*Zspl                                                                                              # Number of effective blades                                               [-]
        w1          = sqrt((w_1_t**2+w_1_h**2)/2)                                                                              # Inlet impeller root-mean square relative velocity                        [m/s]
        Df          = 1-(w2/w1)+((delta_h_0/u_2**2)/((w_1_t/w2)*(Neff/pi)*(1-(D_1_t/D_2))+(2*(D_1_t/D_2))))                    # Diffusion factor                                                         [-]
        delta_hbl   = 0.05*Df**2*u_2**2                                                                                        # Impeller blade loading loss                                              [J/kg]   
    
#Impeller Clearance
        delta_r     = 0.003*D_2                                                                                                # Impeller clearance                                                       [m]                                                                  
        delta_hcl   = 0.6*(delta_r/b2)*cth2*sqrt((4*pi*(r_1_t**2-r_1_h**2)*cth2*C_1_m)/(b2*Neff*(r_2-r_1_t)*(1+rho_2/rho_1)))  # Impeller clearance loss                                                  [J/kg]


#Impeller Mixing
        epsilon     = 0.2                                                                                                      # Mixing loss coefficient                                                  [-]                                     
        delta_hmx   = ((0.5*c_2**2)/(1+tan(alpha2)**2))*(epsilon/(1-epsilon)**2)                                               # Impeller mixing loss                                                     [J/kg]
        
    
#Impeller Friction
        Lax        = D_2*(0.014+0.023*D_2/D1h+2.012*Q_1/u_2*D_2**2)                                                            # Axial impeller length                                                    [m]
        beta_1_t   = atan(u_1_t/C_1_m)                                                                                         # Inlet impeller tip relative flow angle (relative to radial)              [rad.]
        beta_1_h   = atan(u_1_h/C_1_m)                                                                                         # Inlet impeller hub relative flow angle (relative to radial)              [rad.]
        Lb         = pi/8*(2*r_2-(r_1_t-r_1_h)-b2+2*Lax)*(2/(cos(beta_1_t)+cos(beta_1_t)))/2+cos(radians(beta_2_inf))          # Impeller blade lenght                                                    [m]  
        D1         = sqrt((D_1_t**2+D1h**2)/2)                                                                                 # Inlet impeller root mean square diameter                                 [m]
        Dh_1       = pi*(D_1_t**2-D1h**2)/(2*pi*D1)+Neff*(D_1_t-D1h)                                                           # Inlet impellet hub diameter                                              [m]
        c_th1      = sqrt(C_1**2-C_1_m**2)                                                                                     # Inlet impeller tangential component of absolute velocity                 [m/s]
        w_bar      = c_th1+w_1_t+c_2+2*w_1_h+3*w2/8                                                                            # Average velocity                                                         [m/s]
        deltah_sf  = 2*Cf*(w_bar**2)*Lb/Dh_1                                                                                   # Impeller friction loss                                                   [J/kg]
   
    
#Impeller Recirculation
        deltah_rc   = 8*10**(-5)*sinh(3.5*alpha2**3)*((Df)**2)*((u_2)**2)                                                      # Impeller recirculation loss                                              [J/kg]

    
#Impeller Disc Friction
        rho_bar       =  (rho_1+rho_2)/2                                                                                       # Impeller average density                                                 [kg/m^3]
        deltah_disk   = Cf*(rho_bar*r_2**2*u_2**3)/(4*m_f)                                                                     # Impeller disk friction loss                                              [J/kg]
    
    
#Impeller Leakage
        b1          = r_1_t-r_1_h                                                                                               # Inlet impeller blade width                                              [m]
        b_bar       = (b1+b2)/2                                                                                                 # Impeller average blade width                                            [m]
        r_bar       = (r_2+r_1_t)/2                                                                                             # Impeller average radius                                                 [m]
        delta_PL    = m_f*(r_2*cth2)/Neff*r_bar*b_bar*Lax                                                                       # Pressure difference at impeller tip                                     [Pa]
        u_l         = 0.816*sqrt(2*delta_PL/rho_2)                                                                              # Leakage peripheral velocity                                             [m/s]
        delta_rad   = 0.001*D_2                                                                                                 # Leakage clearance                                                       [m]
        m_fl        = rho_2*u_l*delta_rad*Lax*Neff                                                                              # Leakage mass flow rate                                                  [kg/s]
        deltah_lk   = m_fl*u_l*u_2/2*m_f                                                                                        # Impeller leakage loss                                                   [J/kg]
        
    
#=======================================================================#
# Diffuser Losses
#=======================================================================#


#Diffuser Friction
        deltah_df = Cf*r_2*(1-(r_2/r_4)**1.5)*c_2**2/1.5*b2*cos(alpha2)                                                         # Diffuser friction loss                                                   [J/kg]
                        

#=======================================================================#
# Return Channel Losses
#=======================================================================# 


# Return Channel Incidence
        alpha_4_opt     = radians(88)                                                                                            # Inlet return channel optimum angle (relative to radial)                 [rad.]            
        delta_h_0_inc   = 0.6*(c_4**2)*sin(alpha4-alpha_4_opt)**2                                                                # Return channel incidence loss                                           [J/kg]

                    
#=======================================================================#
# Real Efficiency with Losses
#=======================================================================#        

    
        Ext_loss   = deltah_rc+deltah_disk+deltah_lk   # External loss   [J/kg]
    
    
        if i == number_stage:
                delta_mer = 0.5*c_4r**2   # Volute meridional loss [J/kg]
        else:
                pass


        delta_h_0_id   = c_p * T_01 * ((beta_tt**((k - 1) / k)) - 1) / 1                                                                     # Ideal impeller total enthalpy drop   [J/kg]
        if i < number_stage:
                tot_loss       = delta_hbl+delta_hcl+deltah_df+delta_hmx+deltah_sf+deltah_rc+deltah_lk+deltah_disk+delta_h_0_inc             # Total loss                           [J/kg]
        else:
                tot_loss       = delta_hbl+delta_hcl+deltah_df+delta_hmx+deltah_sf+deltah_rc+deltah_lk+deltah_disk+delta_h_0_inc+delta_mer   # Total loss                           [J/kg]
        delta_h_real   = delta_h_0_id+tot_loss                                                                                               # Real impeller total enthalpy drop    [J/kg]
        Eta_is_n       = delta_h_0_id/delta_h_real                                                                                           # New value of isentropic efficiency   [-]
        Error_eta      = abs(Eta_is_n-Eta_is)/Eta_is                                                                                         # Isentropic efficiency error          [-]
    print(Eta_is_n)
    

#===============================RESULTS=================================#


    stage_list.append(i)
    
    
    p_01_list.append(p_01)
    T_01_list.append(T_01)
    h_01_list.append(h_01)
    rho_01_list.append(rho_01)
    p_1_list.append(p_1)
    T_1_list.append(T_1)
    h_1_list.append(h_1)
    rho_1_list.append(rho_1)
    
    
    p_02_list.append(p_02)
    T_02_list.append(T_02)
    h_02_list.append(h_02)
    rho_02_list.append(rho_02)
    p_2_list.append(p_2)
    T_2_list.append(T_2)
    h_2   = T_2/c_p   
    h_2_list.append(h_2)
    rho_2_list.append(rho_2)
    
    
    T_04_list.append(T_04)
    h_04_list.append(h_04)
    p_4_list.append(p_4)
    T_4_list.append(T_4)
    h_4_list.append(h_4)
    rho_4_list.append(rho_4)
    
    
    p_04_list.append(p_04)


    c_1_list.append(C_1)
    w_1_h_list.append(w_1_h)
    w_1_t_list.append(w_1_t)
    u_1_h_list.append(u_1_h)
    u_1_t_list.append(u_1_t)
    
    
    c_2_list.append(c_2)
    w_2_list.append(w2)
    u_2_list.append(u_2)
    
    
    c_4_list.append(c_4)
    
    
    alpha_1_list.append(0)    
    beta_1_h_list.append(degrees(beta_1_h))   
    beta_1_t_list.append(degrees(beta_1_t))  
    
    
    alpha_2_list.append(degrees(alpha2))
    beta_2_list.append(degrees(beta2))
    
    
    alpha_4_list.append(degrees(alpha4))
    
    
    Tau_list.append(tau)
    Pot_list.append(P_t)
    Eta_list.append(Eta_is)

    i = i+1   # Stage index update   [-]



#================================PLOT==================================#


plt.plot(stage_list,Eta_list,marker='o',linestyle='-',color='blue',markersize=5,linewidth=1.5)
title_font = {'family':'sans-serif','color':'blue','weight':'normal','size':13}
label_font = {'family':'sans-serif','color':'black','weight':'normal','size':11}
plt.title(r'$Isentropic~Efficiency$',fontdict=title_font)
plt.xlabel('$Stage$ [-]',fontdict=label_font)
plt.ylabel(r'$Eta_{is} [-]$',fontdict=label_font)
plt.grid(True)
plt.show()



























    


