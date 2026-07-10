# Preliminary Design of a Centrifugal Compressor Stage
import math
import numpy as np
from math import asin,atan,degrees,exp,log,pi,radians,sin,tan,sqrt,cos,sinh
#=======================================================================#
# Input
#=======================================================================#
gas       = 'AIR'             # AIR|CO2|CH4|H2|N2
p_01      = 1.03*10**5                                                        # Inlet total pressure		                                 [Pa]
T_01      = 293.00		                                                # Inlet total temperature                                    [K]
m_f       = 25                                                          # Mass Flow Rate                                             [kg/s]
m_w_air   = 0.028965                                                    # Molar Mass of Air                                          [kg/mol]
c_p   = 1006.2                                                          # Specific Heat Capacity at Constant Pressure                [J/kg K]
c_v   = 719.15                                                         # Specific Heat Capacity at Constant Volume                  [J/kg K]
R         = 8.314                                                       # Universal Gas Constant of Air                              [J/kg K]
R_specific = R / m_w_air                                                # Specific Gas Constant of Air                               [-]
k         = c_p/c_v                                                     # Isentropic Expansion Factor                                [-]
rho_01    = p_01/ (R_specific* T_01)                                             # Density                                                    [kg/m^3]
Eta_is_n    = 0.8                                                       # Isentropic Efficiency                                      [-]
beta        = 80                                                        # Total pressure ratio                                       [Pa]
beta_tt   = 80**0.125                                                   # Total pressure ratio                                       [-]
N_s       = 0.12                                                        # Specific Speed                                             [-]
Mu        = 1.789*10**(-5)                                              # Dynamic Viscocity                                          [Pa.s]
h_01      = c_p*T_01                                                    # Stagnation Enthalpy at the Inlet                           [J/kg]
p_02 = p_01*beta_tt                                                     # Stagnation Exit Pressure                                   [Pa]
 
#####################FIRST STAGE#########################################
#=======================================================================#
# Impeller Section
#=======================================================================#
Error_eta=1
iteration_eta=0
max_iter=1000
while Error_eta>0.001 and iteration_eta<max_iter:
    #Eta_is=Eta_is_n
    delta_h_0 = c_p* T_01 * ((beta_tt**((k - 1) / k)) - 1) / Eta_is_n       # Enthalpy Change                                            [J/kg]
    delta_T_0 = delta_h_0/c_p                                               # Temperature Change                                         [K]
    Q_1     = m_f/rho_01                                                    # Volumetric Flow at the inlet                               [m^3/s]
    N = N_s * ( (delta_h_0 * Eta_is_n) ** 0.75 ) / (Q_1 ** 0.5)             # Number of Rotations                                        [rads]
    Omega = 2*math.pi*N/60                                                     # Rotational Speed                                           [rads/s]
    P_t        = m_f*delta_h_0                                              # Total Power                                                [KW]
    T_02      = T_01+delta_T_0                                              # Exit Temperature                                           [K]
    beta_2_inf= 50                                                          # Ideal blade flow angle                                     [degree]
    Z = int((90.0 - beta_2_inf) / 3.0)                                        # Number of Blades                                           [-]
    phi_2     = 0.22                                                        # Exit flow coefficient                                      [-]
    tau_inf = 1 - phi_2*tan(beta_2_inf)                                     # Ideal Work Coefficient                                     [-]                      
    slip_f = 1 - (cos(radians(beta_2_inf)) ** 0.5) / (Z ** 0.7)             # Slip factor                                                [-]
    tau       = slip_f-phi_2*tan(beta_2_inf)                                # Work Coefficient                                           [-]
    u_2        = (delta_h_0/tau)**0.5                                       # Exit peripheral speed                                      [m/s]
    D_2       = (2*u_2/Omega)                                               # Impeller Exit Diameter                                     [mm]
    r_2       = D_2*2                                                       # Impeller Exit Radius                                       [mm]
    cr2      = phi_2*u_2                                                    # Velocity perpendicular to the tangential (Radial)          [m/s]                                                                                                               [m/s]
    cth2_inf = u_2-cr2/tan(radians(beta_2_inf))                             # Ideal Velocity in rotational direction (Tangential)        [m/s]
    cth2     = cth2_inf-u_2*(1-slip_f)                                      # Velocity in rotational direction(Tangential)               [m/s]
    wth2     = u_2-cth2                                                     # Tangential Component of Relative Velocity                  [m/s]
    w2       = (cr2**2+wth2**2)**0.5                                        # Relative Velocity                                          [m/s]
    c_2  = (cr2**2+cth2**2)**0.5                                            # Absolute Velocity                                          [m/s]
    beta2    = degrees(atan(cr2/wth2))                                      # Blade Metal Angle                                          [deg]
    alpha2   = degrees(atan(cr2/c_2))                                        # Absolute Flow Angle                                       [deg]
    t          = 4                                                          # Blade thickness                                            [mm]
    a_02       = (k*R_specific*T_02)**0.5                                   # Exit Impeller Sound Speed                                  [m/s]
    a_01       = (k*R_specific*T_01)**0.5                                   # Inlet Impeller Sound Speed                                 [m/s] 
    M_u        = u_2/a_01                                                   # Peripheral Mach Number                                     [-]
    M_1_t      = 1                                                          # Mach Number of the Inlet at the Tip                        [-]
    D_1_h      = 200                                                        # Hub Diameter at the Inlet                                  [mm]              
    r_1_h      = D_1_h*2                                                    # Hub Radius at the Inlet                                    [mm]
    D_1_t      = D_1_h+100                                                  # Tip Diameter at the Inlet                                  [mm]
    r_1_t      = D_1_t*2                                                    # Tip Radius at the Inlet                                    [mm]
    u_1_h      = Omega*r_1_h                                                # Peripheral Speed at the hub                                [m/s]
    while D_1_t > D_1_h:
        rho_1_h    = p_01/T_01*R_specific                                  # Density at the Hub Diameter                               [kg/m^3]
        Error      = 1  
        while Error>0.0001:
            rho_1=rho_1_h
            A1   = pi*((D_1_t**2)-(D_1_h**2))/4*10**(-6)                        # Flow area in the axial direction                           [mm^2]   
            C_1_m = m_f / (rho_1 * A1)                                        # Meridional Velocity                                        [m/s] 
            C_1 = C_1_m                                                        # Absolute velocity at the inlet                             [m/s]
            h_1 = h_01-(0.5*C_1**2)                                             # Static Enthalpy at the inlet                               [J/kg]
            T_1 = h_1/c_p                                                      # Temperature at the Inlet                                   [K]
            p_1   = p_01-0.5*rho_1*C_1**2                                      # Pressure at the inlet                                      [Pa]
            rho_1_h  = p_1/(T_1*R_specific)                                    # Density at the hub diameter                                [kg/m^3]
            Error  = abs(rho_1_h-rho_1)/rho_1
        u_1_t   = Omega*r_1_t                                                  # Peripheral Speed at the tip                                [m/s]
        w_1_t    = sqrt(C_1**2+u_1_t**2)                                       # Relative velocity at the tip                               [m/s]
        a_1_t    = sqrt(k * R_specific * T_01)                                 # Sound at the tip                                           [m/s]
        M_1_t_n= w_1_t / a_1_t                                                 # New Mach Number at the tip                                 [-]
        if M_1_t_n<M_1_t:
            r_1_t_n=(D_1_t/2)+0.001                                            # Radius at the tip
        else:
            break
      # Shaft diameter based on torque transmission
    y_stress   = 235                                                        # Yield Stress                                                  [MPa]
    safety_factor  = 5                                                      # Safety Factor                                                 [-]
    torsion_stress = (y_stress/safety_factor)/(3**0.5)                      # Torsion Stress                                                [MPa]
    D1h_01 = (10**2)*((16*P_t)/(pi*Omega*torsion_stress))**(0.33)           # Shaft diameter calculation                                    [mm]
      # Shaft diameter based on torsional deformation
    E            = 220                                                      # Young's modulus                                               [GPa]
    vi           = 0.3                                                      # Poisson's ratio
    twist_moment = 0.25                                                     # max allowable twist polar moment of inertia                   [deg/m]
    Gf           = E/(2*(1+vi))                                             # Shear Modulus                                                 [-]
    D1h_02 = (10**2)*((32*P_t)/(pi*Omega*twist_moment*Gf))**0.25            # Shaft Diameter calculation                                    [mm]
      # Shaft diameter based on axial load and fatigue
    D1h_03 = 0.2*D_2                                                        # Shaft Diameter Calculation                                    [mm]
    D1h    = int(round(max(D1h_01,D1h_02,D1h_03),0))                        # Final Shaft Diameter Calculation                              [mm]          
    T_2  = T_02-c_2**2/(2*c_p)                                              # Static Exit Temperature                                        [K]
    rho_02 = p_02/(R_specific*T_02)
    p_2=p_02-(1/2)*rho_02*c_2**2 
    rho_2=p_2/(R_specific*T_2)
    # p_2=p_02                  
    # Error=1
    # while Error>0.0001:
    #     p_2=rho_2*R_specific*T_2
    #     rho_2n=abs(p_02-p_2)*2/c_2**2                                      
    #     Error=abs(rho_2n-rho_02)/rho_02
    Q_2      = m_f/rho_2                                                    # Volumetric Flow at the Outlet                                 [m^3/s]
    f        = 1-Z*t/(pi*D_2*sin(radians(beta_2_inf)))                      # Blade Incumberence Coefficient                                [-]
    b2       = Q_2/(f*pi*D_2*cr2)*10**6                                     # Blade Width                                                   [mm]
    # L_x      = D_2-(D_1_t-D_1_h)/2/(2+b2/2)                                 # Impeller Axial Width                                          [mm]
    # rho_n    = rho_01                                                         # Density at the hub                                            [kg/m^3]
    # rho_delta=1                            
    # while rho_delta>0.0001:
    #    p_2=p_02-0.5*rho_n*w2**2 
    #    rho_n=R_specific*T_2/p_2
    #    rho_delta=abs(rho_n-rho_02)/rho_02
    #    print(rho_delta)
    
    #=======================================================================#
    # Diffuser Section
    #=======================================================================#
    # Physical constants                               
    r_3 = r_2                                                           # Diffuser Inlet Diameter                                          [mm]         
    b_3 = b2                                                            # Blade Width                                                      [mm]
    c3 = c_2                                                            # Absolute velocity at the Diffuser Inlet                          [m/s]
    cp   = 0.65                                                         # Pressure recovery coefficient                                    [-]
    error_d   = 1.00                                                    # Exit diffuser density error                                      
    Deltar    = 0.001                                                   # Increment of exit diffuser radius
    r_4new    = r_2+Deltar                                              # Value of first attempt for exit diffuser radius                  [mm]
    while r_4new > r_2:
        r_4   = r_4new                                                      # Exit diffuser radius
    
    
        rho_3      = rho_2                                                 # Inlet diffuser density
        rho_4new   = rho_3                                                 # Value of first attemot for exit diffuser density
        while error_d > 0.0001:
            D_h = (2 * (r_4 - r_2))*1000                                        # Hydraulic diameter of annular diffuser passage                   [mm]
            A_3 = pi * 2 * r_3 * b_3 * 1e-6                                     # Flow area at diffuser exit                                       [m^2]
            Re = (rho_2 * c_2*D_h)/Mu                                           # Reynolds number at diffuser exit                                 [-]
            Cf = 0.079 / (Re ** 0.25)                                           # Friction factor                                                  [-]
            c_3th= cth2*(r_2 / r_3)+ (2*pi*Cf*rho_2*cth2*(r_3**2-r_2*r_3))/m_f  # Tangential velocity at the diffuser inlet                        [m/s]   
            
            rho_4   = rho_4new                                                                                             # Exit diffuser density                                      [kg/m^3]
            b_3     = b2                                                                                                   # Inlet diffuser height                                      [mm]
            b_4     = b_3                                                                                                  # Exit diffuser height                                       [mm]
            A_4     = 2*math.pi*r_4*b_4                                                                                    # Exit diffuser area                                         [m^2]
            c_4r    = m_f/(rho_4*A_4)                                                                                      # Exit diffuser radial component of absolute velocity        [m/s]                                       
            # c_4th   = (-(r_4/r_3)+math.sqrt((r_4/r_3)**2-4*((2*math.pi*Cf*rho_4*(r_3**2-r_3*r_4))/m_f)*c_3th))/(2*c_3th)   # Exit diffuser tangential component of absolute velocity    [m/s]
            c_4th   = (r_4/r_2)*cth2
            c_4     = sqrt(c_4th**2+c_4r**2)                                                                               # Exit diffuser absolute velocity                            [m/s]                                                                              
            alpha4  = degrees(atan(c_4r/c_4th))                                                                            # Exit diffuser flow angle                                   [-]
            h_02    = h_01+delta_h_0                                                                                       # Exit impeller total enthalpy
            h_04    = h_02                                                                                                 # Exit diffuser total enthalpy
            h_4     = h_04-0.5*c_4**2                                                                                      # Exit diffuser enthalpy
            
            
            T_04    = h_04/c_p             # Exit diffuser total temperature
            T_4     = T_04-c_4**2/(2*c_p)   # Exit diffuser temperature                                  
            p_4     = (p_02-p_2)*c_p+p_2   # Exit diffuser pressure
            P_04    = cp*R_specific*T_04   # Exit diffuser stagnation pressure
            
            
            
            rho4_new   = p_4/(R*T_4)                   # New value of exit diffuser density                                        
            error_d    = abs(rho4_new-rho_4)/rho_4   # Exit diffuser density error
    
        if r_4/r_2 < 0.95*c_2/c_4:
            r_4new   = r_4+Deltar   # New value of exit diffuser raidus
        elif r_4/r_2 > 1.05*c_2/c_4:
            r_4new   = r_4+Deltar   # New value of exit diffuser radius
        else:
            break
   

    #=======================================================================#
    # Impeller Losses
    #=======================================================================#
    #Impeller Incidenece Loss
    D1=sqrt((D_1_t**2+D_1_h**2)/2)                                                                                         # Diameter of the impeller inlet                                                                               # Geometric average diameter
    u1=pi*D1*N/60                                                                                                          # Peripheral Speed at the Inlet of the Impeller
    w1=sqrt(C_1**2+u1**2)                                                                                                  # Relative Velocity at the Inlet of the Impeller                            [-]
    beta1=degrees(atan(u1/C_1_m))                                                                                          # Relative Flow Angle                                                       [-]
    beta1_bl=beta1+degrees(5)                                                                                              # Impeller Inlet Blade Angle                                                         [-]
    finc=0.5                                                                                                               # Incidence Friction Loss Coefficient                                       [-]
    delta_hin=finc*w_1_t**2*sin(beta1-beta1_bl)**2                                                                         # Impeller Incidence Loss                                                   [-]
    
    # Impeller Blade Loading 
    Zspl=0.5*Z                                                                                                             # Number of Splitter Blades                                                 [-]
    Neff=Z+0.75*Zspl                                                                                                       # Number of Effective Blades                                                [-]
    Df=1-(w2/w1)+((delta_h_0/u_2**2)/((w_1_t/w2)*(Neff/pi)*(1-(D_1_t/D_2))+(2*(D_1_t/D_2))))                               # Diffuser Loss Coefficient                                                 [-]
    delta_hbl=0.05*Df**2*u_2**2                                                                                            # Impeller Blade Loading
    
    #Impeller Clearance
    delta_r=0.003*D_2                                                                                                      # Clearance Thickness                                                      [mm]                                                                  [mm]
    delta_hcl=0.6*(delta_r/b2)*cth2*sqrt((4*pi*(r_1_t**2-r_1_h**2)*cth2*C_1_m)/(b2*Neff*(r_2-r_1_t)*(1+rho_2/rho_1)))      # Impeller clearance                                                       [-]
    
    #Impeller Mixing
    epsilon=0.2                                                                                                            # Mixing loss coefficient                                                  [-]                                     
    delta_hmx=(0.5*c_2**2/1+tan(alpha2)**2*(epsilon/1-epsilon)**2)                                                          # Impeller Mixing Loss                                                     [-]
    
    
    #Impeller Friction
    Lax=D_2(0.014+0.023*D_2/D_1_h+2.012*Q_1/u_2*D_2**2)                                                                    # Axial Impeller Length                                                    [mm]
    beta_1_t=atan(u_1_t/C_1_m)                                                                                             # Relative Flow angle at the tip                                           [-]
    beta_1_h=atan(u_1_h/C_1_m)                                                                                             # Relative Flow angle at the hub                                           [-]
    beta2_bl=degrees(90)-(degrees(10)+0.5*N)                                                                             # Blade Metal Angle at the impeller outlet                                 [-]
    Lb=pi/8*(2*r_2-(r_1_t-r_1_h)-b2+2*Lax)*(2/(cos(beta_1_t)+cos(beta_1_t)))/2+cos(beta2_bl)                               # Length of Blade                                                          [mm]  
    Dh_1=pi*(D_1_t**2-D_1_h**2)/(2*pi*D1)+Neff*(D_1_t-D_1_h)                                                               # Diameter of the hub                                                      [mm]
    c_th1=sqrt(C_1**2-C_1_m**2)                                                                                            # Velocity in Tangential Direction                                         [m/s]
    w_1_h = sqrt(C_1**2+u_1_h**2)                                                                                      # Relative Velocity at the Hub                                             [m/s]
    w_bar=c_th1+w_1_t+c_2+2*w_1_h+3*w2/8                                                                                    # Normalized Work Parameter                                                [-]
    deltah_sf=2*Cf(w_bar**2)*Lb/Dh_1                                                                                       # Impeller Skin Friction Coefficient                                       [-]
    
    #Impeller Recirculation
    deltah_rc=8*10**(-5)*sinh(3.5*alpha2**3)*(Df*u_2)**2                                                                   # Impeller Recirculation Coefficient                                       [-]
    
    #Impeller Disc Friction
    rho_bar=(rho_1+rho_2)/2                                                                                                # Average Density                                                          [-]
    deltah_disk=Cf*(rho_bar*r_2**2*u_2**3/4*m_f)                                                                           # Impeller Disk Friction Coefficient                                       [-]
    
    #Impeller Leakage
    Dimp_avg=D1+D_2/2                                                                                                      # Diameter of the average of the Impeller                                  [-]
    b=m_f/rho_1*pi*Dimp_avg*C_1_m                                                                                          # Average Blade Width                                                      [mm]
    b1=2*b-b2                                                                                                              # Blade width of Inlet Impeller                                            [mm]
    b_bar=b1+b2/2                                                                                                          # Average of the blade width                                               [mm]
    r_bar=(r_2+r_1_t)/2                                                                                                    # Average radii of the blade                                               [mm]
    delta_PL=m_f*((r_2*cth2)-(r_1_t*u_1_t))/Neff*r_bar*b_bar*Lax                                                           # TotalPressure Loss                                                       [Pa.s]
    u_l=0.816*sqrt(2*delta_PL/rho_2)                                                                                       # Leakage Velocity                                                         [m/s]
    delta_rad=0.001*D_2                                                                                                    # Normalized Radius                                                        [mm]
    m_fl=rho_2*u_l*delta_rad*Lax*Neff                                                                                      # Mass flow leakage                                                        [kg/s]
    deltah_lk=m_fl*u_l*u_2/2*m_f                                                                                           # Impeller Leakage                                                         [-]
    
    #=======================================================================#
    # Diffuser Losses
    #=======================================================================#
    #Diffuser Friction
    deltah_df=Cf*r_2(1-(r_2/r_3)**1.5)*c_2**2/1.5*b2*cos(alpha2)                                                         # Diffuser Friction                                                         [-]
                                            
    #=======================================================================#
    # Real Efficiency with Losses
    #=======================================================================#                                    
    delta_h_0_id= c_p * T_01 * ((beta_tt**((k - 1) / k)) - 1) / 1
    tot_loss = delta_hin+delta_hbl+delta_hcl+deltah_df+delta_hmx+deltah_sf+deltah_rc+deltah_lk+deltah_disk
    delta_h_real=delta_h_0_id+tot_loss
    Eta_is_n=delta_h_0_id/delta_h_real
    Error=abs(Eta_is_n-Eta_is)/Eta_is
    print(Error)
    Eta_is=Eta_is_n

 #=======================================================================#
 # Return Channel Sizing
 #=======================================================================#  
BM6=0.12                                                               # Blockage factour of the return bend
alpha_6_specified =1                 
alpha6=math.degrees.atan(0.32+(1.7*phi_2))                
b_6=b_4*tan(alpha4)/tan(alpha6/(1-BM6))                                 # Blade height at return bend shroud contour exit
if b_6 < b_4:
    b_6 = b_4  # Minimum constraint
    # Recalculate alpha_6 based on the constrained b_6
    alpha_6 = math.degrees(math.atan(b_4 * math.tan(math.radians(alpha4)) / (b_6 * (1 - BM6))))
elif b_6 > 2 * b_4:
    b_6 = 2 * b_4  # Maximum constraint
    # Recalculate alpha_6 based on the constrained b_6
    alpha_6 = math.degrees(math.atan(b_4 * math.tan(math.radians(alpha4)) / (b_6 * (1 - BM6))))
else:
    alpha_6 = alpha_6_specified
R_ch=b_6+b_4/2                                                            # Radius of the return channel
b8=(R_ch/0.8)+b_6                                                         # Passage width of the eye of the next impeller
Aco=R_ch+(b_4+b_6)/2                                                      # Elliptical contour of the axial semi axes
Bco=R_ch+b_4                                                              # Elliptical contour of the radial semi axes
Rc_exs=b8                                                                 # Radius of the curvature of the shroud contour
Rc_exh=2*b8                                                               # Radius of the curvature of the hub contour

   
 #=======================================================================#
 # Return Channel Sizing Losses
 #=======================================================================# 
 # Incidence Losses
alpha_4_opt=degrees(88)
delta_h_0_inc=c_4*sin**2*(abs(alpha4-alpha_4_opt))/2


######################SECOND STAGE#######################################
#=======================================================================#
# Impeller Section
#=======================================================================#
Error_eta=1
while Error_eta>0.001:
    Eta_is=Eta_is_n
    delta_h_0 = c_p * T_04 * ((beta_tt**((k - 1) / k)) - 1) / Eta_is        # Enthalpy Change                                            [J/kg]
    delta_T_0 = delta_h_0/c_p                                               # Temperature Change                                         [K]
    h_04s=h_01+Eta_is_n*(h_04-h_01)
    T_04s=h_04s/c_p
    p_04s=p_01*(T_04s/T_01)**(k/k-1)
    p_04=p_04s-(deltah_df+delta_h_0_inc)
    rho_04=P_04/R_specific*T_04
    Q_1     = m_f/rho_01                                                    # Volumetric Flow at the inlet                               [m^3/s]
    P_t        = m_f*delta_h_0                                              # Total Power                                                [KW]
    T_05      = T_04+delta_T_0                                              # Exit Temperature                                           [K]
    beta_2_inf= 50                                                          # Ideal blade flow angle                                     [degree]
    Z      = (int(90-beta_2_inf)/3)                                         # Number of Blades                                           [-]
    phi_2     = 0.22                                                        # Exit flow coefficient                                      [-]
    tau_inf = 1 - phi_2*tan(beta_2_inf)                                     # Ideal Work Coefficient                                     [-]                      
    slip_f = 1 - (cos(radians(beta_2_inf)) ** 0.5) / (Z ** 0.7)             # Slip factor                                                [-]
    tau       = slip_f-phi_2*tan(beta_2_inf)                                # Work Coefficient                                           [-]
    beta2    = degrees(atan(c_4r/c_4th))                                      # Blade Metal Angle                                          [deg]
    alpha2   = degrees(atan(c_4r/c_4))                                        # Absolute Flow Angle                                       [deg]
    t          = 4                                                          # Blade thickness                                            [mm]
    while D_1_t > D_1_h:
        rho_1_h    = p_04/T_04*R_specific                                   # Density at the Hub Diameter                               [kg/m^3]
        Error      = 1  
        while Error>0.0001:
            rho_1=rho_1_h
            C_1_m= m_w_air/(rho_1*A1)                                          # Meridional Velocity                                        [m/s] 
            C_1 = C_1_m                                                        # Absolute velocity at the inlet                             [m/s]
            h_4 = h_04-(0.5*C_1**2)                                            # Static Enthalpy at the inlet                               [J/kg]
            T_4 = h_4/c_p                                                      # Temperature at the Inlet                                   [K]
            p_4   = p_04-0.5*rho_1*C_1**2                                      # Pressure at the inlet                                      [Pa]
            rho_1_h  = p_4/(T_4*R_specific)                                    # Density at the hub diameter                                [kg/m^3]
            Error  = abs(rho_1_h-rho_1)/rho_1
        a_1_t    = sqrt(k * R_specific * T_04)                                 # Sound at the tip                                           [m/s]
        M_1_t_n= w_1_t / a_1_t                                                 # New Mach Number at the tip                                 [-]
        if M_1_t_n<M_1_t:
            r_1_t_n=(D_1_t/2)+0.001                                            # Radius at the tip
        else:
            break          
    T_5  = T_05-c_2**2/(2*c_p)                                              # Static Exit Temperature                                        [K]
    p_05=p_04*beta_tt
    rho_05 = p_05/(R_specific*T_05)
    p_5=p_05-(1/2)*rho_05*c_2**2 
    rho_5=p_5/(R_specific*T_5)
    p_5=p_05                  
    Error=1
    while Error>0.0001:
        p_5=rho_5*R_specific*T_5
        rho_5n=abs(p_05-p_5)*2/c_2**2                                      
        Error=abs(rho_5n-rho_05)/rho_05
    b2       = Q_2/(f*pi*D_2*cr2)*10**6                                     # Blade Width                                                   [mm]
    L_x      = D_2-(D_1_t-D_1_h)/2/(2+b2/2)                                 # Impeller Axial Width                                          [mm]
    rho_n  = rho_04                                                         # Density at the hub                                            [kg/m^3]
    rho_delta=abs(rho_n-rho_05)/rho_05                         
    while rho_delta>0.0001:
       p_5=p_05-0.5*rho_05*w2**2 
       rho_n=R_specific*T_5/p_5
       rho_delta=abs(rho_n-rho_05)/rho_05

    #=======================================================================#
    # Diffuser Section
    #=======================================================================#
    # Physical constants                               
    Re = (rho_5 * c_2*D_h)/Mu                                           # Reynolds number at diffuser exit                                 [-]
    Cf = 0.079 / (Re ** 0.25)                                           # Friction factor                                                  [-]
    c_3th= cth2*(r_2 / r_3)+ (2*pi*Cf*rho_2*cth2*(r_3**2-r_2*r_3))/m_f  # Tangential velocity at the diffuser inlet                        [m/s]   
    rho_6    = rho_5                                                # Inlet diffuser density
    while error_d > 0.0001:
        c_4th   = (-(r_4/r_3)+math.sqrt((r_4/r_3)**2-4*((2*math.pi*Cf*rho_4*(r_3**2-r_3*r_4))/m_f)*c_3th))/(2*c_3th)   # Exit diffuser tangential component of absolute velocity    [m/s]
        c_4     = sqrt(c_4th**2+c_4r**2)                                                                               # Exit diffuser absolute velocity                            [m/s]                                                                              
        alpha4  = degrees(atan(c_4r/c_4th))                                                                            # Exit diffuser flow angle                                   [-]
        h_05   = h_04+delta_h_0                                                                                       # Exit impeller total enthalpy
        h_07   = h_05                                                                                                 # Exit diffuser total enthalpy
        
        
        T_07    = h_07/c_p             # Exit diffuser total temperature
        T_7     = T_07-c_4**2/(2*c_p)   # Exit diffuser temperature                                  
        p_7     = (p_05-p_5)*c_p+p_5   # Exit diffuser pressure 
        p_07    = p_7*(T_07/T_7)**(k/k-1)
    #=======================================================================#
    # Impeller Losses
    #=======================================================================#
    #Impeller Incidenece Loss
    w1=sqrt(C_1**2+u1**2)                                                                                                  # Relative Velocity at the Inlet of the Impeller                            [-]
    beta1=degrees(atan(u1/C_1_m))                                                                                          # Relative Flow Angle                                                       [-]
    beta1_bl=beta1+degrees(5)                                                                                              # Impeller Inlet Blade Angle                                                [-]
    delta_hin=finc*w_1_t**2*sin(beta1-beta1_bl)**2                                                                         # Impeller Incidence Loss                                                   [-]
    
    # Impeller Blade Loading 
    Df=1-(w2/w1)+((delta_h_0/u_2**2)/((w_1_t/w2)*(Neff/pi)*(1-(D_1_t/D_2))+(2*(D_1_t/D_2))))                               # Diffuser Loss Coefficient                                                  [-]
    delta_hbl=0.05*Df**2*u_2**2                                                                                            # Impeller Blade Loading
    
    #Impeller Clearance                                                               
    delta_hcl=0.6*(delta_r/b2)*cth2*sqrt((4*pi*(r_1_t**2-r_1_h**2)*cth2*C_1_m)/(b2*Neff*(r_2-r_1_t)*(1+rho_5/rho_4)))      # Impeller clearance                                                       [-]
    
    #Impeller Mixing                                     
    delta_hmx=(0.5*c_2**2/1+tan(alpha2)**2*(epsilon/1-epsilon)**2)                                                          # Impeller Mixing Loss                                                     [-]
    
    
    #Impeller Friction
    Lax=D_2(0.014+0.023*D_2/D_1_h+2.012*Q_1/u_2*D_2**2)                                                                    # Axial Impeller Length                                                    [mm]
    beta_1_t=atan(u_1_t/C_1_m)                                                                                             # Relative Flow angle at the tip                                           [-]
    beta_1_h=atan(u_1_h/C_1_m)                                                                                             # Relative Flow angle at the hub                                           [-]
    beta2_bl=degrees(90)-(degrees(10)+0.5*N)                                                                             # Blade Metal Angle at the impeller outlet                                   [-]
    Lb=pi/8*(2*r_2-(r_1_t-r_1_h)-b2+2*Lax)*(2/(cos(beta_1_t)+cos(beta_1_t)))/2+cos(beta2_bl)                               # Length of Blade                                                          [mm]  
    c_th1=sqrt(C_1**2-C_1_m**2)                                                                                            # Velocity in Tangential Direction                                         [m/s]
    w_1_h = sqrt(C_1**2+u_1_h**2)                                                                                      # Relative Velocity at the Hub                                                 [m/s]
    w_bar=c_th1+w_1_t+c_2+2*w_1_h+3*w2/8                                                                                    # Normalized Work Parameter                                               [-]
    deltah_sf=2*Cf(w_bar**2)*Lb/Dh_1                                                                                       # Impeller Skin Friction Coefficient                                       [-]
    
    #Impeller Recirculation
    deltah_rc=8*10**(-5)*sinh(3.5*alpha2**3)*(Df*u_2)**2                                                                   # Impeller Recirculation Coefficient                                       [-]
    
    #Impeller Disc Friction
    rho_bar=(rho_4+rho_5)/2                                                                                                # Average Density                                                          [-]
    deltah_disk=Cf*(rho_bar*r_2**2*u_2**3/4*m_f)                                                                           # Impeller Disk Friction Coefficient                                       [-]
    
    #Impeller Leakage
    delta_PL=m_f*((r_2*cth2)-(r_1_t*u_1_t))/Neff*r_bar*b_bar*Lax                                                           # TotalPressure Loss                                                       [Pa.s]
    deltah_lk=m_fl*u_l*u_2/2*m_f                                                                                           # Impeller Leakage                                                         [-]
    
    #=======================================================================#
    # Diffuser Losses
    #=======================================================================#
    #Diffuser Friction
    deltah_df=Cf*r_2(1-(r_2/r_3)**1.5)*c_2**2/1.5*b2*cos(alpha2)                                                          # Diffuser Friction                                                         [-]
                                            
    #=======================================================================#
    # Real Efficiency with Losses
    #=======================================================================#                                    
    delta_h_0_id= c_p * T_04 * ((beta_tt**((k - 1) / k)) - 1) / 1
    tot_loss = delta_hin+delta_hbl+delta_hcl+deltah_df+delta_hmx+deltah_sf+deltah_rc+deltah_lk+deltah_disk
    delta_h_real=delta_h_0_id+tot_loss
    Eta_is_n=delta_h_0_id/delta_h_real
    Error_eta=abs(Eta_is_n-Eta_is)/Eta_is
    Eta_is=Eta_is_n

#=======================================================================#
# Volute Design for Second Stage (After Diffuser)
#=======================================================================#
r_3_volute = r_4                           # base circumference
c_th2_volute = c_4th                       # tangential velocity entering the volute
Q_volute = m_f /  rho_6                      # volumetric flow rate entering the volute
alpha_2_volute = alpha4                     # flow angle entering the volute

# Volute parameters
t_v = 0.003  # Tongue thickness (min 3 mm)
k_b = 1.03   

# Base circumference
r_3 = k_b * r_2  

# Initialize volute table
volute_table = []
angles = [math.pi/4, math.pi/2, 3*math.pi/4, math.pi,
          5*math.pi/4, 3*math.pi/2, 7*math.pi/4, 2*math.pi]

# Calculate volute parameters for each azimuthal angle
for theta_i in angles:
    # Solve for A_t using conservation of angular momentum
    # Equation: A_t * (r_2 * c_th2_volute / Q_volute) - sqrt(A_t/π) - (r_3 + t_v) = 0
   
    # Define the function to solve
    def volute_area_eq(A_t):
        return A_t * (r_2 * c_th2_volute / Q_volute) - math.sqrt(A_t/math.pi) - (r_3 + t_v)
   
    # Solve using Newton-Raphson method
    A_t_guess = Q_volute * theta_i / (2 * math.pi * c_th2_volute)
    tolerance = 1e-6
    max_iter = 100
   
    for iter in range(max_iter):
        f_val = volute_area_eq(A_t_guess)
        f_deriv = (r_2 * c_th2_volute / Q_volute) - (1/(2*math.sqrt(math.pi*A_t_guess)))
       
        # Avoid division by zero
        if abs(f_deriv) < 1e-10:
            break
           
        A_t_new = A_t_guess - f_val / f_deriv
       
        if abs(A_t_new - A_t_guess) < tolerance:
            A_t_guess = A_t_new
            break
           
        A_t_guess = A_t_new
   
    # Calculate section radius for circular cross-section
    r_si = math.sqrt(A_t_guess / math.pi)
   
    # Calculate local external radius
    r_i = r_3 * math.exp(theta_i * math.tan(math.radians(alpha_2_volute)))
   
    # Calculate velocity at this section
    c_i = Q_volute / A_t_guess
   
    # Store results
    volute_table.append({
        'theta_i': math.degrees(theta_i),
        'r_i': r_i,
        'A_i': A_t_guess,
        'r_si': r_si,
        'c_i': c_i
    })

# Print volute design table
print("\nVolute Design Parameters for Second Stage:")
print("θ_i [deg] | r_i [m] | A_i [m²] | r_si [m] | c_i [m/s]")
print("-" * 55)
for row in volute_table:
    print(f"{row['theta_i']:8.1f} | {row['r_i']:7.4f} | {row['A_i']:8.6f} | {row['r_si']:7.4f} | {row['c_i']:7.2f}")
# =======================================================================#
# Calculate Volute Losses
# =======================================================================#
delta_mv=rho_6*(c_4r**2)/2                                     #meridional velocity dump losses at the impeller exit                                                                                                                                                                                                                                                                                                        
delta_h_volute =  delta_mv/ rho_6





#=======================================================================#
# Print Results - INCLUDING REQUESTED PARAMETERS
#=======================================================================#
print(f"Mass Flow Rate: {m_f:.2f} kg/s")
print(f"Inlet Pressure: {p_01:.2f} bar")
print(f"Inlet Temperature: {T_01:.2f} K")
print(f"Pressure Ratio: {beta_tt:.3f}")
print(f"Rotational Speed: {N:.2f} rad/s ({N/(2*math.pi)*60:.0f} RPM)")
print(f"Impeller Exit Diameter: {D_2:.2f} mm")
print(f"Impeller Exit Width: {b2:.2f} mm")
print(f"Peripheral Speed: {u_2:.2f} m/s")
print(f"Absolute Velocity: {c_2:.2f} m/s")
print(f"Relative Velocity: {w2:.2f} m/s")
print(f"Blade Metal Angle: {beta2:.2f}°")
print(f"Absolute Flow Angle: {alpha2:.2f}°")
print(f"Isentropic Efficiency: {Eta_is_n:.4f}")
print(f"Power: {P_t/1000:.2f} kW")

# NEWLY ADDED REQUESTED PARAMETERS
print("\n--- REQUESTED ADDITIONAL PARAMETERS ---")
print(f"Tip Mach Number: {M_1_t_n:.4f}")
print(f"Tip Radius: {r_1_t:.4f} mm")
print(f"Blade Angle at Inlet: {beta1_bl:.2f}°")
print(f"Work Coefficient (τ): {tau:.4f}")
print(f"Ideal Work Coefficient (τ_inf): {tau_inf:.4f}")
print(f"Slip Factor: {slip_f:.4f}")

print(f"Diffuser Exit Radius: {r_4:.2f} mm")
print(f"Return Channel Height: {b_6:.2f} mm")
print("="*80)

# Store key variables for second stage
p_04 = P_04/1e5  # [bar]
T_04 = T_04      # [K]  
h_04 = h_04      # [J/kg]
rho_04 = rho_4   # [kg/m³]

print(f"\nExit conditions for second stage:")
print(f"Pressure: {p_04:.3f} bar")
print(f"Temperature: {T_04:.2f} K")
print(f"Enthalpy: {h_04:.2f} J/kg")

#=======================================================================#
# Additional Detailed Output Section
#=======================================================================#
print("\n" + "="*80)
print("DETAILED PERFORMANCE PARAMETERS")
print("="*80)
print(f"Specific Speed (N_s): {N_s:.3f}")
print(f"Flow Coefficient (φ): {phi_2:.3f}")
print(f"Peripheral Mach Number (M_u): {M_u:.3f}")
print(f"Enthalpy Change: {delta_h_0/1000:.2f} kJ/kg")
print(f"Temperature Rise: {delta_T_0:.2f} K")
print(f"Impeller Axial Length: {L_x:.2f} mm")
print(f"Shaft Diameter: {D1h:.2f} mm")
print(f"Number of Blades: {Z}")
print(f"Blade Thickness: {t} mm")
print("="*80)


# ############################################### COMPRESSOR TWO ###########################################################################################

#=======================================================================#
# Input
#=======================================================================#
p_01      = p_07                                                        # Inlet total pressure		                                 [bar]
T_01      = T_07		                                                # Inlet total temperature                                    [K]
rho_01    =(p_01 * 1e5) / (R_specific * T_01)                           # Density                                                    [kg/m^3]
h_01      = c_p*T_01                                                    # Stagnation Enthalpy at the Inlet                           [J/kg]
p_02 = p_01*beta_tt                                                     # Stagnation Exit Pressure                                   [Pa]

#####################FIRST STAGE#########################################
#=======================================================================#
# Impeller Section
#=======================================================================#
Error_eta=1
while Error_eta>0.001:
    Eta_is=Eta_is_n
    delta_h_0 = c_p * T_01 * ((beta_tt**((k - 1) / k)) - 1) / Eta_is        # Enthalpy Change                                            [J/kg]
    delta_T_0 = delta_h_0/c_p                                               # Temperature Change                                         [K]
    Q_1     = m_f/rho_01                                                    # Volumetric Flow at the inlet                               [m^3/s]
    T_02      = T_01+delta_T_0                                              # Exit Temperature                                           [K]
    u_2        = (delta_h_0/tau)**0.5                                       # Exit peripheral speed                                      [m/s]
    cr2      = phi_2*u_2                                                    # Velocity perpendicular to the tangential (Radial)          [m/s]                                                                                                               [m/s]
    cth2_inf = u_2-cr2/tan(radians(beta_2_inf))                             # Ideal Velocity in rotational direction (Tangential)        [m/s]
    cth2     = cth2_inf-u_2*(1-slip_f)                                      # Velocity in rotational direction(Tangential)               [m/s]
    wth2     = u_2-cth2                                                     # Tangential Component of Relative Velocity                  [m/s]
    w2       = (cr2**2+wth2**2)**0.5                                        # Relative Velocity                                          [m/s]
    c_2  = (cr2**2+cth2**2)**0.5                                            # Absolute Velocity                                          [m/s]
    beta2    = degrees(atan(cr2/wth2))                                      # Blade Metal Angle                                          [deg]
    alpha2   = degrees(atan(cr2/c_2))                                        # Absolute Flow Angle                                       [deg    a_02       = (k*R_specific*T_02)**0.5                                   # Exit Impeller Sound Speed                                  [m/s]
    Error      = 1  
    while Error>0.0001:
            rho_1=rho_1_h
            C_1_m= m_w_air/(rho_1*A1)                                          # Meridional Velocity                                        [m/s] 
            C_1 = C_1_m                                                        # Absolute velocity at the inlet                             [m/s]
            h_1 = h_01-(0.5*C_1**2)                                            # Static Enthalpy at the inlet                               [J/kg]
            T_1 = h_1/c_p                                                      # Temperature at the Inlet                                   [K]
            p_1   = p_01-0.5*rho_1*C_1**2                                      # Pressure at the inlet                                      [Pa]
            rho_1_h  = p_1/(T_1*R_specific)                                    # Density at the hub diameter                                [kg/m^3]
            Error  = abs(rho_1_h-rho_1)/rho_1
            u_1_t   = Omega*r_1_t                                                  # Peripheral Speed at the tip                                [m/s]
            w_1_t    = sqrt(C_1**2+u_1_t**2)                                       # Relative velocity at the tip                               [m/s]
            a_1_t    = sqrt(k * R_specific * T_01)                                 # Sound at the tip                                           [m/s]
            M_1_t_n= w_1_t / a_1_t                                                 # New Mach Number at the tip                                 [-]
            if M_1_t_n<M_1_t:
              r_1_t_n=(D_1_t/2)+0.001                                            # Radius at the tip
    else:
            break
      # Shaft diameter based on axial load and fatigue
    D1h_03 = 0.2*D_2                                                        # Shaft Diameter Calculation                                    [mm]
    D1h    = int(round(max(D1h_01,D1h_02,D1h_03),0))                        # Final Shaft Diameter Calculation                              [mm]          
    T_2  = T_02-c_2**2/(2*c_p)                                              # Static Exit Temperature                                        [K]
    rho_02 = p_02/(R_specific*T_02)
    p_2=p_02-(1/2)*rho_02*c_2**2 
    rho_2=p_2/(R_specific*T_2)
    p_2=p_02                  
    Error=1
    while Error>0.0001:
        p_2=rho_2*R_specific*T_2
        rho_2n=abs(p_02-p_2)*2/c_2**2                                      
        Error=abs(rho_2n-rho_02)/rho_02
    Q_2      = m_f/rho_2                                                    # Volumetric Flow at the Outlet                                 [m^3/s]
    b2       = Q_2/(f*pi*D_2*cr2)*10**6                                     # Blade Width                                                   [mm]
    rho_n  = rho_01                                                         # Density at the hub                                            [kg/m^3]
    rho_delta=abs(rho_n-rho_02)/rho_02                         
    while rho_delta>0.0001:
        p_2=p_02-0.5*rho_02*w2**2 
        rho_n=R_specific*T_2/p_2
        rho_delta=abs(rho_n-rho_02)/rho_02
    
    #=======================================================================#
    # Diffuser Section
    #=======================================================================#
    # Physical constants                               
    Re = (rho_2 * c_2*D_h)/Mu                                           # Reynolds number at diffuser exit                                 [-]
    Cf = 0.079 / (Re ** 0.25)                                           # Friction factor                                                  [-]
    c_3th= cth2*(r_2 / r_3)+ (2*pi*Cf*rho_2*cth2*(r_3**2-r_2*r_3))/m_f  # Tangential velocity at the diffuser inlet                        [m/s]   
    rho_4   = rho_4new                                                                                             # Exit diffuser density                                      [kg/m^3]
    A_4     = 2*math.pi*r_4*b_4                                                                                    # Exit diffuser area                                         [m^2]
    c_4r    = m_f/(rho_4*A_4)                                                                                      # Exit diffuser radial component of absolute velocity        [m/s]                                       
    c_4th   = (-(r_4/r_3)+math.sqrt((r_4/r_3)**2-4*((2*math.pi*Cf*rho_4*(r_3**2-r_3*r_4))/m_f)*c_3th))/(2*c_3th)   # Exit diffuser tangential component of absolute velocity    [m/s]
    c_4     = sqrt(c_4th**2+c_4r**2)                                                                               # Exit diffuser absolute velocity                            [m/s]                                                                              
    alpha4  = degrees(atan(c_4r/c_4th))                                                                            # Exit diffuser flow angle                                   [-]
    h_02    = h_01+delta_h_0                                                                                       # Exit impeller total enthalpy
    h_04    = h_02                                                                                                 # Exit diffuser total enthalpy
    h_4     = h_04-0.5*c_4**2                                                                                      # Exit diffuser enthalpy
    T_04    = h_04/c_p             # Exit diffuser total temperature
    T_4     = T_04-c_4**2/(2*c_p)   # Exit diffuser temperature                                  
    p_4     = (p_02-p_2)*c_p+p_2   # Exit diffuser pressure
    P_04    = cp*R_specific*T_04   # Exit diffuser stagnation pressure
        
        
        
    rho4_new   = p_4/(R*T_4)                   # New value of exit diffuser density                                        
    error_d    = abs(rho4_new-rho_4)/rho_4   # Exit diffuser density error

    if r_4/r_2 < 0.95*c_2/c_4:
        r_4new   = r_4+Deltar   # New value of exit diffuser raidus
    elif r_4/r_2 > 1.05*c_2/c_4:
        r_4new   = r_4+Deltar   # New value of exit diffuser radius
    else:
        break
   

    #=======================================================================#
    # Impeller Losses
    #=======================================================================#
    #Impeller Incidenece Loss
    D1=sqrt((D_1_t**2+D_1_h**2)/2)                                                                                         # Diameter of the impeller inlet                                                                               # Geometric average diameter
    u1=pi*D1*N/60                                                                                                          # Peripheral Speed at the Inlet of the Impeller
    w1=sqrt(C_1**2+u1**2)                                                                                                  # Relative Velocity at the Inlet of the Impeller                            [-]
    beta1=degrees(atan(u1/C_1_m))                                                                                          # Relative Flow Angle                                                       [-]
    beta1_bl=beta1+degrees(5)                                                                                              # Impeller Inlet Blade Angle                                                         [-]
    finc=0.5                                                                                                               # Incidence Friction Loss Coefficient                                       [-]
    delta_hin=finc*w_1_t**2*sin(beta1-beta1_bl)**2                                                                         # Impeller Incidence Loss                                                   [-]
    
    # Impeller Blade Loading 
    Zspl=0.5*Z                                                                                                             # Number of Splitter Blades                                                 [-]
    Neff=Z+0.75*Zspl                                                                                                       # Number of Effective Blades                                                [-]
    Df=1-(w2/w1)+((delta_h_0/u_2**2)/((w_1_t/w2)*(Neff/pi)*(1-(D_1_t/D_2))+(2*(D_1_t/D_2))))                               # Diffuser Loss Coefficient                                                 [-]
    delta_hbl=0.05*Df**2*u_2**2                                                                                            # Impeller Blade Loading
    
    #Impeller Clearance
    delta_r=0.003*D_2                                                                                                      # Clearance Thickness                                                      [mm]                                                                  [mm]
    delta_hcl=0.6*(delta_r/b2)*cth2*sqrt((4*pi*(r_1_t**2-r_1_h**2)*cth2*C_1_m)/(b2*Neff*(r_2-r_1_t)*(1+rho_2/rho_1)))      # Impeller clearance                                                       [-]
    
    #Impeller Mixing
    epsilon=0.2                                                                                                            # Mixing loss coefficient                                                  [-]                                     
    delta_hmx=(0.5*c_2**2/1+tan(alpha2)**2*(epsilon/1-epsilon)**2)                                                          # Impeller Mixing Loss                                                     [-]
    
    
    #Impeller Friction
    Lax=D_2(0.014+0.023*D_2/D_1_h+2.012*Q_1/u_2*D_2**2)                                                                    # Axial Impeller Length                                                    [mm]
    beta_1_t=atan(u_1_t/C_1_m)                                                                                             # Relative Flow angle at the tip                                           [-]
    beta_1_h=atan(u_1_h/C_1_m)                                                                                             # Relative Flow angle at the hub                                           [-]
    beta2_bl=degrees(90)-(degrees(10)+0.5*N)                                                                             # Blade Metal Angle at the impeller outlet                                 [-]
    Lb=pi/8*(2*r_2-(r_1_t-r_1_h)-b2+2*Lax)*(2/(cos(beta_1_t)+cos(beta_1_t)))/2+cos(beta2_bl)                               # Length of Blade                                                          [mm]  
    Dh_1=pi*(D_1_t**2-D_1_h**2)/(2*pi*D1)+Neff*(D_1_t-D_1_h)                                                               # Diameter of the hub                                                      [mm]
    c_th1=sqrt(C_1**2-C_1_m**2)                                                                                            # Velocity in Tangential Direction                                         [m/s]
    w_1_h = sqrt(C_1**2+u_1_h**2)                                                                                      # Relative Velocity at the Hub                                             [m/s]
    w_bar=c_th1+w_1_t+c_2+2*w_1_h+3*w2/8                                                                                    # Normalized Work Parameter                                                [-]
    deltah_sf=2*Cf(w_bar**2)*Lb/Dh_1                                                                                       # Impeller Skin Friction Coefficient                                       [-]
    
    #Impeller Recirculation
    deltah_rc=8*10**(-5)*sinh(3.5*alpha2**3)*(Df*u_2)**2                                                                   # Impeller Recirculation Coefficient                                       [-]
    
    #Impeller Disc Friction
    rho_bar=(rho_1+rho_2)/2                                                                                                # Average Density                                                          [-]
    deltah_disk=Cf*(rho_bar*r_2**2*u_2**3/4*m_f)                                                                           # Impeller Disk Friction Coefficient                                       [-]
    
    #Impeller Leakage
    Dimp_avg=D1+D_2/2                                                                                                      # Diameter of the average of the Impeller                                  [-]
    b=m_f/rho_1*pi*Dimp_avg*C_1_m                                                                                          # Average Blade Width                                                      [mm]
    b1=2*b-b2                                                                                                              # Blade width of Inlet Impeller                                            [mm]
    b_bar=b1+b2/2                                                                                                          # Average of the blade width                                               [mm]
    r_bar=(r_2+r_1_t)/2                                                                                                    # Average radii of the blade                                               [mm]
    delta_PL=m_f*((r_2*cth2)-(r_1_t*u_1_t))/Neff*r_bar*b_bar*Lax                                                           # TotalPressure Loss                                                       [Pa.s]
    u_l=0.816*sqrt(2*delta_PL/rho_2)                                                                                       # Leakage Velocity                                                         [m/s]
    delta_rad=0.001*D_2                                                                                                    # Normalized Radius                                                        [mm]
    m_fl=rho_2*u_l*delta_rad*Lax*Neff                                                                                      # Mass flow leakage                                                        [kg/s]
    deltah_lk=m_fl*u_l*u_2/2*m_f                                                                                           # Impeller Leakage                                                         [-]
    
    #=======================================================================#
    # Diffuser Losses
    #=======================================================================#
    #Diffuser Friction
    deltah_df=Cf*r_2(1-(r_2/r_3)**1.5)*c_2**2/1.5*b2*cos(alpha2)                                                         # Diffuser Friction                                                         [-]
                                            
    #=======================================================================#
    # Real Efficiency with Losses
    #=======================================================================#                                    
    delta_h_0_id= c_p * T_01 * ((beta_tt**((k - 1) / k)) - 1) / 1
    tot_loss = delta_hin+delta_hbl+delta_hcl+deltah_df+delta_hmx+deltah_sf+deltah_rc+deltah_lk+deltah_disk
    delta_h_real=delta_h_0_id+tot_loss
    Eta_is_n=delta_h_0_id/delta_h_real
    Error=abs(Eta_is_n-Eta_is)/Eta_is
    Eta_is=Eta_is_n

  #=======================================================================#
  # Return Channel Sizing
  #=======================================================================#  
BM6=0.12                                                               # Blockage factour of the return bend
alpha_6_specified =1                 
alpha6=math.degrees.atan(0.32+(1.7*phi_2))                
b_6=b_4*tan(alpha4)/tan(alpha6/(1-BM6))                                 # Blade height at return bend shroud contour exit
if b_6 < b_4:
    b_6 = b_4  # Minimum constraint
    # Recalculate alpha_6 based on the constrained b_6
    alpha_6 = math.degrees(math.atan(b_4 * math.tan(math.radians(alpha4)) / (b_6 * (1 - BM6))))
elif b_6 > 2 * b_4:
    b_6 = 2 * b_4  # Maximum constraint
    # Recalculate alpha_6 based on the constrained b_6
    alpha_6 = math.degrees(math.atan(b_4 * math.tan(math.radians(alpha4)) / (b_6 * (1 - BM6))))
else:
    alpha_6 = alpha_6_specified
R_ch=b_6+b_4/2                                                            # Radius of the return channel
b8=(R_ch/0.8)+b_6                                                         # Passage width of the eye of the next impeller
Aco=R_ch+(b_4+b_6)/2                                                      # Elliptical contour of the axial semi axes
Bco=R_ch+b_4                                                              # Elliptical contour of the radial semi axes
Rc_exs=b8                                                                 # Radius of the curvature of the shroud contour
Rc_exh=2*b8                                                               # Radius of the curvature of the hub contour

   
  #=======================================================================#
  # Return Channel Sizing Losses
  #=======================================================================# 
  # Incidence Losses
alpha_4_opt=degrees(88)
delta_h_0_inc=c_4*sin(abs(alpha4-alpha_4_opt))**2/2


######################SECOND STAGE#######################################
#=======================================================================#
# Impeller Section
#=======================================================================#
Error_eta=1
while Error_eta>0.001:
    Eta_is=Eta_is_n
    delta_h_0 = c_p * T_04 * ((beta_tt**((k - 1) / k)) - 1) / Eta_is        # Enthalpy Change                                            [J/kg]
    delta_T_0 = delta_h_0/c_p                                               # Temperature Change                                         [K]
    h_04s=h_01+Eta_is_n*(h_04-h_01)
    T_04s=h_04s/c_p
    p_04s=p_01*(T_04s/T_01)**(k/k-1)
    p_04=p_04s-(deltah_df+delta_h_0_inc)
    rho_04=P_04/R_specific*T_04
    Q_1     = m_f/rho_01                                                    # Volumetric Flow at the inlet                               [m^3/s]
    P_t        = m_f*delta_h_0                                              # Total Power                                                [KW]
    T_05      = T_04+delta_T_0                                              # Exit Temperature                                           [K]
    beta_2_inf= 50                                                          # Ideal blade flow angle                                     [degree]
    Z      = (int(90-beta_2_inf)/3)                                         # Number of Blades                                           [-]
    phi_2     = 0.22                                                        # Exit flow coefficient                                      [-]
    tau_inf = 1 - phi_2*tan(beta_2_inf)                                     # Ideal Work Coefficient                                     [-]                      
    slip_f = 1 - (cos(radians(beta_2_inf)) ** 0.5) / (Z ** 0.7)             # Slip factor                                                [-]
    tau       = slip_f-phi_2*tan(beta_2_inf)                                # Work Coefficient                                           [-]
    beta2    = degrees(atan(c_4r/c_4th))                                      # Blade Metal Angle                                          [deg]
    alpha2   = degrees(atan(c_4r/c_4))                                        # Absolute Flow Angle                                       [deg]
    t          = 4                                                          # Blade thickness                                            [mm]
    while D_1_t > D_1_h:
        rho_1_h    = p_04/T_04*R_specific                                   # Density at the Hub Diameter                               [kg/m^3]
        Error      = 1  
        while Error>0.0001:
            rho_1=rho_1_h
            C_1_m= m_w_air/(rho_1*A1)                                          # Meridional Velocity                                        [m/s] 
            C_1 = C_1_m                                                        # Absolute velocity at the inlet                             [m/s]
            h_4 = h_04-(0.5*C_1**2)                                            # Static Enthalpy at the inlet                               [J/kg]
            T_4 = h_4/c_p                                                      # Temperature at the Inlet                                   [K]
            p_4   = p_04-0.5*rho_1*C_1**2                                      # Pressure at the inlet                                      [Pa]
            rho_1_h  = p_4/(T_4*R_specific)                                    # Density at the hub diameter                                [kg/m^3]
            Error  = abs(rho_1_h-rho_1)/rho_1
        a_1_t    = sqrt(k * R_specific * T_04)                                 # Sound at the tip                                           [m/s]
        M_1_t_n= w_1_t / a_1_t                                                 # New Mach Number at the tip                                 [-]
        if M_1_t_n<M_1_t:
            r_1_t_n=(D_1_t/2)+0.001                                            # Radius at the tip
        else:
            break         
    T_5  = T_05-c_2**2/(2*c_p)                                              # Static Exit Temperature                                        [K]
    p_05=p_04*beta_tt
    rho_05 = p_05/(R_specific*T_05)
    p_5=p_05-(1/2)*rho_05*c_2**2 
    rho_5=p_5/(R_specific*T_5)
    p_5=p_05                  
    Error=1
    while Error>0.0001:
        p_5=rho_5*R_specific*T_5
        rho_5n=abs(p_05-p_5)*2/c_2**2                                      
        Error=abs(rho_5n-rho_05)/rho_05
    Q_2      = m_f/rho_5                                                    # Volumetric Flow at the Outlet                                 [m^3/s]
    b2       = Q_2/(f*pi*D_2*cr2)*10**6                                     # Blade Width                                                   [mm]
    L_x      = D_2-(D_1_t-D_1_h)/2/(2+b2/2)                                 # Impeller Axial Width                                          [mm]
    rho_n  = rho_04                                                         # Density at the hub                                            [kg/m^3]
    rho_delta=abs(rho_n-rho_05)/rho_05                         
    while rho_delta>0.0001:
        p_5=p_05-0.5*rho_05*w2**2 
        rho_n=R_specific*T_5/p_5
        rho_delta=abs(rho_n-rho_05)/rho_05

    #=======================================================================#
    # Diffuser Section
    #=======================================================================#
    Re = (rho_5 * c_2*D_h)/Mu                                           # Reynolds number at diffuser exit                                 [-]
    Cf = 0.079 / (Re ** 0.25)                                           # Friction factor                                                  [-]
    c_3th= cth2*(r_2 / r_3)+ (2*pi*Cf*rho_2*cth2*(r_3**2-r_2*r_3))/m_f  # Tangential velocity at the diffuser inlet                        [m/s]   


    rho_6    = rho_5                                                # Inlet diffuser density
    while error_d > 0.0001:
                            
        c_4th   = (-(r_4/r_3)+math.sqrt((r_4/r_3)**2-4*((2*math.pi*Cf*rho_4*(r_3**2-r_3*r_4))/m_f)*c_3th))/(2*c_3th)   # Exit diffuser tangential component of absolute velocity    [m/s]
        c_4     = sqrt(c_4th**2+c_4r**2)                                                                               # Exit diffuser absolute velocity                            [m/s]                                                                              
        alpha4  = degrees(atan(c_4r/c_4th))                                                                            # Exit diffuser flow angle                                   [-]
        h_05   = h_04+delta_h_0                                                                                       # Exit impeller total enthalpy
        h_07   = h_05                                                                                                 # Exit diffuser total enthalpy
        
        
        T_07    = h_07/c_p             # Exit diffuser total temperature
        T_7     = T_07-c_4**2/(2*c_p)   # Exit diffuser temperature                                  
        p_7     = (p_05-p_5)*c_p+p_5   # Exit diffuser pressure      
        
   

    #=======================================================================#
    # Impeller Losses
    #=======================================================================#
    #Impeller Incidenece Loss
    w1=sqrt(C_1**2+u1**2)                                                                                                  # Relative Velocity at the Inlet of the Impeller                            [-]
    beta1=degrees(atan(u1/C_1_m))                                                                                          # Relative Flow Angle                                                       [-]
    beta1_bl=beta1+degrees(5)                                                                                              # Impeller Inlet Blade Angle                                                         [-]
    delta_hin=finc*w_1_t**2*sin(beta1-beta1_bl)**2                                                                         # Impeller Incidence Loss                                                   [-]
    
    # Impeller Blade Loading 
    Df=1-(w2/w1)+((delta_h_0/u_2**2)/((w_1_t/w2)*(Neff/pi)*(1-(D_1_t/D_2))+(2*(D_1_t/D_2))))                               # Diffuser Loss Coefficient                                                 [-]
    delta_hbl=0.05*Df**2*u_2**2                                                                                            # Impeller Blade Loading
    
    #Impeller Clearance
    delta_hcl=0.6*(delta_r/b2)*cth2*sqrt((4*pi*(r_1_t**2-r_1_h**2)*cth2*C_1_m)/(b2*Neff*(r_2-r_1_t)*(1+rho_5/rho_4)))      # Impeller clearance                                                       [-]
    
    #Impeller Mixing              
    delta_hmx=(0.5*c_2**2/1+tan(alpha2)**2*(epsilon/1-epsilon)**2)                                                          # Impeller Mixing Loss                                                     [-]
    
    
    #Impeller Friction
    Lax=D_2(0.014+0.023*D_2/D_1_h+2.012*Q_1/u_2*D_2**2)                                                                    # Axial Impeller Length                                                    [mm]
    beta_1_t=atan(u_1_t/C_1_m)                                                                                             # Relative Flow angle at the tip                                           [-]
    beta_1_h=atan(u_1_h/C_1_m)                                                                                             # Relative Flow angle at the hub                                           [-]
    beta2_bl=degrees(90)-(degrees(10)+0.5*N)                                                                             # Blade Metal Angle at the impeller outlet                                 [-]
    Lb=pi/8*(2*r_2-(r_1_t-r_1_h)-b2+2*Lax)*(2/(cos(beta_1_t)+cos(beta_1_t)))/2+cos(beta2_bl)                               # Length of Blade                                                          [mm]  
    c_th1=sqrt(C_1**2-C_1_m**2)                                                                                            # Velocity in Tangential Direction                                         [m/s]
    w_1_h = sqrt(C_1**2+u_1_h**2)                                                                                      # Relative Velocity at the Hub                                             [m/s]
    w_bar=c_th1+w_1_t+c_2+2*w_1_h+3*w2/8                                                                                    # Normalized Work Parameter                                                [-]
    deltah_sf=2*Cf(w_bar**2)*Lb/Dh_1                                                                                       # Impeller Skin Friction Coefficient                                       [-]
    
    #Impeller Recirculation
    deltah_rc=8*10**(-5)*sinh(3.5*alpha2**3)*(Df*u_2)**2                                                                   # Impeller Recirculation Coefficient                                       [-]
    
    #Impeller Disc Friction
    rho_bar=(rho_4+rho_5)/2                                                                                                # Average Density                                                          [-]
    deltah_disk=Cf*(rho_bar*r_2**2*u_2**3/4*m_f)                                                                           # Impeller Disk Friction Coefficient                                       [-]
    
    #Impeller Leakage
    delta_PL=m_f*((r_2*cth2)-(r_1_t*u_1_t))/Neff*r_bar*b_bar*Lax                                                           # TotalPressure Loss                                                       [Pa.s]
    deltah_lk=m_fl*u_l*u_2/2*m_f                                                                                           # Impeller Leakage                                                         [-]
    
    #=======================================================================#
    # Diffuser Losses
    #=======================================================================#
    #Diffuser Friction
    deltah_df=Cf*r_2(1-(r_2/r_3)**1.5)*c_2**2/1.5*b2*cos(alpha2)                                                          # Diffuser Friction                                                         [-]
                                            
    #=======================================================================#
    # Real Efficiency with Losses
    #=======================================================================#                                    
    delta_h_0_id= c_p * T_04 * ((beta_tt**((k - 1) / k)) - 1) / 1
    tot_loss = delta_hin+delta_hbl+delta_hcl+deltah_df+delta_hmx+deltah_sf+deltah_rc+deltah_lk+deltah_disk
    delta_h_real=delta_h_0_id+tot_loss
    Eta_is_n=delta_h_0_id/delta_h_real
    Error=abs(Eta_is_n-Eta_is)/Eta_is
    Eta_is=Eta_is_n

#=======================================================================#
# Volute Design for Second Stage (After Diffuser)
#=======================================================================#
r_3_volute = r_4                           # base circumference
c_th2_volute = c_4th                       # tangential velocity entering the volute
Q_volute = m_f /  rho_6                      # volumetric flow rate entering the volute
alpha_2_volute = alpha4                     # flow angle entering the volute

# Volute parameters
t_v = 0.003  # Tongue thickness (min 3 mm)
k_b = 1.03   

# Base circumference
r_3 = k_b * r_2  

# Initialize volute table
volute_table = []
angles = [math.pi/4, math.pi/2, 3*math.pi/4, math.pi,
          5*math.pi/4, 3*math.pi/2, 7*math.pi/4, 2*math.pi]

# Calculate volute parameters for each azimuthal angle
for theta_i in angles:
    # Solve for A_t using conservation of angular momentum
    # Equation: A_t * (r_2 * c_th2_volute / Q_volute) - sqrt(A_t/π) - (r_3 + t_v) = 0
   
    # Define the function to solve
    def volute_area_eq(A_t):
        return A_t * (r_2 * c_th2_volute / Q_volute) - math.sqrt(A_t/math.pi) - (r_3 + t_v)
   
    # Solve using Newton-Raphson method
    A_t_guess = Q_volute * theta_i / (2 * math.pi * c_th2_volute)
    tolerance = 1e-6
    max_iter = 100
   
    for iter in range(max_iter):
        f_val = volute_area_eq(A_t_guess)
        f_deriv = (r_2 * c_th2_volute / Q_volute) - (1/(2*math.sqrt(math.pi*A_t_guess)))
       
        # Avoid division by zero
        if abs(f_deriv) < 1e-10:
            break
           
        A_t_new = A_t_guess - f_val / f_deriv
       
        if abs(A_t_new - A_t_guess) < tolerance:
            A_t_guess = A_t_new
            break
           
        A_t_guess = A_t_new
   
    # Calculate section radius for circular cross-section
    r_si = math.sqrt(A_t_guess / math.pi)
   
    # Calculate local external radius
    r_i = r_3 * math.exp(theta_i * math.tan(math.radians(alpha_2_volute)))
   
    # Calculate velocity at this section
    c_i = Q_volute / A_t_guess
   
    # Store results
    volute_table.append({
        'theta_i': math.degrees(theta_i),
        'r_i': r_i,
        'A_i': A_t_guess,
        'r_si': r_si,
        'c_i': c_i
    })

# Print volute design table
print("\nVolute Design Parameters for Second Stage:")
print("θ_i [deg] | r_i [m] | A_i [m²] | r_si [m] | c_i [m/s]")
print("-" * 55)
for row in volute_table:
    print(f"{row['theta_i']:8.1f} | {row['r_i']:7.4f} | {row['A_i']:8.6f} | {row['r_si']:7.4f} | {row['c_i']:7.2f}")
# =======================================================================#
# Calculate Volute Losses
# =======================================================================#
delta_mv=rho_6*(c_4r**2)/2                                     #meridional velocity dump losses at the impeller exit                                                                                                                                                                                                                                                                                                        
delta_h_volute =  delta_mv/ rho_6
