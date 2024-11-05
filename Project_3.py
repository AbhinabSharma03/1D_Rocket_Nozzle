#Quasi One-dimensional Fluid Flow through a Converging-diverging rocket Nozzle
#Abhinab Sharma

import numpy as np
import matplotlib.pyplot as plt

isentropic_subsonic = False # True for the subsonic condition and False for the Subsonic-Supersonic Condition

if isentropic_subsonic:

    # Constants and initialization
    dx = 0.01  
    c = 0.64  #courant number
    gamma = 1.4 
    R = 287.0

    # Discretization of nondimensional distance along the nozzle
    x = np.arange(0, 3, dx)
    k = len(x)
    n = 25  # Number of vertical points

    # Area distribution along the nozzle
    A = np.zeros(k)

    # Parabolic area distribution along the nozzle
    for i in range(k):
        if x[i] < 1.5:
            A[i] = 1 + 2.2 * (x[i] - 1.5) ** 2
        else:
            A[i] = 1 + 0.2 * (x[i] - 1.5) ** 2

    # Grid distribution in the nozzle
    xx = np.zeros((k, n))
    yy = np.zeros((k, n))

    for i in range(k):
        y = np.linspace(-A[i] / 2, A[i] / 2, n)  # Y distribution
        xx[i, :] = x[i]  # x positions
        yy[i, :] = y  ## y positions

    ## Initialization of variables
    rho = np.zeros(k)
    V = np.zeros(k)
    T = np.zeros(k)
    
    ## Time Derivatives of the variables
    d_rho_dt = np.zeros(k)
    d_V_dt = np.zeros(k)
    d_T_dt = np.zeros(k)
    dt_arr = np.zeros(k)  # Time step sizes

    # Initial Conditions for Subsonic Flow
    P0 = 1.0  # inlet pressure
    P_out_ratio = 0.928 # Outlet pressure ratio
    P = np.zeros(k) 

    # Set initial conditions
    P[0] = P0  # Inlet pressure
    P[-1] = P0 * P_out_ratio  # Outlet pressure
    T[0] = 1.0
    rho[0] = P[0] / (R * T[0])  # Density using ideal gas law

    # Initialize density, temperature, and velocity using linear gradients
    for i in range(1, k):
        if i < k - 1:
            # Linear interpolation for pressure
            P[i] = P[0] - (P[0] - P[-1]) * (i / (k - 1))
            T[i] = T[0] * (P[i] / P[0])  # Adjust temperature based on pressure
            rho[i] = P[i] / (R * T[i])  # Density from ideal gas law
        else:
            P[i] = P[-1]
            T[i] = T[i - 1] * (P[i] / P[i - 1])  # Final temperature adjustment
            rho[i] = P[i] / (R * T[i])  # Final density adjustment

    
    V = 0.05 + 0.11 * x  # Initial velocity distribution

    
    a = np.sqrt(gamma * R * T)  # Speed of sound
    a[a < 1e-3] = 1e-3  # Preventing speed of sound from being zero

   
    M = V / a  # Mach number

    
    mf = rho * A * V  # Mass flow rate

    # Prepare for plotting results
    max_steps = 1000  # Max time steps
    res_max = 1e-8  # Residual threshold
    res = 1  # Initial residual
    t = 0  # Current time
    step_count = 0  # Time step counter

    # Time-stepping loop
    rho_vals = np.zeros((max_steps, k))  # Density values
    T_vals = np.zeros((max_steps, k))  # Temperature values
    M_vals = np.zeros((max_steps, k))  # Mach number values
    P_vals = np.zeros((max_steps, k))  # Pressure values

    while res > res_max and step_count < max_steps:
        prev_rho = rho.copy()  # Previous density
        prev_V = V.copy()  # Previous velocity
        prev_T = T.copy()  # Previous temperature

        # Calculate time step sizes
        for i in range(k):
            dt_arr[i] = (c * dx) / (T[i] ** 0.5 + np.abs(V[i]))
        dt = np.min(dt_arr)  # Minimum time step size

        # Predictor Step
        d_rho_pred = np.zeros(k) 
        d_V_pred = np.zeros(k)  
        d_T_pred = np.zeros(k)  

        for i in range(1, k - 1):
            dV = (V[i + 1] - V[i]) / dx  # Velocity gradient
            dA = (np.log(A[i + 1]) - np.log(A[i])) / dx  # Area gradient
            drho = (rho[i + 1] - rho[i]) / dx  # Density gradient
            dT = (T[i + 1] - T[i]) / dx  # Temperature gradient
            d_rho_pred[i] = (-rho[i] * dV - rho[i] * V[i] * dA - V[i] * drho)   # Density equation
            d_V_pred[i] = (-V[i] * dV - (1 / gamma) * dT - (1 / gamma) * drho * T[i] / rho[i])  # Velocity equation
            d_T_pred[i] = (-V[i] * dT - (gamma - 1) * T[i] * dV - (gamma - 1) * T[i] * V[i] * dA)  # Temperature equation

        # Predicted values
        rho_p = rho.copy() 
        V_p = V.copy() 
        T_p = T.copy()  

        for i in range(1, k - 1):
            rho_p[i] += d_rho_pred[i] * dt  # Update density prediction
            V_p[i] += d_V_pred[i] * dt  # Update velocity prediction
            T_p[i] += d_T_pred[i] * dt  # Update temperature prediction

        t += dt  # Update current time
        step_count += 1  # Increment time step counter

        # Boundary conditions at first node
        rho_p[0] = rho[0]  
        V_p[0] = V[0] 
        T_p[0] = T[0] 

        # Corrector Step
        d_rho_corr = np.zeros(k) 
        d_V_corr = np.zeros(k) 
        d_T_corr = np.zeros(k)  

        for i in range(1, k - 1):
            dV = (V_p[i] - V_p[i - 1]) / dx  # Velocity gradient for correction
            dA = (np.log(A[i]) - np.log(A[i - 1])) / dx  # Area gradient for correction
            drho = (rho_p[i] - rho_p[i - 1]) / dx  # Density gradient for correction
            dT = (T_p[i] - T_p[i - 1]) / dx  # Temperature gradient for correction
            d_rho_corr[i] = (-rho_p[i] * dV - rho_p[i] * V_p[i] * dA - V_p[i] * drho)  # Corrected density equation
            d_V_corr[i] = (-V_p[i] * dV - (1 / gamma) * dT - (1 / gamma) * drho * T_p[i] / rho_p[i])  # Corrected velocity equation
            d_T_corr[i] = (-V_p[i] * dT - (gamma - 1) * T_p[i] * dV - (gamma - 1) * T_p[i] * V_p[i] * dA)  # Corrected temperature equation

        # Average time derivatives
        for i in range(1, k - 1):
            d_rho_dt[i] = 0.5 * (d_rho_pred[i] + d_rho_corr[i])  
            d_V_dt[i] = 0.5 * (d_V_pred[i] + d_V_corr[i])  
            d_T_dt[i] = 0.5 * (d_T_pred[i] + d_T_corr[i])  

        # Update the variables
        for i in range(1, k - 1):
            rho[i] += d_rho_dt[i] * dt 
            V[i] += d_V_dt[i] * dt  
            T[i] += d_T_dt[i] * dt  

    
        T[T < 1e-3] = 1e-3  # Preventing temperature from going below threshold

        # Recalculate speed of sound
        a = np.sqrt(gamma * R * T)
        a[a < 1e-3] = 1e-3  # Preventing speed of sound from being zero

        # Calculating Mach number again, ensuring no division by zero
        M = np.where(a != 0, V / a, 0)

        # Store results for plotting
        if step_count < max_steps:
            rho_vals[step_count, :] = rho  #density
            T_vals[step_count, :] = T  #temperature
            M_vals[step_count, :] = M  #Mach number
            P_vals[step_count, :] = rho * T  #pressure

        # Residual calculations
        res_rho = np.max(np.abs(rho - prev_rho))  # for density
        res_V = np.max(np.abs(V - prev_V))  # for velocity
        res_T = np.max(np.abs(T - prev_T))  # for temperature
        res = max(res_rho, res_V, res_T)  # Overall residual

    # Final values after simulation
    rho = rho_vals[step_count - 1, :]  #density
    T = T_vals[step_count - 1, :]  #temperature
    V = V  #velocity
    a = np.sqrt(gamma * R * T)  #speed of sound
    M = np.where(a != 0, V / a, 0)  #Mach number
    P = rho * T  #pressure

    # Normalize pressure for plotting
    P_norm = P / P[0]

    # Plotting Pressure, Density, Temperature, and Mach Number for the steady state
    plt.figure(figsize=(12, 10))

    # Pressure Plot
    plt.subplot(4, 1, 1)
    plt.plot(x, P_norm, label='Steady State', color='blue')
    plt.title('Pressure Distribution along the Nozzle (Subsonic)')
    plt.xlabel('Length of the Nozzle (x)')
    plt.ylabel('Pressure (P/P0)')
    plt.grid()

    # Density Plot
    plt.subplot(4, 1, 2)
    plt.plot(x, rho / rho[0], label='Steady State', color='green')
    plt.title('Density Distribution along the Nozzle (Subsonic)')
    plt.xlabel('Length of the Nozzle (x)')
    plt.ylabel('Density (rho/rho0)')
    plt.grid()

    # Temperature Plot
    plt.subplot(4, 1, 3)
    plt.plot(x, T, label='Steady State', color='red')
    plt.title('Temperature Distribution along the Nozzle (Subsonic)')
    plt.xlabel('Length of the Nozzle (x)')
    plt.ylabel('Temperature (T/T0)')
    plt.grid()

    # Mach Number Plot
    plt.subplot(4, 1, 4)
    plt.plot(x, M, label='Steady State', color='purple')
    plt.title('Mach Number Distribution along the Nozzle (Subsonic)')
    plt.xlabel('Length of the Nozzle (x)')
    plt.ylabel('Mach Number (M)')
    plt.grid()

    plt.tight_layout()
    plt.show()

    # Print final values
    print("Final Pressure (P/P0):", P / P[0])
    print("Final Density (rho/rho0):", rho / rho[0])
    print("Final Temperature (T/T0):", T)
    print("Final Velocity (V/a0):", V / V[0])
    print("Final Mach Number (M):", M)


else: 

    # Constants and initialization
    dx = 0.01  
    c = 0.64  # Courant number
    gamma = 1.4 

    # Discretization of nondimensional distance along the nozzle
    x = np.arange(0, 3, dx) 
    k = len(x) 
    n = 25 

    # Parabolic area distribution along the nozzle
    A = 1 + 2.2 * (x - 1.5) ** 2 

    # Grid distribution in the nozzle
    xx = np.zeros((k, n)) 
    yy = np.zeros((k, n))  

    for i in range(k):
        y = np.linspace(-A[i] / 2, A[i] / 2, n)  # Y distribution
        xx[i, :] = x[i]  # Assign x positions
        yy[i, :] = y  # Assign y positions

    # Initialization of variables
    d_rho_dt = np.zeros(k)  # Density time derivative
    d_V_dt = np.zeros(k)  # Velocity time derivative
    d_T_dt = np.zeros(k)  # Temperature time derivative
    dt_arr = np.zeros(k)  # Time step sizes

    
    # Initial Conditions
    # Supersonic initial conditions
    rho = 1 - 0.3146 * x  
    T = 1 - 0.2314 * x  
    V = (0.1 + 1.09 * x) * T ** 0.5

    P = rho * T  # Pressure (P/P0)
    M = V * T ** -0.5  # Mach number
    mf = rho * A * V  # Mass flow rate

    # Mass flow rates at different time steps
    max_steps = 1000  # Max time steps
    rho_vals = np.zeros((max_steps, k))  # Density values
    T_vals = np.zeros((max_steps, k))  # Temperature values
    M_vals = np.zeros((max_steps, k))  # Mach number values
    P_vals = np.zeros((max_steps, k))  # Pressure values

    # Prepare for plotting results
    res_max = 1e-8 
    res = 1 
    t = 0 
    step_count = 0 

    # Time-stepping loop
    while res > res_max and step_count < max_steps:
        prev_rho = d_rho_dt.copy()
        prev_V = d_V_dt.copy() 
        prev_T = d_T_dt.copy() 

        # Calculate time step sizes
        for i in range(k):
            dt_arr[i] = (c * dx) / (T[i] ** 0.5 + V[i])
        dt = np.min(dt_arr)  # Minimum time step size

        # Predictor Step
        d_rho_pred = np.zeros(k)
        d_V_pred = np.zeros(k)
        d_T_pred = np.zeros(k)

        #gradients of variables  
        for i in range(1, k - 1):
            dV = (V[i + 1] - V[i]) / dx 
            dA = (np.log(A[i + 1]) - np.log(A[i])) / dx  
            drho = (rho[i + 1] - rho[i]) / dx 
            dT = (T[i + 1] - T[i]) / dx  
            d_rho_pred[i] = (-rho[i] * dV - rho[i] * V[i] * dA - V[i] * drho)  
            d_V_pred[i] = (-V[i] * dV - (1 / gamma) * dT - (1 / gamma) * drho * T[i] / rho[i])  
            d_T_pred[i] = (-V[i] * dT - (gamma - 1) * T[i] * dV - (gamma - 1) * T[i] * V[i] * dA) 

        # Predicted values
        rho_p = rho.copy()
        V_p = V.copy()  
        T_p = T.copy()

        #update variables
        for i in range(1, k - 1):
            rho_p[i] += d_rho_pred[i] * dt  
            V_p[i] += d_V_pred[i] * dt  
            T_p[i] += d_T_pred[i] * dt  

        t += dt  # Update current time
        step_count += 1  # Increment time step counter

        # Boundary conditions at first node
        rho_p[0] = rho[0] 
        V_p[0] = V[0] 
        T_p[0] = T[0]  

        # Corrector Step
        d_rho_corr = np.zeros(k)
        d_V_corr = np.zeros(k) 
        d_T_corr = np.zeros(k)

        #gradients of variables 

        for i in range(1, k - 1):
            dV = (V_p[i] - V_p[i - 1]) / dx  
            dA = (np.log(A[i]) - np.log(A[i - 1])) / dx 
            drho = (rho_p[i] - rho_p[i - 1]) / dx 
            dT = (T_p[i] - T_p[i - 1]) / dx 
            d_rho_corr[i] = (-rho_p[i] * dV - rho_p[i] * V_p[i] * dA - V_p[i] * drho) 
            d_V_corr[i] = (-V_p[i] * dV - (1 / gamma) * dT - (1 / gamma) * drho * T_p[i] / rho_p[i])
            d_T_corr[i] = (-V_p[i] * dT - (gamma - 1) * T_p[i] * dV - (gamma - 1) * T_p[i] * V_p[i] * dA)  
        
        # Average time derivatives
        for i in range(1, k - 1):
            d_rho_dt[i] = 0.5 * (d_rho_pred[i] + d_rho_corr[i]) 
            d_V_dt[i] = 0.5 * (d_V_pred[i] + d_V_corr[i])
            d_T_dt[i] = 0.5 * (d_T_pred[i] + d_T_corr[i])  

        # Update the variables
        for i in range(1, k - 1):
            rho[i] += d_rho_dt[i] * dt
            V[i] += d_V_dt[i] * dt
            T[i] += d_T_dt[i] * dt  

        # Calculate pressure and store results for all positions
        P = rho * T  
        if step_count < max_steps:
            rho_vals[step_count, :] = rho  
            T_vals[step_count, :] = T  
            M_vals[step_count, :] = M  
            P_vals[step_count, :] = P  

        # Residual calculations
        res_rho = np.max(np.abs(rho - prev_rho)) 
        res_V = np.max(np.abs(V - prev_V))  
        res_T = np.max(np.abs(T - prev_T))  
        res = max(res_rho, res_V, res_T)  

    # Final values after simulation
    rho = rho_vals[step_count - 1, :]  
    T = T_vals[step_count - 1, :]  
    V = V  
    M = V * T ** -0.5 
    P = rho * T  

    # Plotting Pressure, Density, Temperature, and Mach Number for the steady state
    plt.figure(figsize=(12, 10))

    # Pressure Plot
    plt.subplot(4, 1, 1)
    plt.plot(x, P, label='Steady State', color='blue')
    plt.title('Pressure Distribution along the Nozzle')
    plt.xlabel('Length of the Nozzle (x)')
    plt.ylabel('Pressure (P/P0)')
    plt.grid()

    # Density Plot
    plt.subplot(4, 1, 2)
    plt.plot(x, rho, label='Steady State', color='green')
    plt.title('Density Distribution along the Nozzle')
    plt.xlabel('Length of the Nozzle (x)')
    plt.ylabel('Density (rho/rho0)')
    plt.grid()

    # Temperature Plot
    plt.subplot(4, 1, 3)
    plt.plot(x, T, label='Steady State', color='red')
    plt.title('Temperature Distribution along the Nozzle')
    plt.xlabel('Length of the Nozzle (x)')
    plt.ylabel('Temperature (T/T0)')
    plt.grid()

    # Mach Number Plot
    plt.subplot(4, 1, 4)
    plt.plot(x, M, label='Steady State', color='purple')
    plt.title('Mach Number Distribution along the Nozzle')
    plt.xlabel('Length of the Nozzle (x)')
    plt.ylabel('Mach Number (M)')
    plt.grid()

    plt.tight_layout()
    plt.show()

    # Print final values
    print("Final Pressure (P/P0):", P / P[0])
    print("Final Density (rho/rho0):", rho / rho[0])
    print("Final Temperature (T/T0):", T)
    print("Final Velocity (V/a0):", V / V[0])
    print("Final Mach Number (M):", M)
