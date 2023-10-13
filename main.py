#Rocket Sim 3 - By: Reese Kimmel - 10/12/2023

import math

# Define the stage 1 variables
Stage1_area = 0.0027  # m^2
drag_coefficient = 1.05  # drag coefficient of the rocket
hm = 0  # max height, m
ht = 0  # time of max height, s
vm = 0  # max velocity, m/s
pa = 0  # max acceleration, m/s^2
mt = 0  # max (estimated) temperature of the rocket, F
ad = 0  # density of the atmosphere at apogee, kg/m^3
vmt = 0  # time of max velocity, s
mtt = 0  # time of max temperature, s
mva = 0  # altitude of max velocity, m
pmva = 0  # Density at altitude of max velocity, kg/m^3
tmv = 0  # temperature of max velocity, f
mtpv = 0  # velocity at max temperature, m/s
mtpa = 0  # altitude of max temperature, m
ma = 0  # max acceleration, m/s^2
ter = 0  # time of engine run-out, s
pd = 0  # peak drag, N
loop = False
drag_force = 0  # N

fuel_mass = 0.020  # kg
starting_fuel_mass = fuel_mass  # kg, defined to calculate the decreasing fuel mass
payload_mass = 0.0  # kg, 1190 grams of payload (the second stage)
misc_mass = 0.0  # kg, miscellaneous bits of extra mass, fuel tubes, nozzles, etc...
body_mass = 0.22  # kg, mass of the body of the rocket

payload_mass += misc_mass
payload_mass += body_mass
total_mass = payload_mass  # kg
total_mass += fuel_mass

print("Total mass: " + str(total_mass) + " kg")
Ns = 11.5  # The amount of fuel spent per second, N
impulse = 30  # The total impulse of the engine, N*s
start_impulse = impulse  # Ns, defined to calculate the fuel mass

# Other miscellaneous variables
Parachute_delay = 8.7  # time delay of parachute deployment, s
Parachute_diameter = 10  # diameter of parachute, inches
Parachute_area = math.pi * ((0.5 * (Parachute_diameter / 39.4)) ** 2)  # area of parachute, m^2
Parachute_drag_coefficient = 0.2  # drag coefficient of parachute
Parachute_deployment_velocity = 0  # velocity at which the parachute deploys, m/s
Parachute_jolt_force = 0  # force of the parachute jolt, N
Parachute_deployment_altitude = 0  # altitude at which the parachute deploys, m
Parachute_deployed = False  # whether the parachute has been deployed
Odds_of_parachute_failure = 0  # odds of the parachute failing to deploy, %

t = 0  # time, s
vx = 0  # velocity on the X-Axis, m/s
vz = 0  # velocity on the Z-Axis, m/s
tilt = 0  # tilt of the rocket, degrees
h = 0  # height, m
dx = 0  # distance on the X-Axis, m
dy = 0  # distance on the Y-Axis, m
stage1_loss_percentage = 0
rocket_temp = 0  # temperature of rocket
gravity = 9.81  # m/s^2
dt = 0.01  # time step, s
acceleration = 0  # m/s^2
max_tilt = 0  # max tilt of the rocket, degrees
impact_time = 0  # time of rocket impact with ground
height_offset = 250  # starting altitude above sea level
wind_speed = 1  # wind speed, m/s
sea_level_speed = wind_speed  # wind speed at sea level, m/s
z_percent = 0  # percentage of that the rocket is from vertical
rail_length = 1  # length of the launch rail, m
forced_tilt = 0  # the tilt of the rocket forced by the launch rail, degrees
start_mass = total_mass  # kg

# account for losses
impulse = impulse * (1 - (stage1_loss_percentage / 100))


# define function to calculate aerodynamic heating in degrees F given the velocity in m/s and the density in kg/m^3
def aerodynamic_heating(velocity, altitude):
    altitude = altitude + height_offset
    ambient_temp = 0
    # the aerodynamic heating is modeled as an exponential decay
    velocity = math.sqrt(velocity ** 2)
    velocity = velocity * 3.28084  # convert to ft/s
    heating = 0.00000
    ambient = 72 - (11.7 * (altitude / 1000))  # approximate ambient temperature at altitude, F
    start_density = 0.0752  # approximate density of atmosphere at sea level, kg/m^3
    # the atmospheric thinning is modeled as an exponential decay
    alt_density = start_density * (2.71828 ** (-0.00012 * altitude))
    percentage = alt_density / start_density
    heating = ambient + ((velocity / 8) * percentage)
    if ambient + ((velocity / 8) * percentage) > 0 and heating < ambient + ((velocity / 8) * percentage):
        heating += (ambient + ((velocity / 8) * percentage) / 100)
    if 0 > ambient + ((velocity / 8) * percentage) > heating:
        heating -= (ambient + ((velocity / 8) * percentage) / 100)
    return heating

# define function to calculate the thickness in kg/m^3 of the atmosphere at a given height
def atmosphere(height):
    start_density = 1.225  # approximate density of atmosphere at sea level, kg/m^3
    # the atmospheric thinning is modeled as an exponential decay
    alt_density = start_density * (2.71828 ** (-0.00012 * (height + height_offset)))
    return alt_density

def tilt2(windspeed, velocity):
    if windspeed > 0:
        tilt = (math.atan(windspeed / velocity) * 0.2)
    return tilt

def wind(height, wind_speed):
    # the wind speed is modeled as an exponential decay
    wind_speed = wind_speed * 2.71828 ** (-0.00012 * (height + height_offset))
    return wind_speed

# the main loop is defined, where all the main physics is done
while h > -1:
    density = atmosphere(h)  # calculate the density of the atmosphere at the current height

    wind_speed = wind(h, sea_level_speed)

    if forced_tilt == 0:
        # calculate the tilt of the rocket
        if h > rail_length:
            tilt = tilt2(wind_speed, vz)
        if tilt > max_tilt:
            max_tilt = tilt
        tilt = max_tilt
    else:
        tilt = forced_tilt

    # find the x velocity of the rocket
    if impulse > 0:
        vx += math.sin(math.radians(tilt)) * (Ns / total_mass) * dt

    # calculate the distance the rocket has traveled
    dx += (vx * dt)
    z_percent = math.cos(math.radians(tilt))
    rho = density  # kg/m^3
    fin_authority = 0.01 * rho * vz ** 2  # drag force in N
    print("Rocket drifting parameters:", "Tilt:", "{:.2f}".format(tilt), "Velocity:", "{:.2f}".format(vx),
          "Distance:", str(int(dx)) + "m,", "Fin Authority:", "{:.2f}".format(fin_authority))

    # calculate for the rocket engine spending fuel
    if impulse > 0:  # if the rocket has fuel, subtract the fuel
        impulse = impulse - Ns * dt
        fuel_mass = starting_fuel_mass * (impulse / start_impulse)  # calculate the fuel mass
        total_mass = fuel_mass + payload_mass
        vz = vz + ((Ns / total_mass) * dt) * z_percent  # calculate the velocity
    if impulse < 0:  # fix a bug where the impulse would go under 0
        impulse = 0

    # calculate for the drag force
    rho = density  # kg/m^3
    Cd = drag_coefficient  # drag coefficient
    A = Stage1_area  # M^2

    straightV = math.sqrt(vx**2 + vz**2)

    # calculate the drag force
    drag_force_N = 0.5 * rho * Cd * A * (straightV ** 2)  # drag force in newtons
    if straightV > 0:  # if the rocket is moving up, subtract the drag force
        straightV -= ((drag_force_N / total_mass) * dt)
    if vx < 0:  # if the rocket is moving down, add the drag force
        straightV += ((drag_force_N / total_mass) * dt)

    if vz < 0:
        vz = -1 * straightV * math.cos(math.radians(tilt))  # calculate the vertical velocity
    else:
        vz = straightV * math.cos(math.radians(tilt))
    if vx < 0:
        vx = -1 * straightV * math.sin(math.radians(tilt))  # calculate the vertical velocity
    else:
        vx = straightV * math.sin(math.radians(tilt))


    # calculate for the velocity and altitude
    if h > 1:  # if the rocket is above 1m, add the gravity
        vz = vz - (gravity * dt)
    h = h + (vz * dt)  # calculate the height
    t = t + dt  # increment time

    temp = aerodynamic_heating(vz, h)  # calculate the aerodynamic heating

    # find time of engine run-out
    if (impulse < 1 and loop == False):
        ter = t
        loop = True
    # find the max altitude
    if h > hm:  # if the current height is greater than the max height, set the max height to the current height
        hm = h
        ht = t
        ad = density
    # find the max temperature
    if temp > mt:
        mt = temp
        mtt = t
        mtpv = vz
        mtpa = h
    # find the max velocity
    if vz > vm:  # if the current velocity is greater than the max velocity, set the max velocity to the current velocity
        vm = vz
        vmt = t
        mva = h
        pmva = density
        tmv = temp
    # find the max drag force
    if pd < drag_force_N:
        pd = drag_force_N
        pdd = t
        pda = h
        pddt = t

    # Output the conditions of the simulation each cycle
    print("Stage 1:", "Time:", "{:.2f}".format(t) + "s", " Velocity:", "{:.1f}".format(vz) + "m/s", " Height:",
          str(int(h)) + "m", " Impulse:", str(int(impulse)) + "Ns",
          "(" + "{:.2f}".format((impulse / start_impulse) * 100) + "%)" +
          " Drag Force:", "{:.2f}".format(drag_force_N) + "N ", "Density:", "{:.7f}".format(density) + "kg/m^3",
          "Temp:",
          str(int(temp)) + "F", end=" ")

    if impulse > 0:  # Calculate the Thrust to Weight Ratio
        print("TWR:", "{:.2f}".format((Ns / (total_mass * gravity))), "Theoretical Acceleration:",
              "{:.1f}".format(((Ns / total_mass) - 9.8)) + "m/s^2", "True Acceleration:",
              "{:.1f}".format((((Ns - drag_force_N) / total_mass) - 9.8)) + "m/s^2")
    if impulse == 0:  # If the rocket is out of fuel, then there will be no thrust
        print("TWR: 0", "Theoretical Acceleration: 0m/s^2", "True Acceleration: 0m/s^2")

    if (((Ns - drag_force_N) / total_mass) - 9.8) > pa:
        pa = (((Ns - drag_force_N) / total_mass) - 9.8)
        pat = t

    # check if it is time to deploy the parachute
    if t >= Parachute_delay:
        # Only run this once
        Parachute_deployed = True
        Stage1_area = Parachute_area
        drag_coefficient_1 = Parachute_drag_coefficient
        Parachute_deployment_velocity = vz
        # odds of parachute failure go up exponentially as the parachute is deployed at higher speeds

        Parachute_jolt_force = drag_force_N  # jolt force upon parachute deployment
        Odds_of_parachute_failure = (drag_force_N * 0.8) - 0.5  # odds of parachute failure are proportional to the drag force
        if Odds_of_parachute_failure > 100:
            Odds_of_parachute_failure = 100
            Parachute_deployed = False
        if Odds_of_parachute_failure < 0:
            Odds_of_parachute_failure = 0
        Parachute_deployment_altitude = h
        time_deployed = t - dt
        Parachute_delay = 999999999999999999
    if Parachute_deployed:
        vx = wind_speed

iv = vz  # velocity of impact, m/s
impact_time = t  # time that the rocket hits the ground, s
end_mass = total_mass  # kg

# Give the most important information, Apogee, max Velocity, and max Temperature for both stages
print("")
print("")
print("Stage 1:")
print("Max Altitude is:", str(int(hm)) + "m,", str(int(hm * 3.28)) + "ft", "at time:", "{:.3f}".format(ht) + ", Density:",
      "{:.7f}".format(ad) + "kg/m^3")
print("Max Velocity is:", str(int(vm)) + "m/s", "at time:", "{:.3f}".format(vmt) + ", Altitude:",
      str(int(mva)) + "m" + ", Density:", "{:.7f}".format(pmva) + "kg/m^3", "Temperature:", str(int(tmv)) + "F")
print("Max Temperature is:", "{:.2f}".format(mt) + "F" + " at time:", "{:.3f}".format(mtt) + ", Velocity:",
      str(int(mtpv)) + "m/s" + ", Altitude:", str(int(mtpa)) + "m")
print("Max Drag Force is:", "{:.2f}".format(pd) + "N" + " at time:", "{:.3f}".format(pddt) + ", Altitude:", str(int(pda)) +  "m")
print("Max Acceleration is:", "{:.2f}".format(pa) + " m/s^2,", "{:.2f}".format(pa/9.81) + " g's", "at time:", "{:.2f}".format(pat))
print("Starting Mass:", "{:.3f}".format(start_mass) + "kg" + ", End Mass:", "{:.3f}".format(end_mass) + "kg")
print("Time of engine fuel run-out:", "{:.3f}".format(ter) + "s")
print("Time of impact with ground:", "{:.1f}".format(impact_time) + "s, " + "{:.1f}".format(impact_time / 60),
      "minutes")
print("Velocity of impact: " + "{:.1f}".format(iv) + "m/s")
print("")
print("")
print("Parachute Deployment:")
if Parachute_deployed:
    print("Parachute Deployment Altitude:", "{:.2f}".format(Parachute_deployment_altitude) + "m")
    print("Parachute Deployment Time:", "{:.2f}".format(time_deployed) + "s")
    print("Parachute Deployment Velocity:", "{:.2f}".format(Parachute_deployment_velocity) + "m/s")
    print("Parachute Odds of Failure:", "{:.2f}".format(Odds_of_parachute_failure) + "%," + " Jolt Force:", "{:.2f}".format(Parachute_jolt_force) + "N")
else:
    print("Parachute failed to deploy")