#Rocket Sim 3 - By: Reese Kimmel - 10/12/2023

import math

# Rocket Parameters
Stage1_area = 0.0027  # m^2
drag_coefficient = 0.85  # drag coefficient of the rocket

fuel_mass = 0.020  # kg
payload_mass = 0.0  # kg, 1190 grams of payload (the second stage)
misc_mass = 0.0  # kg, miscellaneous bits of extra mass, fuel tubes, nozzles, etc...
body_mass = 0.22  # kg, mass of the body of the rocket
fin_coefficient = 0.01

impulse = 20  # The total impulse of the engine, Ns

Parachute_delay = 8.7  # Delay between motor ignition and parachute deployment, s
Parachute_drag_coefficient = 4.05
Parachute_area = 0.027

stage1_loss_percentage = 0  # percent of impulse lost to inefficiencies (good for approximating misc. losses)

#Stage 0 Parameters
rail_length = 1  # length of the launch rail, m
forced_tilt = 0  # the tilt of the rocket forced by the launch rail, degrees
h = 0  # height, m


#Simulation Parameters
gravity = 9.81  # m/s^2
dt = 0.0001  # time step, s
height_offset = 250  # starting altitude above sea level
wind_speed = 1  # wind speed, m/s
output_frequency = 1000


#  Misc. variables
t = 0  # time, s
vx = 0  # velocity on the X-Axis, m/s
vz = 0  # velocity on the Z-Axis, m/s
tilt = 0  # tilt of the rocket, degrees
dx = 0  # distance on the X-Axis, m
rocket_temp = 0  # temperature of rocket
payload_mass += misc_mass
payload_mass += body_mass
total_mass = payload_mass  # kg
total_mass += fuel_mass
counter = 0  # used for the output frequency
impulse = impulse * (1 - (stage1_loss_percentage / 100))  # account for losses
sea_level_speed = wind_speed  # wind speed at sea level, m/s
acceleration = 0  # m/s^2
max_tilt = 0  # max tilt of the rocket, degrees
impact_time = 0  # time of rocket impact with ground
start_mass = total_mass  # kg
start_impulse = impulse  # Ns, defined to calculate the fuel mass
starting_fuel_mass = fuel_mass  # kg, defined to calculate the decreasing fuel mass
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
pdv = 0  # velocity at parachute deployment, m/s
mtpa = 0  # altitude of max temperature, m
ma = 0  # max acceleration, m/s^2
ter = 0  # time of engine run-out, s
pd = 0  # peak drag, N
drag_force = 0  # N
Deployed = False
Jolt_Force = 0  # N

print("Total mass: " + str(total_mass) + " kg")


# Define functions
def aerodynamic_heating(velocity, altitude, temp):
    altitude = altitude + height_offset
    # the aerodynamic heating is modeled as an exponential decay
    velocity = math.sqrt(velocity ** 2)
    ambient = 72 - (11.7 * (altitude / 1000))  # approximate ambient temperature at altitude, F
    density = atmosphere(altitude)
    temp_diff = temp - ambient
    heating = ((velocity * 1.225) / 100) - (temp_diff * 0.1)
    return heating


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

def drag(density, drag_coefficient, Stage1_area, straightV):
    # calculate for the drag force
    rho = density  # kg/m^3
    Cd = drag_coefficient  # drag coefficient
    A = Stage1_area  # M^2

    drag_force = 0.5 * rho * Cd * A * (straightV ** 2)  # drag force in newtons
    return drag_force

def Celcius(temp):
    temp -= 32
    temp *= (5/9)
    return temp

def E12(time):
    if time <= 0.25:
        return (293 * (time * time)) + (45.4 * time) + 0.0357
    if time > 0.25 and time <= 0.55:
        return (-1778 * (time * time * time)) + (2638 * (time * time)) - (1307 * time) + 228
    if time > 0.55 and time <= 2.9:
        return 10
    else:
        return 0
def D12(time):
    if time <= 0.25:
        return (293 * (time * time)) + (45.4 * time) + 0.0357
    if time > 0.25 and time <= 0.55:
        return (-1778 * (time * time * time)) + (2638 * (time * time)) - (1307 * time) + 228
    if time > 0.55 and time <= 1.65:
        return 10
    else:
        return 0


# the main loop is defined, where all the main physics is done
while h >= 0:
    N = D12(t)  # calculate the thrust of the rocket at the current time

    counter += 1

    density = atmosphere(h)  # calculate the density of the atmosphere at the current height

    wind_speed = wind(h, sea_level_speed)

    if forced_tilt == 0:
        if h > rail_length:
            tilt = tilt2(wind_speed, vz)
        max_tilt = max(max_tilt, tilt)
        tilt = max_tilt
    else:
        tilt = forced_tilt

    # find the x velocity of the rocket
    if impulse > 0:
        vx += math.sin(math.radians(tilt)) * (N / total_mass) * dt

    # calculate the distance the rocket has traveled
    dx += (vx * dt)
    z_percent = math.cos(math.radians(tilt))
    rho = density  # kg/m^3
    fin_authority = fin_coefficient * vz ** 2  # drag force in N

    # calculate for the rocket engine spending fuel
    if impulse > 0:  # if the rocket has fuel, subtract the fuel
        impulse = impulse - N * dt
        fuel_mass = starting_fuel_mass * (impulse / start_impulse)  # calculate the fuel mass
        total_mass = fuel_mass + payload_mass
        vz = vz + ((N / total_mass) * dt) * z_percent  # calculate the velocity


    if impulse < 0:  # fix a bug where the impulse would go under 0
        impulse = 0

    # calculate the drag force
    straightV = math.sqrt(vx ** 2 + vz ** 2)
    drag_force_N = drag(density, drag_coefficient, Stage1_area, straightV)


    if straightV > 0:  # if the rocket is moving up, subtract the drag force
        straightV -= ((drag_force_N / total_mass) * dt)
    if vx < 0:  # if the rocket is moving down, add the drag force
        straightV += ((drag_force_N / total_mass) * dt)


    if vz < 0:
        vz = -1 * straightV * math.cos(math.radians(tilt))  # calculate the vertical velocity
    else:
        vz = straightV * math.cos(math.radians(tilt))
    if vx < 0:
        vx = -1 * straightV * math.sin(math.radians(tilt))  # calculate the horizontal velocity
    else:
        vx = straightV * math.sin(math.radians(tilt))


    # calculate for the velocity and altitude
    if h > 1:  # if the rocket is above 1m, add the gravity
        vz = vz - (gravity * dt)

    h += (vz * dt)  # calculate the height
    t += dt  # increment time
    rocket_temp += aerodynamic_heating(vz, h, rocket_temp)  # calculate the aerodynamic heating

    #  output simulation conditions
    if counter == output_frequency:
        # Output the conditions of the simulation each cycle
        print("Stage 1:", "Time:", "{:.2f}".format(t) + "s", " Velocity:",
              "{:.1f}".format(vz) + "m/s (" + "{:.1f}".format(vz * 2.2) + "mi/h)", " Height:",
              str(int(h)) + "m (" + "{:.1f}".format(h * 3.28) + "ft)", " Impulse:", str(int(impulse)) + "Ns",
              "(" + "{:.2f}".format((impulse / start_impulse) * 100) + "%)" +
              " Drag Force:", "{:.2f}".format(drag_force_N) + "N ", "Density:", "{:.3f}".format(density) + "kg/m^3",
              "Temperature:",
              str(int(Celcius(rocket_temp))) + "°C, " + str(int(rocket_temp)) + "°F")

        print("Rocket drifting parameters:", "Tilt:", "{:.2f}".format(tilt) + "°", "Velocity:", "{:.2f}".format(vx) + "m/s (" + "{:.1f}".format(vx * 2.2) + "mi/h)",
              "Distance:", str(int(dx)) + "m (" + "{:.1f}".format(dx * 3.28) + "ft)", "Fin Authority:", "{:.1f}".format(fin_authority * 50) + "%")

        if impulse > 0:  # Calculate the Thrust to Weight Ratio
            print("TWR:", "{:.2f}".format((N / (total_mass * gravity))), "Theoretical Acceleration:",
                  "{:.1f}".format(((N / total_mass) - 9.8)) + "m/s^2", "True Acceleration:",
                  "{:.1f}".format((((N - drag_force_N) / total_mass) - 9.8)) + "m/s^2")
        if impulse == 0:  # If the rocket is out of fuel, then there will be no thrust
            print("TWR: 0", "Theoretical Acceleration: 0m/s^2", "True Acceleration: 0m/s^2")
        counter = 0
        print("")














    # find time of engine run-out
    if (impulse < 1 and ter == 0):
        ter = t
    # find the max altitude
    if h > hm:  # if the current height is greater than the max height, set the max height to the current height
        hm = h
        ht = t
        ad = density
    # find the max temperature
    if rocket_temp > mt:
        mt = rocket_temp
        mtt = t
        mtpv = vz
        mtpa = h
    # find the max velocity
    if vz > vm:  # if the current velocity is greater than the max velocity, set the max velocity to the current velocity
        vm = vz
        vmt = t
        mva = h
        pmva = density
        tmv = rocket_temp

    if (((N - drag_force_N) / total_mass) - 9.8) > pa:
        pa = (((N - drag_force_N) / total_mass) - 9.8)
        pat = t

    if t >= Parachute_delay:
        drag_coefficient = Parachute_drag_coefficient
        Stage1_area = Parachute_area
        pdv = straightV
        Parachute_delay = 999999999999
        Jolt_Force = drag(density, drag_coefficient, Stage1_area, straightV)
        Deployed = True
    # find the max drag force
    if pd < drag_force_N:
        if not Deployed:
            pd = drag_force_N
            pda = h
            pddt = t












iv = vz  # velocity of impact, m/s
impact_time = t  # time that the rocket hits the ground, s
end_mass = total_mass  # kg

# Give the most important information, Apogee, max Velocity, and max Temperature for both stages
print("")
print("")
print("Stage 1:")
print("Max Altitude is:", str(int(hm)) + "m,", str(int(hm * 3.28)) + "ft", "at time:", "{:.3f}".format(ht) + "s, Density:",
      "{:.3f}".format(ad) + "kg/m^3")
print("Max Velocity is:", str(int(vm)) + "m/s (" + "{:.1f}".format(vm * 2.2) + "mi/h)", "at time:", "{:.3f}".format(vmt) + "s, Altitude:",
      str(int(mva)) + "m" + ", Density:", "{:.4f}".format(pmva) + "kg/m^3,", "Temperature:", str(int(Celcius(tmv))) + "°C, " + str(int(tmv)) + "°F" )
print("Max Temperature is:", str(int(Celcius(mt))) + "°C, " + str(int(mt)) + "°F" + " at time:", "{:.3f}".format(mtt) + "s, Altitude:", str(int(mtpa)) + "m")
print("Max Drag Force is:", "{:.2f}".format(pd) + "N" + " at time:", "{:.3f}".format(pddt) + "s, Altitude:", str(int(pda)) +  "m")
print("Max Acceleration is:", "{:.2f}".format(pa) + " m/s^2,", "{:.2f}".format(pa/9.81) + " g's", "at time:", "{:.2f}".format(pat) + "s")
print("Starting Mass:", "{:.3f}".format(start_mass) + "kg" + ", End Mass:", "{:.3f}".format(end_mass) + "kg")
print("Time of engine fuel run-out:", "{:.1f}".format(ter) + "s")
print("Time of impact with ground:", "{:.1f}".format(impact_time) + "s, " + "{:.1f}".format(impact_time / 60),
      "minutes")
print("Velocity of impact: " + "{:.1f}".format(iv) + "m/s (" + "{:.1f}".format(iv * 2.2) + "mi/h)")
print("")
print("")
print("Parachute Deployment:")
if Deployed:
    print("Velocity at Parachute Deployment: " + "{:.1f}".format(pdv) + "m/s (" + "{:.1f}".format(pdv * 2.2) + "mi/h)")
    print("Jolt Force: " + "{:.1f}".format(Jolt_Force) + "N")
else:
    print("Parachute Failed to Deploy")
