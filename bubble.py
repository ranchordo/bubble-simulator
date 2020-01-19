import matplotlib.pyplot as plt
import numpy as np
dens_water=1000 #Density of water, KG/M^3
dens_air_init=1.225 #Initial uncompressed density of air, KG/M^3
air_spec_gas_constant=287.058 #Specific gas constant of air
drag_coeff=0.57 #Drag coefficient of a sphere
bub_radius_init=0.03 #Initial bubble radius #M
g=9.80665 #Gravitational acceleration, M/S^2
temperature=25 #C
tank_depth=100 #M
period=0.0001 #Iterative period, Smaller --> more accurate
print("Reciprocal of Period: "+str(1/period))

plt.style.use('dark_background')

#CALC STUFF
temp=temperature+273.15 #Calculate temperature in K
bub_v_init=(4/3)*np.pi*(bub_radius_init**3) ##calculate initial bubble volume
bub_m_init=dens_air_init*bub_v_init #Calculate initial bubble mass
p=(dens_water*g*tank_depth) #Calculate pressure at bottom
dens_air_bottom=p/(air_spec_gas_constant*temp)#Density of air at the bottom of the tank

#Initialize some arrays
hs=[] #Height
drs=[] #Drag force
ds=[] #Air density
ts=[] #Time
rs=[] #Radius
fs=[] #Bubble acceleration
vs=[] #Velocity
bs=[] #Buoyancy force
#Extract sign of a value
def sign(a):
    if a<0:
        f=-1
    if a>0:
        f=1
    if a==0:
        f=0
    return f
#Main function:
def main():
    #Initialize some variables
    h=0
    v=0
    t=0
    drag=0
    dens_air=dens_air_bottom
    bub_radius=bub_radius_init
    bub_acc=0
    buoyancy=0
    #While bubble is submerged:
    while h<tank_depth:
        #Write data to lists
        hs.append(h)
        drs.append(drag)
        ds.append(dens_air)
        ts.append(t)
        rs.append(bub_radius)
        fs.append(bub_acc)
        vs.append(v)
        bs.append(buoyancy)
        depth=tank_depth-h
        p=(dens_water*g*depth)+101325
        dens_air=p/(air_spec_gas_constant*temp)
        bub_v=bub_m_init/dens_air
        bub_radius=((3*bub_v)/(4*np.pi))**(1/3)
        bub_A=np.pi*(bub_radius**2)
        buoyancy=dens_water*bub_v*g
        grav=g*bub_m_init
        drag=0.5*bub_A*dens_water*drag_coeff*((v**2)*sign(v))
        total_force=buoyancy-(grav+drag)
        bub_acc=total_force/bub_m_init

        delta_h=((1/2)*bub_acc*(period**2))+(v*period)
        delta_v=bub_acc*period

        if delta_h < 0:
            print("DELTA H < 0 AT TIME: "+str(t)) #Bubble is going down
            print("VAR DUMP:") #Dump variables
            print("  Velocity: "+str(v+delta_v))
            if(v+delta_v > 3e8):
                print("First object to break light speed!") #Joke message if greater than lightspeed
            print("  Bubble Total Force: "+str(total_force))
            print("  Total Force Breakdown:")
            print("    Buoyancy Force: "+str(buoyancy))
            print("    Grav Force: "+str(grav))
            print("    Drag Force: "+str(drag))
            print("  Air Density at "+str(depth)+" Meters: "+str(dens_air))
            print("  Pressure at "+str(depth)+" Meters: "+str(p))
            quit()
        #Update time, height, and velocity
        t+=period
        h+=delta_h
        v+=delta_v
print("Definitions done, running main loop...")
main() #Run main
print("Finshed. Length of data was "+str(len(hs)))

print("Starting postprocessing...")
#Convert the data to np array
hsa=np.array(hs)
drsa=np.array(drs)
dsa=np.array(ds)
tsa=np.array(ts)
rsa=np.array(rs)
fsa=np.array(fs)
vsa=np.array(vs)
bsa=np.array(bs)

#Display the plots
fig=plt.figure()
fig.suptitle('Height / Time')
plt.xlim(0,max(tsa))
plt.ylim(0,max(hsa))
plt.plot(tsa,hsa)

fig2=plt.figure()
fig2.suptitle('Bubble Radius / Time')
plt.xlim(0,max(tsa))
plt.ylim(0,max(rsa))
plt.plot(tsa,rsa)

fig3=plt.figure()
fig3.suptitle('Air Density / Time')
plt.xlim(0,max(tsa))
plt.ylim(0,max(dsa))
plt.plot(tsa,dsa)

fig4=plt.figure()
fig4.suptitle('Bubble Velocity / Time')
plt.xlim(0,max(tsa))
plt.ylim(0,max(vsa))
plt.plot(tsa,vsa)

fig5=plt.figure()
fig5.suptitle('Buoyancy Force / Time')
plt.xlim(0,max(tsa))
plt.ylim(0,max(bsa))
plt.plot(tsa,bsa)

fig6=plt.figure()
fig6.suptitle('Drag Force / Time')
plt.xlim(0,max(tsa))
plt.ylim(0,max(drsa))
plt.plot(tsa,drsa)

fig7=plt.figure()
fig7.suptitle('Bubble Upward Acceleration / Time (Zoom into left side)')
plt.xlim(0,max(tsa))
plt.ylim(0,max(fsa))
plt.plot(tsa,fsa)
plt.show()
