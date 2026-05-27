# import mesa_reader 
import mesa_reader as mr 
import matplotlib.pyplot as plt 



# load entire LOG directory information
l = mr.MesaLogDir('./LOGS')
# grab the last profile
p = l.profile_data()

# this works even if you only have logRho and logT!
plt.loglog(p.Rho, p.T)
plt.xlabel("Density")
plt.ylabel("Temperature")
plt.show() 