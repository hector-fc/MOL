
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D
import scipy.interpolate 
from mpl_toolkits.mplot3d import Axes3D


# def  funcion solution
# Parameters
#style.use('dark_backgound')
##  densities ....

## solitons 
soli = np.genfromtxt('5000.txt',skip_header=0)
solc = np.genfromtxt('5001.txt',skip_header=0)
solf = np.genfromtxt('5002.txt',skip_header=0)

qdeni = np.genfromtxt('7000.txt',skip_header=0)
qdenc = np.genfromtxt('7001.txt',skip_header=0)
qdenf = np.genfromtxt('7002.txt',skip_header=0)

gdeni = np.genfromtxt('8000.txt',skip_header=0)
gdenc = np.genfromtxt('8001.txt',skip_header=0)
gdenf = np.genfromtxt('8002.txt',skip_header=0)

ener = np.genfromtxt('energy.txt',skip_header=0)   

ax = -20.0
bx =  20.0

ay = -10.0
by =  10.0


tmax = ener[:,0].max()


fig1 = plt.figure(1)
#====================================
#Soliton 
#
plt.subplot(211)
plt.plot(soli[:,0],soli[:,1],'-g',lw =1.0,label='$t_i$')
plt.plot(solc[:,0],solc[:,1],'-b',lw =1.0,label='$t_c$')
plt.plot(solf[:,0],solf[:,1],'-r',lw =1.0,label='$t_f$')

plt.axis([ax,bx,ay,by])
plt.grid(True)
# legend
plt.legend(bbox_to_anchor=(1.0,1.0),loc=1,fontsize=10, borderaxespad = 0.)


#====================================
#

plt.subplot(212)
'''
plt.plot(ener[:,0],ener[:,1],'-r',lw =1.0,label=r'$E(t)\approx %0.2f$'% np.mean(ener[:,1]))
plt.plot(ener[:,0],ener[:,2],'-b',lw =1.0,label=r'$P(t)\approx %0.2f$'% np.mean(ener[:,2]))
plt.plot(ener[:,0],ener[:,3],'-g',lw =1.0,label=r'$a_{+}(t)\approx %0.2f$'% np.mean(ener[:,3]))
plt.plot(ener[:,0],ener[:,4],'-k',lw =1.0,label=r'$a_{-}(t)\approx %0.2f$'% np.mean(ener[:,4]))
'''

plt.plot(ener[:,0],ener[:,1],'-r',lw =1.0,label=r'$E(t)$' )
plt.plot(ener[:,0],ener[:,2],'-b',lw =1.0,label=r'$P(t)$')

plt.axis([0.0,tmax,-100,100])
plt.grid(True)
# legend
plt.legend(bbox_to_anchor=(1.0,1.0),loc=1,fontsize=10, borderaxespad = 0.)

#====================================
fig1.subplots_adjust(hspace=.7)
fig1.show()
#====================================
#===================================

fig2 = plt.figure(2)
#====================================
#Soliton 
#
plt.subplot(211)
plt.plot(qdeni[:,0],qdeni[:,1],'-g',lw =1.0,label='$t_i$')
plt.plot(qdenc[:,0],qdenc[:,1],'-b',lw =1.0,label='$t_c$')
plt.plot(qdenf[:,0],qdenf[:,1],'-r',lw =1.0,label='$t_f$')

plt.axis([ax,bx,-50,50])
plt.grid(True)
# legend
plt.legend(bbox_to_anchor=(1.0,1.0),loc=1,fontsize=10, borderaxespad = 0.)

#====================================
#

plt.subplot(212)

'''
plt.plot(ener[:,0],ener[:,1],'-r',lw =1.0,label=r'$E(t)\approx %0.2f$'% np.mean(ener[:,1]))
plt.plot(ener[:,0],ener[:,2],'-b',lw =1.0,label=r'$P(t)\approx %0.2f$'% np.mean(ener[:,2]))
plt.plot(ener[:,0],ener[:,3],'-g',lw =1.0,label=r'$a_{+}(t)\approx %0.2f$'% np.mean(ener[:,3]))
plt.plot(ener[:,0],ener[:,4],'-k',lw =1.0,label=r'$a_{-}(t)\approx %0.2f$'% np.mean(ener[:,4]))
'''

plt.plot(ener[:,0],ener[:,3],'-g',lw =1.0,label=r'$a_{+}(t)$')

plt.axis([0.0,tmax,-100,100])
plt.grid(True)
# legend
plt.legend(bbox_to_anchor=(1.0,1.0),loc=1,fontsize=10, borderaxespad = 0.)

#====================================
fig2.subplots_adjust(hspace=.7)
fig2.show()


fig3 = plt.figure(3)
#====================================
#Soliton 
#
plt.subplot(211)
plt.plot(qdeni[:,0],qdeni[:,2],'-g',lw =1.0,label='$t_i$')
plt.plot(qdenc[:,0],qdenc[:,2],'-b',lw =1.0,label='$t_c$')
plt.plot(qdenf[:,0],qdenf[:,2],'-r',lw =1.0,label='$t_f$')

plt.axis([ax,bx,-200,200])
plt.grid(True)
# legend
plt.legend(bbox_to_anchor=(1.0,1.0),loc=1,fontsize=10, borderaxespad = 0.)

#====================================
#

plt.subplot(212)

plt.plot(ener[:,0],ener[:,4],'-k',lw =1.0,label=r'$a_{-}(t)$')

plt.axis([0.0,tmax,-100,100])
plt.grid(True)
# legend
plt.legend(bbox_to_anchor=(1.0,1.0),loc=1,fontsize=10, borderaxespad = 0.)

#====================================
fig3.subplots_adjust(hspace=.7)
fig3.show()


fig4 = plt.figure(4)
#====================================
#Soliton 
#
plt.subplot(211)
plt.plot(gdeni[:,0],gdeni[:,1],'-g',lw =1.0,label='$t_i$')
plt.plot(gdenc[:,0],gdenc[:,1],'-b',lw =1.0,label='$t_c$')
plt.plot(gdenf[:,0],gdenf[:,1],'-r',lw =1.0,label='$t_f$')

plt.axis([ax,bx,-200,200])
plt.grid(True)
# legend
plt.legend(bbox_to_anchor=(1.0,1.0),loc=1,fontsize=10, borderaxespad = 0.)


#====================================
#

plt.subplot(212)
plt.plot(ener[:,0],ener[:,5],'-k',lw =1.0,label=r'$\tilde{q}_{2}(t)$')
plt.axis([0.0,tmax,-500,500])
plt.grid(True)
# legend
plt.legend(bbox_to_anchor=(1.0,1.0),loc=1,fontsize=10, borderaxespad = 0.)

#====================================
fig4.subplots_adjust(hspace=.7)
fig4.show()


fig5 = plt.figure(5)
#====================================
#Soliton 
#
plt.subplot(211)
plt.plot(gdeni[:,0],gdeni[:,2],'-g',lw =1.0,label='$t_i$')
plt.plot(gdenc[:,0],gdenc[:,2],'-b',lw =1.0,label='$t_c$')
plt.plot(gdenf[:,0],gdenf[:,2],'-r',lw =1.0,label='$t_f$')

plt.axis([ax,bx,-200,200])
plt.grid(True)
# legend
plt.legend(bbox_to_anchor=(1.0,1.0),loc=1,fontsize=10, borderaxespad = 0.)

#====================================
#
plt.subplot(212)
plt.plot(ener[:,0],ener[:,6],'-k',lw =1.0,label=r'$q_{2}(t)$')
plt.axis([0.0,tmax,-100,100])
plt.grid(True)
# legend
plt.legend(bbox_to_anchor=(1.0,1.0),loc=1,fontsize=10, borderaxespad = 0.)

#====================================
fig5.subplots_adjust(hspace=.7)
fig5.show()


fig1.savefig('kisol.eps', dpi=300, format='eps')
#fig2.savefig('qp.eps')
#fig3.savefig('qn.eps')
fig4.savefig('q2p.eps')
fig5.savefig('q2n.eps')

enter = input("enter")
