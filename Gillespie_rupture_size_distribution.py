import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import math
import random
from termcolor import colored

lam = 0.643111 #Vegetative bacteria replication rate
mu = 1.637989 #Vegetative bacteria death rate
gamma = 0.043792 #Cell rupture rate
ngb_death = 0.003502 #Newly germinated bacteria death rate
ga = 0.894274 #Type A germination rate
   
xrange = 100

def rup_dist(lam,mu,gamma,ngb_death,ga,nsim):
   #Number of Gillespie simualtions
   numsim = nsim
   
   #Number of initial spores
   init_spore = 1
   
   num_ngb_death = 0
   num_vb_death = 0
   num_rupture = 0
   num_continue = 0
   time_to_rupture = []
   time_to_recovery = []
   rupture_size = []
   for run in range(numsim):
      
      #val = np.random.uniform()
      
      x1 = init_spore #Number of initial spores
      x1t = [init_spore]
      x2 = 0 #Number of initial newly germinated bacteria
      x2t = [0] 
      x3 = 0 #Number of initial vegetative bacteria
      x3t = [0]
      
      #Time to run process until (should be large enough that all Gillespie simulations end)
      t_last = 1000
      t = 0 #Initial time
      tt = [0] #Initialise list of times
         
      #Gillespie algorithm
      while t < t_last:
            
         germin1 = ga*x1 #Rate of spore to NGB
         germin2 = ga*x2 #Rate of NGB to vegetative bacteria
         death1 = ngb_death*x2 #Rate of NGB death
         birth = lam*x3 #Rate of vegetative bacteria replication
         death2 = mu*x3 #Rate of vegetative bacterial death
         rupture = gamma*x3 #Rate of cell rupture
      
         ratesum = germin1 + germin2 + death1 + birth + death2 + rupture
         urv=random.random()
         #Choose the event which occurs
         if urv < germin1/ratesum:
            x1-=1
            x2+=1
         elif urv < (germin1+germin2)/ratesum:
            x2-=1
            x3+=1
         elif urv < (germin1+germin2+death1)/ratesum:
            x2-=1   
            num_ngb_death += 1 #Count number of sims that end with NGB death
            t+=-math.log(random.random())/ratesum            
            time_to_recovery.append(t)    
            rupture_size.append(0)                    
            break
         elif urv < (germin1+germin2+death1+birth)/ratesum:
            x3+=1
         elif urv < (germin1+germin2+death1+birth+death2)/ratesum:
            x3-=1
            if x3 == 0 and x1 != init_spore and x2 != 1:
               num_vb_death += 1 #Count number of sims that end with vegetative bacterial death
               time_to_recovery.append(t)       
               rupture_size.append(0)
               break #Break the loop when no bacteria remain
         else:
            num_rupture += 1 #Count number of sims that end rupture
            t+=-math.log(random.random())/ratesum
            time_to_rupture.append(t)            
            rupture_size.append(x3)
            break
      
         t+=-math.log(random.random())/ratesum
         tt.append(t)
         x1t.append(x1)
         x2t.append(x2)
         x3t.append(x3)   
      
      if tt[-1] > t_last:
         num_continue += 1 #Count number of sims that are ongoing at time t_last
   
   if num_continue > 0:
      print(colored('Gillespie simulations still running - increase t_last', 'red', attrs=['bold']))
      print(' ')
      
   #Compute expected rupture size from Gillespie
   mrs = max(rupture_size)
   rupture_dist_gillespie = []
   for nn in range(mrs+1):
      rupture_dist_gillespie.append(rupture_size.count(nn)/numsim)

   return rupture_dist_gillespie
   
plot = rup_dist(lam,mu,gamma,ngb_death,ga,100)
if len(plot) < 10:
   for i in range(10-len(plot)):
      plot.append(0)

#############################################################################

fig, ax = plt.subplots(figsize=(5.5,6.5))

ax.bar(range(10), plot[0:10], color='gray', alpha=0.8)
plt.subplots_adjust(bottom=0.55,top=0.95)

plt.title('Gillespie rupture size distribution', fontsize=12)
plt.ylabel('Probability', fontsize=12)
plt.xlabel('Number of bacteria released', fontsize=12)

yvals = [0.42,0.34,0.26,0.18,0.1,0.02]
labels = ["Nsim",r"$\lambda$",r"$\mu$",r"$\gamma$",r"$\tilde{\mu}$",r"$g$"]
valmins = [1,-1,0,-2,-3,-1]
valmaxs = [4,1,1,0,-2,1]
valinits = [200,0.643111,1.637989,0.043792,0.003502,0.894274]
cols = ['C5','C0','C1','C2','C3','C4']

pdict = {}
for ii in range(6):
   sax = plt.axes([0.16, yvals[ii], 0.65, 0.04]) #left, bottom, width, height
   pdict["key%s" %ii] = Slider(ax=sax,label=labels[ii],valmin=valmins[ii],valmax=valmaxs[ii],valinit=np.log10(valinits[ii]),color=cols[ii])
   pdict["key%s" %ii].vline.set_color('black')
   pdict["key%s" %ii].vline.set_linewidth(3)
   pdict["key%s" %ii].label.set_size(14)
   if ii == 0:
      pdict["key%s" %ii].valtext.set_text(int(valinits[ii]))
   else:
      pdict["key%s" %ii].valtext.set_text(np.round(valinits[ii],3))

def update(val):
      
   pars = []
   for ii in range(6):
      amp = pdict["key%s" %ii].val
      pars.append(amp)
      if ii == 0:
         pdict["key%s" %ii].valtext.set_text(int(10**amp))
      else:
         pdict["key%s" %ii].valtext.set_text(np.round(10**amp,3))
   plot = rup_dist(10**pars[1],10**pars[2],10**pars[3],10**pars[4],10**pars[5],int(10**pars[0]))
   if len(plot) < 10:
      for i in range(10-len(plot)):
         plot.append(0)   
   xind = 10
   for xx in range(len(plot)):
      if plot[xx] > 0.001:
         xind = xx
   if xind <= 10:
      xrange = 10
   elif xind <= 100:
      xrange = xind
   else:
      xrange = 100
   ax.clear()
   ax.bar(range(xrange), plot[0:xrange], color='gray', alpha=0.8)
   ax.set_title('Gillespie rupture size distribution', fontsize=12)
   ax.set_ylabel('Probability', fontsize=12)
   ax.set_xlabel('Number of bacteria released', fontsize=12)   
   fig.canvas.draw_idle()
   
for ii in range(6):   
   pdict["key%s" %ii].on_changed(update)

resetax = plt.axes([0.79, 0.9, 0.1, 0.04])
button = Button(resetax, 'Reset', hovercolor='0.975')

def reset(event):
   for ii in range(6):
      pdict["key%s" %ii].reset()

button.on_clicked(reset)
button.label.set_fontsize(12)

plt.show()

