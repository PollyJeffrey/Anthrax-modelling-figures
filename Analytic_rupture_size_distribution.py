import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

lam = 0.643111 #Vegetative bacteria replication rate
mu = 1.637989 #Vegetative bacteria death rate
gamma = 0.043792 #Cell rupture rate
ngb_death = 0.003502 #Newly germinated bacteria death rate
ga = 0.894274 #Type A germination rate
   
xrange = 100

def rup_dist(lam,mu,gamma,ngb_death,ga):
   #Compute analytic expected rupture size
   aa = ((lam+mu+gamma)-np.sqrt((lam+mu+gamma)**2-4*mu*lam))/(2*lam)
   bb = ((lam+mu+gamma)+np.sqrt((lam+mu+gamma)**2-4*mu*lam))/(2*lam)
   
   def R1n(n):
      if n == 0:
         r1n = (ngb_death + (ga*aa))/(ngb_death + ga)
      else:
         r1n = ((1-aa)*(bb-1))/(bb**n)
      return r1n

   rupt_dist_analytic = []
   for nn in range(xrange):
      first = ga/(ngb_death+ga)
      if nn == 0:
         value = R1n(nn)
      else:
         value = R1n(nn)*first
      rupt_dist_analytic.append(value)   
   return rupt_dist_analytic
   
plot = rup_dist(lam,mu,gamma,ngb_death,ga)

#############################################################################

fig, ax = plt.subplots(figsize=(5.5,6.5))

ax.bar(range(10), plot[0:10], color='gray', alpha=0.8)
plt.subplots_adjust(bottom=0.5,top=0.95)

plt.title('Analytic rupture size distribution', fontsize=12)
plt.ylabel('Probability', fontsize=12)
plt.xlabel('Number of bacteria released', fontsize=12)

yvals = [0.34,0.26,0.18,0.1,0.02]
labels = [r"$\lambda$",r"$\mu$",r"$\gamma$",r"$\tilde{\mu}$",r"$g$"]
valmins = [-1,0,-2,-3,-1]
valmaxs = [1,1,0,-2,1]
valinits = [0.643111,1.637989,0.043792,0.003502,0.894274]
cols = ['C0','C1','C2','C3','C4']

pdict = {}
for ii in range(5):
   sax = plt.axes([0.16, yvals[ii], 0.65, 0.04]) #left, bottom, width, height
   pdict["key%s" %ii] = Slider(ax=sax,label=labels[ii],valmin=valmins[ii],valmax=valmaxs[ii],valinit=np.log10(valinits[ii]),color=cols[ii])
   pdict["key%s" %ii].vline.set_color('black')
   pdict["key%s" %ii].vline.set_linewidth(3)
   pdict["key%s" %ii].label.set_size(14)
   pdict["key%s" %ii].valtext.set_text(np.round(valinits[ii],3))

def update(val):
   pars = []
   for ii in range(5):
      amp = pdict["key%s" %ii].val
      pars.append(amp)
      pdict["key%s" %ii].valtext.set_text(np.round(10**amp,3))
   plot = rup_dist(10**pars[0],10**pars[1],10**pars[2],10**pars[3],10**pars[4])
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
   ax.set_title('Analytic rupture size distribution', fontsize=12)
   ax.set_ylabel('Probability', fontsize=12)
   ax.set_xlabel('Number of bacteria released', fontsize=12)   
   fig.canvas.draw_idle()
   
for ii in range(5):   
   pdict["key%s" %ii].on_changed(update)

resetax = plt.axes([0.79, 0.9, 0.1, 0.04])
button = Button(resetax, 'Reset', hovercolor='0.975')

def reset(event):
   for ii in range(5):
      pdict["key%s" %ii].reset()

button.on_clicked(reset)
button.label.set_fontsize(12)

plt.show()