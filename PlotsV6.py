# %%
## setting the plot parameters
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import sem 

#setting figure size
plt.rcParams['figure.figsize'] = [8,6]
plt.rcParams['figure.dpi'] = 120
#setting x and y ticks to true
plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = 1
plt.rcParams['ytick.left'] = plt.rcParams['ytick.labelleft'] = 1
plt.rcParams["axes.linewidth"] = 1.2
plt.rcParams['font.family']='serif'
plt.rcParams['font.size']=15
plt.rcParams["figure.autolayout"] = True

# %%
## reading the files and setting the mass and sfr variables
results = pd.read_csv('/Users/nischal/Desktop/PhD/FinalCigaleChecks/CigaleOut/XrayReffV4Lx.csv')
# results = results[results['best.universe.redshift']>0.15]
# results = results.loc[(results['M_criteria']<5) & (results['M_criteria']>0.2)]
controls = pd.read_csv('/Users/nischal/Desktop/PhD/FinalCigaleChecks/CigaleOut/ControlReffV3.csv')
# controls = controls.loc[(controls['M_criteria']<5) & (controls['M_criteria']>0.2)]
#controls = controls[controls['best.reduced_chi_square']<10]
#results = results[results['best.reduced_chi_square']<10]


#control_class = pd.read_csv('/Users/nischal/Desktop/PhD/FinalCigaleChecks/CigaleOut/ControlBPT_WHAN.csv')

control_class = pd.read_csv('/Users/nischal/Desktop/PhD/FinalCigaleChecks/CigaleOut/ControlBPTv4.csv')
control_class = control_class[control_class['major_arcsec']>2]
control_class = control_class[control_class['rMAGsdss']<20]

control_composite = control_class[control_class['sf_class']=='Composite']
control_sf = control_class[control_class['sf_class']=='SF']

compositeids = ('|'.join(control_composite['id']))
sfids = ('|'.join(control_sf['id']))

composites = controls[controls['id'].str.contains(compositeids)]
starformings = controls[controls['id'].str.contains(sfids)]

controls = pd.concat([composites,starformings])

shells1_xray = pd.read_csv('/Users/nischal/Desktop/PhD/FinalCigaleChecks/CigaleOut/CentralXrayBins.csv')
shell1xraymass = np.log10(shells1_xray['bayes.stellar.m_star'])[~np.isnan(shells1_xray['bayes.stellar.m_star'])]
shell1xraysfr = np.log10(shells1_xray['bayes.sfh.sfr']).dropna()
#shell1xraysfr = shell1xraysfr[shell1xraysfr>-3.5]

# %% defining the shells
##setting shell mass and sfr variables

#for xrays
shell1 = results[results.id.str.endswith('-0.25')]
shell2 = results[results.id.str.endswith('-0.5')]
shell3 = results[results.id.str.endswith('-0.75')]
shell4 = results[results.id.str.endswith('-1.0')]
shell5 = results[results.id.str.endswith('-1.25')]
shell6 = results[results.id.str.endswith('-1.5')]
shell7 = results[results.id.str.endswith('-1.75')]
shell8 = results[results.id.str.endswith('-2.0')]

#for control
control1 = controls[controls.id.str.endswith('-0.25')]
control2 = controls[controls.id.str.endswith('-0.5')]
control3 = controls[controls.id.str.endswith('-0.75')]
control4 = controls[controls.id.str.endswith('-1.0')]
control5 = controls[controls.id.str.endswith('-1.25')]
control6 = controls[controls.id.str.endswith('-1.5')]
control7 = controls[controls.id.str.endswith('-1.75')]
control8 = controls[controls.id.str.endswith('-2.0')]

#for composite
composite1 = composites[composites.id.str.endswith('-0.25')]
composite2 = composites[composites.id.str.endswith('-0.5')]
composite3 = composites[composites.id.str.endswith('-0.75')]
composite4 = composites[composites.id.str.endswith('-1.0')]
composite5 = composites[composites.id.str.endswith('-1.25')]
composite6 = composites[composites.id.str.endswith('-1.5')]
composite7 = composites[composites.id.str.endswith('-1.75')]
composite8 = composites[composites.id.str.endswith('-2.0')]

#for starforming
starforming1 = starformings[starformings.id.str.endswith('-0.25')]
starforming2 = starformings[starformings.id.str.endswith('-0.5')]
starforming3 = starformings[starformings.id.str.endswith('-0.75')]
starforming4 = starformings[starformings.id.str.endswith('-1.0')]
starforming5 = starformings[starformings.id.str.endswith('-1.25')]
starforming6 = starformings[starformings.id.str.endswith('-1.5')]
starforming7 = starformings[starformings.id.str.endswith('-1.75')]
starforming8 = starformings[starformings.id.str.endswith('-2.0')]

# %% areas

reff = [0.25,0.50,0.75,1.00,1.25,1.50,1.75]

shell1area = shell1["kpc_area"]
shell2area = shell2["kpc_area"]
shell3area = shell3["kpc_area"]
shell4area = shell4["kpc_area"]
shell5area = shell5["kpc_area"]
shell6area = shell6["kpc_area"]
shell7area = shell7["kpc_area"]

control1area = control1["kpc_area"]
control2area = control2["kpc_area"]
control3area = control3["kpc_area"]
control4area = control4["kpc_area"]
control5area = control5["kpc_area"]
control6area = control6["kpc_area"]
control7area = control7["kpc_area"]

#composite area
composite1area = composite1["kpc_area"]
composite2area = composite2["kpc_area"]
composite3area = composite3["kpc_area"]
composite4area = composite4["kpc_area"]
composite5area = composite5["kpc_area"]
composite6area = composite6["kpc_area"]
composite7area = composite7["kpc_area"]

#starforming area
starforming1area = starforming1["kpc_area"]
starforming2area = starforming2["kpc_area"]
starforming3area = starforming3["kpc_area"]
starforming4area = starforming4["kpc_area"]
starforming5area = starforming5["kpc_area"]
starforming6area = starforming6["kpc_area"]
starforming7area = starforming7["kpc_area"]

# %% mass and sfr

#xray masses
shell1mass = np.log10(shell1['bayes.stellar.m_star']/shell1area).dropna()
shell2mass = np.log10(shell2['bayes.stellar.m_star']/shell2area).dropna()
shell3mass = np.log10(shell3['bayes.stellar.m_star']/shell3area).dropna()
shell4mass = np.log10(shell4['bayes.stellar.m_star']/shell4area).dropna()
shell5mass = np.log10(shell5['bayes.stellar.m_star']/shell5area).dropna()
shell6mass = np.log10(shell6['bayes.stellar.m_star']/shell6area).dropna()
shell7mass = np.log10(shell7['bayes.stellar.m_star']/shell7area).dropna()
#shell8mass = np.log10(shell8['bayes.stellar.m_star'])[~np.isnan(shell8['bayes.stellar.m_star'])]/areas[7]

#control masses
control1mass = np.log10(control1['bayes.stellar.m_star']/control1area).dropna()
control2mass = np.log10(control2['bayes.stellar.m_star']/control2area).dropna()
#control2mass = np.log10(control2['bayes.stellar.m_star']/control2area)[control2['bayes.stellar.m_star']>0]/areas[1]
control3mass = np.log10(control3['bayes.stellar.m_star']/control3area).dropna()
control4mass = np.log10(control4['bayes.stellar.m_star']/control4area).dropna()
control5mass = np.log10(control5['bayes.stellar.m_star']/control5area).dropna()
control6mass = np.log10(control6['bayes.stellar.m_star']/control6area).dropna()
control7mass = np.log10(control7['bayes.stellar.m_star']/control7area).dropna()
#control8mass = np.log10(control8[control8['bayes.stellar.m_star']>0]['bayes.stellar.m_star'])

#composite masses
composite1mass = np.log10(composite1['bayes.stellar.m_star']/composite1area).dropna()
composite2mass = np.log10(composite2['bayes.stellar.m_star']/composite2area).dropna()
composite3mass = np.log10(composite3['bayes.stellar.m_star']/composite3area).dropna()
composite4mass = np.log10(composite4['bayes.stellar.m_star']/composite4area).dropna()
composite5mass = np.log10(composite5['bayes.stellar.m_star']/composite5area).dropna()
composite6mass = np.log10(composite6['bayes.stellar.m_star']/composite6area).dropna()
composite7mass = np.log10(composite7['bayes.stellar.m_star']/composite7area).dropna()

#starforming masses
starforming1mass = np.log10(starforming1['bayes.stellar.m_star']/starforming1area).dropna()
starforming2mass = np.log10(starforming2['bayes.stellar.m_star']/starforming2area).dropna()
starforming3mass = np.log10(starforming3['bayes.stellar.m_star']/starforming3area).dropna()
starforming4mass = np.log10(starforming4['bayes.stellar.m_star']/starforming4area).dropna()
starforming5mass = np.log10(starforming5['bayes.stellar.m_star']/starforming5area).dropna()
starforming6mass = np.log10(starforming6['bayes.stellar.m_star']/starforming6area).dropna()
starforming7mass = np.log10(starforming7['bayes.stellar.m_star']/starforming7area).dropna()

#xray sfr
shell1sfr = np.log10(shell1['bayes.sfh.sfr']/shell1area).dropna()
shell2sfr = np.log10(shell2['bayes.sfh.sfr']/shell2area).dropna()
shell3sfr = np.log10(shell3['bayes.sfh.sfr']/shell3area).dropna()
shell4sfr = np.log10(shell4['bayes.sfh.sfr']/shell4area).dropna()
shell5sfr = np.log10(shell5['bayes.sfh.sfr']/shell5area).dropna()
shell6sfr = np.log10(shell6['bayes.sfh.sfr']/shell6area).dropna()
shell7sfr = np.log10(shell7['bayes.sfh.sfr']/shell7area).dropna()
#shell8sfr = np.log10(shell8['bayes.sfh.sfr'])[~np.isnan(shell8['bayes.stellar.m_star'])]

#control sfr
control1sfr = np.log10(control1['bayes.sfh.sfr']/control1area).dropna()
control2sfr = np.log10(control2['bayes.sfh.sfr']/control2area).dropna()
control3sfr = np.log10(control3['bayes.sfh.sfr']/control3area).dropna()
control4sfr = np.log10(control4['bayes.sfh.sfr']/control4area).dropna()
control5sfr = np.log10(control5['bayes.sfh.sfr']/control5area).dropna()
control6sfr = np.log10(control6['bayes.sfh.sfr']/control6area).dropna()
control7sfr = np.log10(control7['bayes.sfh.sfr']/control7area).dropna()

#composite sfr
composite1sfr = np.log10(composite1['bayes.sfh.sfr']/composite1area).dropna()
composite2sfr = np.log10(composite2['bayes.sfh.sfr']/composite2area).dropna()
composite3sfr = np.log10(composite3['bayes.sfh.sfr']/composite3area).dropna()
composite4sfr = np.log10(composite4['bayes.sfh.sfr']/composite4area).dropna()
composite5sfr = np.log10(composite5['bayes.sfh.sfr']/composite5area).dropna()
composite6sfr = np.log10(composite6['bayes.sfh.sfr']/composite6area).dropna()
composite7sfr = np.log10(composite7['bayes.sfh.sfr']/composite7area).dropna()

#starforming sfr
starforming1sfr = np.log10(starforming1['bayes.sfh.sfr']/starforming1area).dropna()
starforming2sfr = np.log10(starforming2['bayes.sfh.sfr']/starforming2area).dropna()
starforming3sfr = np.log10(starforming3['bayes.sfh.sfr']/starforming3area).dropna()
starforming4sfr = np.log10(starforming4['bayes.sfh.sfr']/starforming4area).dropna()
starforming5sfr = np.log10(starforming5['bayes.sfh.sfr']/starforming5area).dropna()
starforming6sfr = np.log10(starforming6['bayes.sfh.sfr']/starforming6area).dropna()
starforming7sfr = np.log10(starforming7['bayes.sfh.sfr']/starforming7area).dropna()

# %% variables for plotting
#------FOR MASS PLOTS-----
mass_to_plot_xray = [shell1mass, shell2mass, shell3mass, shell4mass, shell5mass, shell6mass, shell7mass]
mass_to_plot_control = [control1mass, control2mass, control3mass, control4mass, control5mass, control6mass, control7mass]
mass_to_plot_composite = [composite1mass, composite2mass, composite3mass, composite4mass, composite5mass, composite6mass, composite7mass]
mass_to_plot_starforming = [starforming1mass, starforming2mass, starforming3mass, starforming4mass, starforming5mass, starforming6mass, starforming7mass]


mass_average_xray = [np.nanmedian(i) for i in(mass_to_plot_xray)]
mass_average_control = [np.nanmedian(i) for i in(mass_to_plot_control)]
mass_average_composite = [np.nanmedian(i) for i in(mass_to_plot_composite)]
mass_average_starforming = [np.nanmedian(i) for i in(mass_to_plot_starforming)]

### errors on mass for errorbar plots
mass_xray_err = [sem(ii) for ii in mass_to_plot_xray]
mass_control_err = [sem(jj) for jj in mass_to_plot_control]
mass_composite_err = [sem(jj) for jj in mass_to_plot_composite]
mass_starforming_err = [sem(jj) for jj in mass_to_plot_starforming]


# %% violin plots mass

### half violin plots

#%matplotlib inline

fig, ax = plt.subplots(figsize=(8, 6))
data1 = np.array(mass_to_plot_xray)
data2 = np.array(mass_to_plot_control)

v1 = ax.violinplot(data1, points=100, positions=reff,
               showmeans=False, showextrema=False, showmedians=False, widths=0.25)
for b in v1['bodies']:
    # get the center
    m = np.mean(b.get_paths()[0].vertices[:, 0])
    # modify the paths to not go further right than the center
    b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], -np.inf, m)
    b.set_color('#30a2da')

v2 = ax.violinplot(data2, points=100, positions=reff, 
               showmeans=False, showextrema=False, showmedians=False, widths=0.25)

for b in v2['bodies']:
    # get the center
    m = np.mean(b.get_paths()[0].vertices[:, 0])
    # modify the paths to not go further left than the center
    b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], m, np.inf)
    b.set_color('#fc4f30')


ax.plot(reff,mass_average_xray, color='#30a2da', linestyle='--')
ax.fill_between(reff,np.subtract(mass_average_xray,mass_xray_err),np.add(mass_average_xray,mass_xray_err), color='#30a2da', alpha =0.4)

ax.plot(reff,mass_average_control, color='#fc4f30', linestyle ='--')
ax.fill_between(reff,np.subtract(mass_average_control,mass_control_err),np.add(mass_average_control,mass_control_err), color='#fc4f30', alpha =0.4)


comass = ax.plot(reff,mass_average_composite, color='#4B0092', linestyle ='solid')
ax.fill_between(reff,np.subtract(mass_average_composite,mass_composite_err),np.add(mass_average_composite,mass_composite_err), color='#4B0092', 
                label='Composites', alpha =0.4)

sfmass = ax.plot(reff,mass_average_starforming, color='black', linestyle =':')
ax.fill_between(reff,np.subtract(mass_average_starforming,mass_starforming_err),np.add(mass_average_starforming,mass_starforming_err), color='black', 
                label='Star Forming', alpha =0.4)


ax.legend([v1['bodies'][0],v2['bodies'][0],comass[0],sfmass[0]],['X-Ray', 'Control', 'Composites', 'Star Forming'])
plt.xlabel('R/Re')
plt.ylabel("Stellar Mass Density log[M$_\u2609$kpc$^{_-2}$]")
plt.xticks(reff)
plt.show()


# %% sfr plots variables
#------FOR SFR PLOTS-----
sfr_to_plot_xray = [shell1sfr, shell2sfr, shell3sfr, shell4sfr, shell5sfr, shell6sfr, shell7sfr]
sfr_to_plot_control = [control1sfr, control2sfr, control3sfr, control4sfr, control5sfr, control6sfr, control7sfr]
sfr_to_plot_composite = [composite1sfr, composite2sfr, composite3sfr, composite4sfr, composite5sfr, composite6sfr, composite7sfr]
sfr_to_plot_starforming = [starforming1sfr, starforming2sfr, starforming3sfr, starforming4sfr, starforming5sfr, starforming6sfr, starforming7sfr]

sfr_average_xray = [np.nanmedian(i) for i in(sfr_to_plot_xray)]
sfr_average_control = [np.nanmedian(i) for i in(sfr_to_plot_control)]
sfr_average_composite = [np.nanmedian(i) for i in(sfr_to_plot_composite)]
sfr_average_starforming = [np.nanmedian(i) for i in(sfr_to_plot_starforming)]

### errors on sfr for errorbar plots
sfr_xray_err = [sem(ii) for ii in sfr_to_plot_xray]
sfr_control_err = [sem(jj) for jj in sfr_to_plot_control]
sfr_composite_err = [sem(jj) for jj in sfr_to_plot_composite]
sfr_starforming_err = [sem(jj) for jj in sfr_to_plot_starforming]


# %% violin plots sfr
### half violin plots

fig, ax = plt.subplots(figsize=(8, 6))
data1 = np.array(sfr_to_plot_xray)
data2 = np.array(sfr_to_plot_control)

v1 = ax.violinplot(data1, points=100, positions=reff,
               showmeans=False, showextrema=False, showmedians=False, widths=0.25)
for b in v1['bodies']:
    # get the center
    m = np.mean(b.get_paths()[0].vertices[:, 0])
    # modify the paths to not go further right than the center
    b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], -np.inf, m)
    b.set_color('#30a2da')

v2 = ax.violinplot(data2, points=100, positions=reff, 
               showmeans=False, showextrema=False, showmedians=False, widths=0.25)

for b in v2['bodies']:
    # get the center
    m = np.mean(b.get_paths()[0].vertices[:, 0])
    # modify the paths to not go further left than the center
    b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], m, np.inf)
    b.set_color('#fc4f30')

ax.plot(reff,sfr_average_xray, color='#30a2da', linestyle='--')
ax.fill_between(reff,np.subtract(sfr_average_xray,sfr_xray_err),np.add(sfr_average_xray,sfr_xray_err), color='#30a2da', alpha =0.4)

ax.plot(reff,sfr_average_control, color='#fc4f30', linestyle ='--')
ax.fill_between(reff,np.subtract(sfr_average_control,sfr_control_err),np.add(sfr_average_control,sfr_control_err), color='#fc4f30', alpha =0.4)

cosfr = ax.plot(reff,sfr_average_composite, color='#4B0092', linestyle ='solid')
ax.fill_between(reff,np.subtract(sfr_average_composite,sfr_composite_err),np.add(sfr_average_composite,sfr_composite_err), color='#4B0092', 
                label='Composites', alpha =0.4)

sfsfr = ax.plot(reff,sfr_average_starforming, color='black', linestyle =':')
ax.fill_between(reff,np.subtract(sfr_average_starforming,sfr_starforming_err),np.add(sfr_average_starforming,sfr_starforming_err), color='black', 
                label='Star Forming', alpha =0.4)

ax.legend([v1['bodies'][0],v2['bodies'][0],cosfr[0],sfsfr[0]],['X-Ray', 'Control', 'Composites', 'Star Forming'])
#ax.violinplot(shell1xraymass,showmedians=True, positions=[0.25], widths = 0.15)

#plt.ylim(-3,ax.get_ylim()[1])
plt.xlabel('R/Re')
plt.ylabel("SFR Density log[M$_\u2609$yr$^{_-1}$kpc$^{_-2}$]")
plt.xticks(reff)
plt.show()

# %% ssfr
#xray SSFR
shell1ssfr = np.log10(shell1['bayes.sfh.sfr']/shell1['bayes.stellar.m_star']).dropna()
shell2ssfr = np.log10(shell2['bayes.sfh.sfr']/shell2['bayes.stellar.m_star']).dropna()
shell3ssfr = np.log10(shell3['bayes.sfh.sfr']/shell3['bayes.stellar.m_star']).dropna()
shell4ssfr = np.log10(shell4['bayes.sfh.sfr']/shell4['bayes.stellar.m_star']).dropna()
shell5ssfr = np.log10(shell5['bayes.sfh.sfr']/shell5['bayes.stellar.m_star']).dropna()
shell6ssfr = np.log10(shell6['bayes.sfh.sfr']/shell6['bayes.stellar.m_star']).dropna()
shell7ssfr = np.log10(shell7['bayes.sfh.sfr']/shell7['bayes.stellar.m_star']).dropna()

shell1xrayssfr = np.log10(shells1_xray['bayes.sfh.sfr']/shells1_xray['bayes.stellar.m_star']).dropna()
#shell1xrayssfr = shell1xrayssfr[shell1xrayssfr>-12]


#control SFR
control1ssfr = np.log10(control1['bayes.sfh.sfr']/control1['bayes.stellar.m_star']).dropna()
control2ssfr = np.log10(control2['bayes.sfh.sfr']/control2['bayes.stellar.m_star']).dropna()
control3ssfr = np.log10(control3['bayes.sfh.sfr']/control3['bayes.stellar.m_star']).dropna()
control4ssfr = np.log10(control4['bayes.sfh.sfr']/control4['bayes.stellar.m_star']).dropna()
control5ssfr = np.log10(control5['bayes.sfh.sfr']/control5['bayes.stellar.m_star']).dropna()
control6ssfr = np.log10(control6['bayes.sfh.sfr']/control6['bayes.stellar.m_star']).dropna()
control7ssfr = np.log10(control7['bayes.sfh.sfr']/control7['bayes.stellar.m_star']).dropna()

#composite SSFR
composite1ssfr = np.log10(composite1['bayes.sfh.sfr']/composite1['bayes.stellar.m_star']).dropna()
composite2ssfr = np.log10(composite2['bayes.sfh.sfr']/composite2['bayes.stellar.m_star']).dropna()
composite3ssfr = np.log10(composite3['bayes.sfh.sfr']/composite3['bayes.stellar.m_star']).dropna()
composite4ssfr = np.log10(composite4['bayes.sfh.sfr']/composite4['bayes.stellar.m_star']).dropna()
composite5ssfr = np.log10(composite5['bayes.sfh.sfr']/composite5['bayes.stellar.m_star']).dropna()
composite6ssfr = np.log10(composite6['bayes.sfh.sfr']/composite6['bayes.stellar.m_star']).dropna()
composite7ssfr = np.log10(composite7['bayes.sfh.sfr']/composite7['bayes.stellar.m_star']).dropna()

#starforming SSFR
starforming1ssfr = np.log10(starforming1['bayes.sfh.sfr']/starforming1['bayes.stellar.m_star']).dropna()
starforming2ssfr = np.log10(starforming2['bayes.sfh.sfr']/starforming2['bayes.stellar.m_star']).dropna()
starforming3ssfr = np.log10(starforming3['bayes.sfh.sfr']/starforming3['bayes.stellar.m_star']).dropna()
starforming4ssfr = np.log10(starforming4['bayes.sfh.sfr']/starforming4['bayes.stellar.m_star']).dropna()
starforming5ssfr = np.log10(starforming5['bayes.sfh.sfr']/starforming5['bayes.stellar.m_star']).dropna()
starforming6ssfr = np.log10(starforming6['bayes.sfh.sfr']/starforming6['bayes.stellar.m_star']).dropna()
starforming7ssfr = np.log10(starforming7['bayes.sfh.sfr']/starforming7['bayes.stellar.m_star']).dropna()



#----SSFR PLOTS---- 
ssfr_to_plot_xray = [shell1ssfr, shell2ssfr, shell3ssfr, shell4ssfr, shell5ssfr, shell6ssfr, shell7ssfr]
ssfr_to_plot_control = [control1ssfr, control2ssfr, control3ssfr, control4ssfr, control5ssfr, control6ssfr, control7ssfr]
ssfr_to_plot_composite = [composite1ssfr, composite2ssfr, composite3ssfr, composite4ssfr, composite5ssfr, composite6ssfr, composite7ssfr]
ssfr_to_plot_starforming = [starforming1ssfr, starforming2ssfr, starforming3ssfr, starforming4ssfr, starforming5ssfr, starforming6ssfr, starforming7ssfr]

ssfr_average_xray = [np.nanmedian(i) for i in(ssfr_to_plot_xray)]
ssfr_average_control = [np.nanmedian(i) for i in(ssfr_to_plot_control)]
ssfr_average_composite = [np.nanmedian(i) for i in(ssfr_to_plot_composite)]
ssfr_average_starforming = [np.nanmedian(i) for i in(ssfr_to_plot_starforming)]

# %% ssfr errors
### errors on SSFR for errorbar plots
ssfr_xray_err = [sem(ii) for ii in ssfr_to_plot_xray]
ssfr_control_err = [sem(jj) for jj in ssfr_to_plot_control]
ssfr_composite_err = [sem(jj) for jj in ssfr_to_plot_composite]
ssfr_starforming_err = [sem(jj) for jj in ssfr_to_plot_starforming]


# %%
### half violin plots

fig, ax = plt.subplots(figsize=(8, 6))
data1 = np.array(ssfr_to_plot_xray)
data2 = np.array(ssfr_to_plot_control)

v1 = ax.violinplot(data1, points=100, positions=reff,
               showmeans=False, showextrema=False, showmedians=False, widths=0.25)
for b in v1['bodies']:
    # get the center
    m = np.mean(b.get_paths()[0].vertices[:, 0])
    # modify the paths to not go further right than the center
    b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], -np.inf, m)
    b.set_color('#30a2da')

v2 = ax.violinplot(data2, points=100, positions=reff, 
               showmeans=False, showextrema=False, showmedians=False, widths=0.25)

for b in v2['bodies']:
    # get the center
    m = np.mean(b.get_paths()[0].vertices[:, 0])
    # modify the paths to not go further left than the center
    b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], m, np.inf)
    b.set_color('#fc4f30')

ax.plot(reff,ssfr_average_xray, color='#30a2da', linestyle='--')
ax.fill_between(reff,np.subtract(ssfr_average_xray,ssfr_xray_err),np.add(ssfr_average_xray,ssfr_xray_err), color='#30a2da', alpha =0.4)

ax.plot(reff,ssfr_average_control, color='#fc4f30', linestyle='--')
ax.fill_between(reff,np.subtract(ssfr_average_control,ssfr_control_err),np.add(ssfr_average_control,ssfr_control_err), color='#fc4f30', alpha=0.4)

cosfr = ax.plot(reff,ssfr_average_composite, color='#4B0092', linestyle ='solid')
ax.fill_between(reff,np.subtract(ssfr_average_composite,ssfr_composite_err),np.add(ssfr_average_composite,ssfr_composite_err), color='#4B0092', 
                label='Composites', alpha =0.4)

sfsfr = ax.plot(reff,ssfr_average_starforming, color='black', linestyle =':')
ax.fill_between(reff,np.subtract(ssfr_average_starforming,ssfr_starforming_err),np.add(ssfr_average_starforming,ssfr_starforming_err), color='black', 
                label='Star Forming', alpha =0.4)

ax.legend([v1['bodies'][0],v2['bodies'][0],cosfr[0],sfsfr[0]],['X-Ray', 'Control', 'Composites', 'Star Forming'])

plt.xlabel('R/Re')
plt.ylabel("sSFR log[yr$^{_-1}$]")
plt.xticks(reff)
plt.show()


