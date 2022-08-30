import numpy as np

def lum(flux, distance):
    return 4. * np.pi * (distance * 1000 * 3.0857 * 10 ** 18) ** 2 * 10 ** flux


def spec_plotting(Plot, AllData, product_dir, part, rebin, rescale_params, separate, visible = 0):
    Plot.commands = ()
    Plot.device = f'{product_dir}xspec_{part}.ps/cps'
    Plot.xAxis = "keV"
    Plot.xLog = True
    Plot.yLog = True
    Plot.addCommand('time off')
    Plot.addCommand('wi 1')
    if rebin:
    	Plot.setRebin(rescale_params[0], rescale_params[1])
    	Plot.addCommand(f'rescale y {rescale_params[2]} {rescale_params[3]}')
    Plot.addCommand('la t ')
    Plot.addCommand(r'cs 1.4')
    Plot.addCommand(r'font RO')
    Plot.addCommand(r'la y Counts s\u-1 \dkeV\u-1')
    Plot.addCommand(r'la pos y 2.5')
    if not separate:
        for xi in range(1, AllData.nGroups + 1):
            Plot.addCommand(f'lw 5 on {2 * xi}') #change linewidth of model prediction
            #Plot.addCommand(f'Marker {5 + xi} on {2 * xi - 1}')
            #Plot.addCommand(f'lw 5 on {2 * xi + AllData.nGroups}')
    else:
        for xi in range(1, AllData.nGroups + 1):
            if xi == visible:
                Plot.addCommand(f'lw 5 on {2 * xi}') #change linewidth of model prediction
                Plot.addCommand(f'lw 1 on {2 * xi - 1}')
                Plot.addCommand(f'color 1 on {2 * xi}')
                Plot.addCommand(f'color 1 on {2 * xi - 1}')
                #Plot.addCommand(f'Marker {5 + xi} on {2 * xi - 1}')
                #Plot.addCommand(f'lw 5 on {2 * xi + AllData.nGroups}')
            else:
                #continue
                Plot.addCommand(f'color off {2 * xi}') #change linewidth of model prediction
                Plot.addCommand(f'color off {2 * xi - 1}') #change linewidth of model prediction
    #Plot.addCommand(r'lw 2 on 2') #change linewidth only for the model/data (not sure if 2 corresponds to model or data, need to try)
    #Plot.addCommand(r'la x Energy (keV)')
    Plot.addCommand('wi 2')
    Plot.addCommand(r'la y \gx')
    #Plot.addCommand(r'LAB 2 COL 1 lw 1 LIN 0 100 JUS Lef') #to change color of 0 line in delchi to black
    if not separate:
        for xi in range(1, AllData.nGroups + 1):
            #Plot.addCommand(f'lw 3 on {2 * AllData.nGroups + 2 * xi - 1}')
            #Plot.addCommand(f'Marker {5 + xi} on {2 * AllData.nGroups + 2 * xi - 1}')
            continue
    else:
        off = ''
        for xi in range(1, AllData.nGroups + 1):
            if xi != visible:
                if off == '':
                    off += f'{2 * AllData.nGroups + 2 * xi - 1}'
                else:
                    off += f',{2 * AllData.nGroups + 2 * xi - 1}'
            else:
                Plot.addCommand(f'color 1 on {2 * AllData.nGroups + 2 * xi - 1}')
    Plot.addCommand(r'la x Energy (keV)')
    if rebin:
    	Plot.addCommand(f'rescale y {rescale_params[4]} {rescale_params[5]}')
    Plot.addCommand(f'wenv xspec_{part}')
    if separate:
        Plot.addCommand(f'color off {off}')
    Plot.device = f'{product_dir}xspec_{part}.ps/cps'
    Plot("ldata", "delchi")

def spec_model(Xset, AllModels, AllData, Model, Fit, Plot, product_dir, table, Z, distance, skip_varabs, epoch, absorption, visible, rebin, rescale_params):
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    """
    Created on Wed Dec  8 10:11:33 2021

    @author: kald
    """

    Xset.abund = 'wilm'

    try:
        Z
    except NameError:
        Z = 0.49 #factor to solar metallicity for elements heavier than He; ~0.49 for LMC, 0.25 for SMC

    try:
        skip_varabs
    except NameError:
        skip_varabs = False

    write_table = True
    try:
        epoch
    except NameError:
        write_table = False

    try:
        absorption
    except NameError:
        absorption = 6 * 10 ** (-2)
        
    separate = True
    try:
        visible
    except NameError:
        separate = False

    #def plotting(part):
    #    Plot.device(f'{product_dir}xspec_{part}.ps/cps')
    #    Plot.xAxis("keV")
    #    Plot.xLog = True
    #    Plot.yLog = True
    #    Plot("ldata", "delchi")
    #    #Plot("data", "resid")
    #    return

    #______________________________________________________________________________
    #Part 1: Most general spectral fit (only galactic absorption + powerlaw)
    Xset.openLog(f'{product_dir}xspec_part1.log')
    AllData.show()
    Model("tbabs*powerlaw")
    AllModels(1).TBabs.nH.values = absorption #galactic absorption in the direction of MC in units of 10^22cm^-2
    AllModels(1).TBabs.nH.frozen = True
    if AllData.nGroups > 1:
        for i in range(2, AllData.nGroups + 1):
            AllModels(i).powerlaw.norm.untie()


    Fit.query = "yes"
    Fit.statMethod = "cstat"
    Fit.renorm()
    Fit.perform()
    print("Uncertainty Powerlaw Index")
    Fit.error("2")
    Xset.closeLog()

    model_tbabs = {}
    model_tbabs['PhoIndex'] = AllModels(1).powerlaw.PhoIndex.values #all equal since they are tied
    for i in range(1, AllData.nGroups + 1):
        model_tbabs[f'm{i}_norm'] = AllModels(i).powerlaw.norm.values

    part = "part1"
    spec_plotting(Plot, AllData, product_dir, part, rebin, rescale_params, separate, visible)
    #plotting("part1")


    #______________________________________________________________________________
    #Part 2: Include possible variable absorption from Magellanic Clouds
    if not skip_varabs:
        Xset.openLog(f'{product_dir}xspec_part2.log')
        AllData.show()
        Model("tbabs*tbvarabs*powerlaw")
        AllModels(1).TBabs.nH.values = absorption #galactic absorption in the direction of MC in units of 10^22cm^-2
        AllModels(1).TBabs.nH.frozen = True
        AllModels(1).powerlaw.PhoIndex.values = model_tbabs['PhoIndex']
        for i in range(1, AllData.nGroups + 1):
            AllModels(i).TBvarabs.C.values = Z
            AllModels(i).TBvarabs.N.values = Z
            AllModels(i).TBvarabs.O.values = Z
            AllModels(i).TBvarabs.Ne.values = Z
            AllModels(i).TBvarabs.Na.values = Z
            AllModels(i).TBvarabs.Mg.values = Z
            AllModels(i).TBvarabs.Al.values = Z
            AllModels(i).TBvarabs.Si.values = Z
            AllModels(i).TBvarabs.S.values = Z
            AllModels(i).TBvarabs.Cl.values = Z
            AllModels(i).TBvarabs.Ar.values = Z
            AllModels(i).TBvarabs.Ca.values = Z
            AllModels(i).TBvarabs.Cr.values = Z
            AllModels(i).TBvarabs.Fe.values = Z
            AllModels(i).TBvarabs.Co.values = Z
            AllModels(i).TBvarabs.Ni.values = Z
        if AllData.nGroups > 1:
            for i in range(2, AllData.nGroups + 1):
                AllModels(i).powerlaw.norm.untie()
        for i in range(1, AllData.nGroups + 1):
            AllModels(i).powerlaw.norm.values = model_tbabs[f'm{i}_norm']

        #Fit.nIterations = 100 #in combination with Fit.query = "no"
        Fit.query = "yes"   
        Fit.statMethod = "cstat"
        Fit.renorm()
        Fit.perform()
        print("Uncertainty variable absorption (LMC)")
        Fit.error("2")
        print("Uncertainty Powerlaw Index")
        Fit.error("44")
        Xset.closeLog()

        model_tbvarabs = {}
        model_tbvarabs['PhoIndex'] = AllModels(1).powerlaw.PhoIndex.values #all equal since they are tied
        model_tbvarabs['var_nH'] = AllModels(1).TBvarabs.nH.values
        for i in range(1, AllData.nGroups + 1):
            model_tbvarabs[f'm{i}_norm'] = AllModels(i).powerlaw.norm.values

        part = "part2"
        spec_plotting(Plot, AllData, product_dir, part, rebin, rescale_params, separate, visible)
        #plotting("part2")

    #______________________________________________________________________________
    #Part 3: Flux determination from fitting (absorbed and unabsorbed with 
    #without variable absorption)

    ########################################
    #Part 3.1: absorbed Flux without varabs#
    ########################################
    Xset.openLog(f'{product_dir}xspec_part3_1.log')
    AllData.show()
    Model("cflux*tbabs*powerlaw")
    AllModels(1).TBabs.nH.values = absorption #galactic absorption in the direction of MC in units of 10^22cm^-2
    AllModels(1).TBabs.nH.frozen = True
    AllModels(1).powerlaw.PhoIndex.values = model_tbabs['PhoIndex']
    AllModels(1).powerlaw.norm.values = model_tbabs[f'm1_norm']
    AllModels(1).powerlaw.norm.frozen = True
    if AllData.nGroups > 1:
        for i in range(2, AllData.nGroups + 1):
            AllModels(i).powerlaw.norm.untie()
            AllModels(i).powerlaw.norm.values = model_tbabs[f'm{i}_norm']
            AllModels(i).powerlaw.norm.frozen = True
            AllModels(i).cflux.lg10Flux.untie()
    AllModels(1).cflux.Emin.values = 0.2
    AllModels(1).cflux.Emax.values = 8.0

    #Fit.nIterations = 100 #in combination with Fit.query = "no"
    Fit.query = "yes"
    Fit.statMethod = "cstat"
    #Fit.renorm()
    Fit.perform()
    print("Uncertainty Flux")
    Fit.error("3")
    if AllData.nGroups > 1:
            for i in range(2, AllData.nGroups + 1):
                    index = (i - 1) * 6 + 3
                    Fit.error(f"{index}")
    print("Uncertainty Powerlaw Index")
    Fit.error("5")
    Xset.closeLog()

    #collecting table entries
    if write_table:
        powerlaw_index = [AllModels(1).powerlaw.PhoIndex.values[0], AllModels(1).powerlaw.PhoIndex.error[0], AllModels(1).powerlaw.PhoIndex.error[1]]
        flux_absorbed = []
        for i in range(1, AllData.nGroups + 1):
            flux_absorbed.append([AllModels(i).cflux.lg10Flux.values[0], AllModels(i).cflux.lg10Flux.error[0], AllModels(i).cflux.lg10Flux.error[1]])

    part = "part3_1"
    spec_plotting(Plot, AllData, product_dir, part, rebin, rescale_params, separate, visible)
    #plotting("part3.1")

    ##########################################
    #Part 3.2: unabsorbed Flux without varabs#
    ##########################################
    Xset.openLog(f'{product_dir}xspec_part3_2.log')
    AllData.show()
    Model("tbabs*cflux*powerlaw")
    AllModels(1).TBabs.nH.values = absorption #galactic absorption in the direction of MC in units of 10^22cm^-2
    AllModels(1).TBabs.nH.frozen = True
    AllModels(1).powerlaw.PhoIndex.values = model_tbabs['PhoIndex']
    AllModels(1).powerlaw.norm.values = model_tbabs[f'm1_norm']
    AllModels(1).powerlaw.norm.frozen = True
    if AllData.nGroups > 1:
        for i in range(2, AllData.nGroups + 1):
            AllModels(i).powerlaw.norm.untie()
            AllModels(i).powerlaw.norm.values = model_tbabs[f'm{i}_norm']
            AllModels(i).powerlaw.norm.frozen = True
            AllModels(i).cflux.lg10Flux.untie()
    AllModels(1).cflux.Emin.values = 0.2
    AllModels(1).cflux.Emax.values = 8.0

    #Fit.nIterations = 100 #in combination with Fit.query = "no"
    Fit.query = "yes"
    Fit.statMethod = "cstat"
    #Fit.renorm()
    Fit.perform()
    print("Uncertainty Flux")
    Fit.error("4")
    if AllData.nGroups > 1:
            for i in range(2, AllData.nGroups + 1):
                    index = (i - 1) * 6 + 4
                    Fit.error(f"{index}")
    print("Uncertainty Powerlaw Index")
    Fit.error("5")
    Xset.closeLog()

    #collecting table entries
    if write_table:
        flux_unabsorbed = []
        for i in range(1, AllData.nGroups + 1):
            flux_unabsorbed.append([AllModels(i).cflux.lg10Flux.values[0], AllModels(i).cflux.lg10Flux.error[0], AllModels(i).cflux.lg10Flux.error[1]])

    part = "part3_2"
    spec_plotting(Plot, AllData, product_dir, part, rebin, rescale_params, separate, visible)
    #plotting("part3.2")

    #writing table
    if write_table:
        for i in range(AllData.nGroups):
            table.write(f'{epoch} & {i + 1} & {round(powerlaw_index[0], 3)}$^{{+{round(powerlaw_index[2] - powerlaw_index[0], 3)}}}_{{-{round(powerlaw_index[0] - powerlaw_index[1], 3)}}}$ & -- & {10 ** flux_absorbed[i][0]}$^{{+{10 ** flux_absorbed[i][2] - 10 ** flux_absorbed[i][0]}}}_{{-{10 ** flux_absorbed[i][0] - 10 ** flux_absorbed[i][1]}}}$ & {lum(flux_unabsorbed[i][0], distance)}$^{{+{lum(flux_unabsorbed[i][2], distance) - lum(flux_unabsorbed[i][0], distance)}}}_{{-{lum(flux_unabsorbed[i][0], distance) - lum(flux_unabsorbed[i][1], distance)}}}$\\\\\n')
            table.write(f'& & & & & \\\\ \n')

    #####################################
    #Part 3.3: absorbed Flux with varabs#
    #####################################
    if not skip_varabs:
        Xset.openLog(f'{product_dir}xspec_part3_3.log')
        AllData.show()
        Model("cflux*tbabs*tbvarabs*powerlaw")
        AllModels(1).TBabs.nH.values = absorption #galactic absorption in the direction of MC in units of 10^22cm^-2
        AllModels(1).TBabs.nH.frozen = True
        AllModels(1).powerlaw.PhoIndex.values = model_tbvarabs['PhoIndex']
        for i in range(1, AllData.nGroups + 1):
            AllModels(i).TBvarabs.C.values = Z
            AllModels(i).TBvarabs.N.values = Z
            AllModels(i).TBvarabs.O.values = Z
            AllModels(i).TBvarabs.Ne.values = Z
            AllModels(i).TBvarabs.Na.values = Z
            AllModels(i).TBvarabs.Mg.values = Z
            AllModels(i).TBvarabs.Al.values = Z
            AllModels(i).TBvarabs.Si.values = Z
            AllModels(i).TBvarabs.S.values = Z
            AllModels(i).TBvarabs.Cl.values = Z
            AllModels(i).TBvarabs.Ar.values = Z
            AllModels(i).TBvarabs.Ca.values = Z
            AllModels(i).TBvarabs.Cr.values = Z
            AllModels(i).TBvarabs.Fe.values = Z
            AllModels(i).TBvarabs.Co.values = Z
            AllModels(i).TBvarabs.Ni.values = Z
        AllModels(1).TBvarabs.nH.values = model_tbvarabs['var_nH']
        AllModels(1).powerlaw.norm.values = model_tbvarabs[f'm1_norm']
        AllModels(1).powerlaw.norm.frozen = True
        if AllData.nGroups > 1:
            for i in range(2, AllData.nGroups + 1):
                AllModels(i).powerlaw.norm.untie()
                AllModels(i).powerlaw.norm.values = model_tbvarabs[f'm{i}_norm']
                AllModels(i).powerlaw.norm.frozen = True
                AllModels(i).cflux.lg10Flux.untie()
        AllModels(1).cflux.Emin.values = 0.2
        AllModels(1).cflux.Emax.values = 8.0


        #Fit.nIterations = 100 #in combination with Fit.query = "no"
        Fit.query = "yes"
        Fit.statMethod = "cstat"
        #Fit.renorm()
        Fit.perform()
        print("Uncertainty Flux")
        Fit.error("3")
        if AllData.nGroups > 1:
                for i in range(2, AllData.nGroups + 1):
                        index = (i - 1) * 48 + 3
                        Fit.error(f"{index}")
        print("Uncertainty Variable Absorption (LMC)")
        Fit.error("5")
        print("Uncertainty Powerlaw Index")
        Fit.error("47")
        Xset.closeLog()
        
        #collecting table entries
        if write_table:
            powerlaw_index = [AllModels(1).powerlaw.PhoIndex.values[0], AllModels(1).powerlaw.PhoIndex.error[0], AllModels(1).powerlaw.PhoIndex.error[1]]
            nH_varab = [AllModels(1).TBvarabs.nH.values[0], AllModels(1).TBvarabs.nH.error[0], AllModels(1).TBvarabs.nH.error[1]]
            flux_absorbed = []
            for i in range(1, AllData.nGroups + 1):
                flux_absorbed.append([AllModels(i).cflux.lg10Flux.values[0], AllModels(i).cflux.lg10Flux.error[0], AllModels(i).cflux.lg10Flux.error[1]])

        part = "part3_3"
        spec_plotting(Plot, AllData, product_dir, part, rebin, rescale_params, separate, visible)
        #plotting("part3.3")

    #######################################
    #Part 3.4: unabsorbed Flux with varabs#
    #######################################
    if not skip_varabs:
        Xset.openLog(f'{product_dir}xspec_part3_4.log')
        AllData.show()
        Model("tbabs*tbvarabs*cflux*powerlaw")
        AllModels(1).TBabs.nH.values = absorption #galactic absorption in the direction of MC in units of 10^22cm^-2
        AllModels(1).TBabs.nH.frozen = True
        AllModels(1).powerlaw.PhoIndex.values = model_tbvarabs['PhoIndex']
        for i in range(1, AllData.nGroups + 1):
            AllModels(i).TBvarabs.C.values = Z
            AllModels(i).TBvarabs.N.values = Z
            AllModels(i).TBvarabs.O.values = Z
            AllModels(i).TBvarabs.Ne.values = Z
            AllModels(i).TBvarabs.Na.values = Z
            AllModels(i).TBvarabs.Mg.values = Z
            AllModels(i).TBvarabs.Al.values = Z
            AllModels(i).TBvarabs.Si.values = Z
            AllModels(i).TBvarabs.S.values = Z
            AllModels(i).TBvarabs.Cl.values = Z
            AllModels(i).TBvarabs.Ar.values = Z
            AllModels(i).TBvarabs.Ca.values = Z
            AllModels(i).TBvarabs.Cr.values = Z
            AllModels(i).TBvarabs.Fe.values = Z
            AllModels(i).TBvarabs.Co.values = Z
            AllModels(i).TBvarabs.Ni.values = Z
        AllModels(1).TBvarabs.nH.values = model_tbvarabs['var_nH']
        AllModels(1).powerlaw.norm.values = model_tbvarabs[f'm1_norm']
        AllModels(1).powerlaw.norm.frozen = True
        if AllData.nGroups > 1:
            for i in range(2, AllData.nGroups + 1):
                AllModels(i).powerlaw.norm.untie()
                AllModels(i).powerlaw.norm.values = model_tbvarabs[f'm{i}_norm']
                AllModels(i).powerlaw.norm.frozen = True
                AllModels(i).cflux.lg10Flux.untie()
        AllModels(1).cflux.Emin.values = 0.2
        AllModels(1).cflux.Emax.values = 8.0


        #Fit.nIterations = 100 #in combination with Fit.query = "no"
        Fit.query = "yes"
        Fit.statMethod = "cstat"
        #Fit.renorm()
        Fit.perform()
        print("Uncertainty Variable Absorption (LMC)")
        Fit.error("2")
        print("Uncertainty Flux")
        Fit.error("46")
        if AllData.nGroups > 1:
                for i in range(2, AllData.nGroups + 1):
                        index = (i - 1) * 48 + 46
                        Fit.error(f"{index}")
        print("Uncertainty Powerlaw Index")
        Fit.error("47")
        Xset.closeLog()

        #collecting table entries
        if write_table:
            flux_unabsorbed = []
            for i in range(1, AllData.nGroups + 1):
                flux_unabsorbed.append([AllModels(i).cflux.lg10Flux.values[0], AllModels(i).cflux.lg10Flux.error[0], AllModels(i).cflux.lg10Flux.error[1]])

        part = "part3_4"
        spec_plotting(Plot, AllData, product_dir, part, rebin, rescale_params, separate, visible)
        #plotting("part3.4")

        #writing table
        if write_table:
            for i in range(AllData.nGroups):
                table.write(f'{epoch} & {i + 1} & {round(powerlaw_index[0], 3)}$^{{+{round(powerlaw_index[2] - powerlaw_index[0], 3)}}}_{{-{round(powerlaw_index[0] - powerlaw_index[1], 3)}}}$ & {100. * nH_varab[0]}$^{{+{100 * (nH_varab[2] - nH_varab[0])}}}_{{-{100 * (nH_varab[0] - nH_varab[1])}}}$ & {10 ** flux_absorbed[i][0]}$^{{+{10 ** flux_absorbed[i][2] - 10 ** flux_absorbed[i][0]}}}_{{-{10 ** flux_absorbed[i][0] - 10 ** flux_absorbed[i][1]}}}$ & {lum(flux_unabsorbed[i][0], distance)}$^{{+{lum(flux_unabsorbed[i][2], distance) - lum(flux_unabsorbed[i][0], distance)}}}_{{-{lum(flux_unabsorbed[i][0], distance) - lum(flux_unabsorbed[i][1], distance)}}}$\\\\\n')
                table.write(f'& & & & & \\\\ \n')
