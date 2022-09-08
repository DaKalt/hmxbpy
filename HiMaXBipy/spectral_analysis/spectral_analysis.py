import numpy as np


def lum(flux, distance):
    return 4. * np.pi * (distance * 1000 * 3.0857 * 10 ** 18) ** 2 * 10 ** flux


def spec_Plotting(Plot_, AllData_, product_dir, part, rebin, rescale_params,
                  separate, plot_command, title, colors, markers, visible=0):
    Plot_.commands = ()
    Plot_.device = f'{product_dir}/working/xspec_{part}.ps/cps'
    Plot_.xAxis = "keV"
    Plot_.xLog = True
    Plot_.yLog = True
    Plot_.addCommand('time off')
    Plot_.addCommand('wi 1')
    if rebin:
        Plot_.setRebin(rescale_params[0], rescale_params[1])
        Plot_.addCommand(f'rescale y {rescale_params[2]} {rescale_params[3]}')
    Plot_.addCommand(f'la t {title}')
    Plot_.addCommand(r'cs 1.4')
    Plot_.addCommand(r'font RO')
    Plot_.addCommand(r'la y Counts s\u-1 \dkeV\u-1')
    Plot_.addCommand(r'la pos y 2.5')
    if not separate:
        for xi in range(1, AllData_.nGroups + 1):
            # change linewidth of Model prediction
            Plot_.addCommand(f'lw 5 on {2 * xi}')
            #Plot_.addCommand(f'Marker {5 + xi} on {2 * xi - 1}')
            #Plot_.addCommand(f'lw 5 on {2 * xi + AllData_.nGroups}')
            if colors != []:
                Plot_.addCommand(
                    f'color {colors[xi]} on {2 * xi}')
                Plot_.addCommand(
                    f'color {colors[xi]} on {2 * xi - 1}')
            if markers != []:
                Plot_.addCommand(
                    f'Marker {markers[xi]} on {2 * xi - 1}')
    else:
        for xi in range(1, AllData_.nGroups + 1):
            if xi == visible:
                # change linewidth of Model prediction
                Plot_.addCommand(f'lw 5 on {2 * xi}')
                Plot_.addCommand(f'lw 1 on {2 * xi - 1}')
                Plot_.addCommand(f'color 1 on {2 * xi}')
                Plot_.addCommand(f'color 1 on {2 * xi - 1}')
                #Plot_.addCommand(f'Marker {5 + xi} on {2 * xi - 1}')
                #Plot_.addCommand(f'lw 5 on {2 * xi + AllData_.nGroups}')
            else:
                # continue
                # change linewidth of Model prediction
                Plot_.addCommand(f'color off {2 * xi}')
                # change linewidth of Model prediction
                Plot_.addCommand(f'color off {2 * xi - 1}')
    # Plot_.addCommand(r'lw 2 on 2') #change linewidth only for the Model/data (not sure if 2 corresponds to Model or data, need to try)
    #Plot_.addCommand(r'la x Energy (keV)')
    Plot_.addCommand('wi 2')
    Plot_.addCommand(r'la y \gx')
    # Plot_.addCommand(r'LAB 2 COL 1 lw 1 LIN 0 100 JUS Lef') #to change color of 0 line in delchi to black
    if not separate:
        for xi in range(1, AllData_.nGroups + 1):
            #Plot_.addCommand(f'lw 3 on {2 * AllData_.nGroups + 2 * xi - 1}')
            #Plot_.addCommand(f'Marker {5 + xi} on {2 * AllData_.nGroups + 2 * xi - 1}')
            if colors != []:
                Plot_.addCommand(
                    f'color {colors[xi]} on {2 * AllData_.nGroups + 2 * xi - 1}')
            if markers != []:
                Plot_.addCommand(
                    f'Marker {markers[xi]} on {2 * AllData_.nGroups + 2 * xi - 1}')
            continue
    else:
        off = ''
        for xi in range(1, AllData_.nGroups + 1):
            if xi != visible:
                if off == '':
                    off += f'{2 * AllData_.nGroups + 2 * xi - 1}'
                else:
                    off += f',{2 * AllData_.nGroups + 2 * xi - 1}'
            else:
                Plot_.addCommand(
                    f'color 1 on {2 * AllData_.nGroups + 2 * xi - 1}')
    Plot_.addCommand(r'la x Energy (keV)')
    if rebin:
        Plot_.addCommand(f'rescale y {rescale_params[4]} {rescale_params[5]}')
    Plot_.addCommand(f'wenv xspec_{part}')
    if separate:
        Plot_.addCommand(f'color off {off}')
    Plot_.device = f'{product_dir}/working/xspec_{part}.ps/cps'
    Plot_(plot_command[0], plot_command[1])


def spec_model(Xset_, AllModels_, AllData_, Model_, Fit_, Plot_, product_dir,
               table, Z, distance, skip_varabs, epoch, absorption, separate,
               visible, rebin, rescale_params, abund, Emin, Emax,
               varabs_starting_pars, plot_command, title,  colors, markers,
               fit_statistic):
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    """
    Created on Wed Dec  8 10:11:33 2021

    @author: kald
    """

    Xset_.abund = abund

    # ______________________________________________________________________________
    # Part 1: Most general spectral Fit_ (only galactic absorption + powerlaw)
    Xset_.openLog(f'{product_dir}/working/xspec_part1.log')
    AllData_.show()
    Model_("tbabs*powerlaw")
    # galactic absorption in the direction of MC in units of 10^22cm^-2
    AllModels_(1).TBabs.nH.values = absorption
    AllModels_(1).TBabs.nH.frozen = True
    AllModels_(1).powerlaw.PhoIndex.values = varabs_starting_pars[0]
    if AllData_.nGroups > 1:
        for i in range(2, AllData_.nGroups + 1):
            AllModels_(i).powerlaw.norm.untie()
            AllModels_(i).powerlaw.PhoIndex.values = varabs_starting_pars[0]

    Fit_.query = "yes"
    Fit_.statMethod = fit_statistic
    Fit_.renorm()
    Fit_.perform()
    print("Uncertainty Powerlaw Index")
    Fit_.error("2")
    Xset_.closeLog()

    Model_tbabs = {}
    Model_tbabs['PhoIndex'] = AllModels_(
        1).powerlaw.PhoIndex.values  # all equal since they are tied
    for i in range(1, AllData_.nGroups + 1):
        Model_tbabs[f'm{i}_norm'] = AllModels_(i).powerlaw.norm.values

    part = "part1"
    spec_Plotting(Plot_, AllData_, product_dir, part,
                  rebin, rescale_params, separate, plot_command, title,
                  colors, markers, visible)
    # Plot_ting("part1")

    # ______________________________________________________________________________
    # Part 2: Include possible variable absorption from Magellanic Clouds
    if not skip_varabs:
        Xset_.openLog(f'{product_dir}/working/xspec_part2.log')
        AllData_.show()
        Model_("tbabs*tbvarabs*powerlaw")
        # galactic absorption in the direction of MC in units of 10^22cm^-2
        AllModels_(1).TBabs.nH.values = absorption
        AllModels_(1).TBabs.nH.frozen = True
        AllModels_(1).powerlaw.PhoIndex.values = Model_tbabs['PhoIndex']
        AllModels_(1).powerlaw.PhoIndex.values = varabs_starting_pars[0]
        for i in range(1, AllData_.nGroups + 1):
            AllModels_(i).TBvarabs.C.values = Z
            AllModels_(i).TBvarabs.N.values = Z
            AllModels_(i).TBvarabs.O.values = Z
            AllModels_(i).TBvarabs.Ne.values = Z
            AllModels_(i).TBvarabs.Na.values = Z
            AllModels_(i).TBvarabs.Mg.values = Z
            AllModels_(i).TBvarabs.Al.values = Z
            AllModels_(i).TBvarabs.Si.values = Z
            AllModels_(i).TBvarabs.S.values = Z
            AllModels_(i).TBvarabs.Cl.values = Z
            AllModels_(i).TBvarabs.Ar.values = Z
            AllModels_(i).TBvarabs.Ca.values = Z
            AllModels_(i).TBvarabs.Cr.values = Z
            AllModels_(i).TBvarabs.Fe.values = Z
            AllModels_(i).TBvarabs.Co.values = Z
            AllModels_(i).TBvarabs.Ni.values = Z
            AllModels_(i).TBvarabs.nH.values = varabs_starting_pars[1]
            AllModels_(i).powerlaw.PhoIndex.values = varabs_starting_pars[0]
        if AllData_.nGroups > 1:
            for i in range(2, AllData_.nGroups + 1):
                AllModels_(i).powerlaw.norm.untie()
        for i in range(1, AllData_.nGroups + 1):
            AllModels_(i).powerlaw.norm.values = Model_tbabs[f'm{i}_norm']

        # Fit_.nIterations = 100 #in combination with Fit_.query = "no"
        Fit_.query = "yes"
        Fit_.statMethod = "cstat"
        Fit_.renorm()
        Fit_.perform()
        print("Uncertainty variable absorption (LMC)")
        Fit_.error("2")
        print("Uncertainty Powerlaw Index")
        Fit_.error("44")
        Xset_.closeLog()

        Model_tbvarabs = {}
        Model_tbvarabs['PhoIndex'] = AllModels_(
            1).powerlaw.PhoIndex.values  # all equal since they are tied
        Model_tbvarabs['var_nH'] = AllModels_(1).TBvarabs.nH.values
        for i in range(1, AllData_.nGroups + 1):
            Model_tbvarabs[f'm{i}_norm'] = AllModels_(i).powerlaw.norm.values

        part = "part2"
        spec_Plotting(Plot_, AllData_, product_dir, part,
                      rebin, rescale_params, separate, plot_command, title,
                      colors, markers, visible)
        # Plot_ting("part2")

    # ______________________________________________________________________________
    # Part 3: Flux determination from Fit_ting (absorbed and unabsorbed with
    # without variable absorption)

    ########################################
    #Part 3.1: absorbed Flux without varabs#
    ########################################
    Xset_.openLog(f'{product_dir}/working/xspec_part3_1.log')
    AllData_.show()
    Model_("cflux*tbabs*powerlaw")
    # galactic absorption in the direction of MC in units of 10^22cm^-2
    AllModels_(1).TBabs.nH.values = absorption
    AllModels_(1).TBabs.nH.frozen = True
    AllModels_(1).powerlaw.PhoIndex.values = Model_tbabs['PhoIndex']
    AllModels_(1).powerlaw.norm.values = Model_tbabs['m1_norm']
    AllModels_(1).powerlaw.norm.frozen = True
    AllModels_(1).powerlaw.PhoIndex.values = varabs_starting_pars[0]
    if AllData_.nGroups > 1:
        for i in range(2, AllData_.nGroups + 1):
            AllModels_(i).powerlaw.norm.untie()
            AllModels_(i).powerlaw.norm.values = Model_tbabs[f'm{i}_norm']
            AllModels_(i).powerlaw.norm.frozen = True
            AllModels_(i).cflux.lg10Flux.untie()
            AllModels_(i).powerlaw.PhoIndex.values = varabs_starting_pars[0]
    AllModels_(1).cflux.Emin.values = Emin
    AllModels_(1).cflux.Emax.values = Emax

    # Fit_.nIterations = 100 #in combination with Fit_.query = "no"
    Fit_.query = "yes"
    Fit_.statMethod = "cstat"
    # Fit_.renorm()
    Fit_.perform()
    print("Uncertainty Flux")
    Fit_.error("3")
    if AllData_.nGroups > 1:
        for i in range(2, AllData_.nGroups + 1):
            index = (i - 1) * 6 + 3
            Fit_.error(f"{index}")
    print("Uncertainty Powerlaw Index")
    Fit_.error("5")
    Xset_.closeLog()

    # collecting table entries
    powerlaw_index = [AllModels_(1).powerlaw.PhoIndex.values[0], AllModels_(
        1).powerlaw.PhoIndex.error[0], AllModels_(1).powerlaw.PhoIndex.error[1]]
    flux_absorbed = []
    for i in range(1, AllData_.nGroups + 1):
        flux_absorbed.append([AllModels_(i).cflux.lg10Flux.values[0], AllModels_(
            i).cflux.lg10Flux.error[0], AllModels_(i).cflux.lg10Flux.error[1]])

    part = "part3_1"
    spec_Plotting(Plot_, AllData_, product_dir, part,
                  rebin, rescale_params, separate, plot_command, title,
                  colors, markers, visible)
    # Plot_ting("part3.1")

    ##########################################
    #Part 3.2: unabsorbed Flux without varabs#
    ##########################################
    Xset_.openLog(f'{product_dir}/working/xspec_part3_2.log')
    AllData_.show()
    Model_("tbabs*cflux*powerlaw")
    # galactic absorption in the direction of MC in units of 10^22cm^-2
    AllModels_(1).TBabs.nH.values = absorption
    AllModels_(1).TBabs.nH.frozen = True
    AllModels_(1).powerlaw.PhoIndex.values = Model_tbabs['PhoIndex']
    AllModels_(1).powerlaw.norm.values = Model_tbabs['m1_norm']
    AllModels_(1).powerlaw.norm.frozen = True
    AllModels_(1).powerlaw.PhoIndex.values = varabs_starting_pars[0]
    if AllData_.nGroups > 1:
        for i in range(2, AllData_.nGroups + 1):
            AllModels_(i).powerlaw.norm.untie()
            AllModels_(i).powerlaw.norm.values = Model_tbabs[f'm{i}_norm']
            AllModels_(i).powerlaw.norm.frozen = True
            AllModels_(i).cflux.lg10Flux.untie()
            AllModels_(i).powerlaw.PhoIndex.values = varabs_starting_pars[0]
    AllModels_(1).cflux.Emin.values = Emin
    AllModels_(1).cflux.Emax.values = Emax

    # Fit_.nIterations = 100 #in combination with Fit_.query = "no"
    Fit_.query = "yes"
    Fit_.statMethod = "cstat"
    # Fit_.renorm()
    Fit_.perform()
    print("Uncertainty Flux")
    Fit_.error("4")
    if AllData_.nGroups > 1:
        for i in range(2, AllData_.nGroups + 1):
            index = (i - 1) * 6 + 4
            Fit_.error(f"{index}")
    print("Uncertainty Powerlaw Index")
    Fit_.error("5")
    Xset_.closeLog()

    # collecting table entries
    flux_unabsorbed = []
    for i in range(1, AllData_.nGroups + 1):
        flux_unabsorbed.append([AllModels_(i).cflux.lg10Flux.values[0], AllModels_(
            i).cflux.lg10Flux.error[0], AllModels_(i).cflux.lg10Flux.error[1]])

    part = "part3_2"
    spec_Plotting(Plot_, AllData_, product_dir, part,
                  rebin, rescale_params, separate, plot_command, title,
                  colors, markers, visible)
    # Plot_ting("part3.2")

    # writing table
    for i in range(AllData_.nGroups):
        table.write(f'{epoch} & {i + 1} & {round(powerlaw_index[0], 3)}$^{{+{round(powerlaw_index[2] - powerlaw_index[0], 3)}}}_{{-{round(powerlaw_index[0] - powerlaw_index[1], 3)}}}$ & -- & {10 ** flux_absorbed[i][0]}$^{{+{10 ** flux_absorbed[i][2] - 10 ** flux_absorbed[i][0]}}}_{{-{10 ** flux_absorbed[i][0] - 10 ** flux_absorbed[i][1]}}}$ & {lum(flux_unabsorbed[i][0], distance)}$^{{+{lum(flux_unabsorbed[i][2], distance) - lum(flux_unabsorbed[i][0], distance)}}}_{{-{lum(flux_unabsorbed[i][0], distance) - lum(flux_unabsorbed[i][1], distance)}}}$\\\\\n')
        table.write('& & & & & \\\\ \n')

    #####################################
    #Part 3.3: absorbed Flux with varabs#
    #####################################
    if not skip_varabs:
        Xset_.openLog(f'{product_dir}/working/xspec_part3_3.log')
        AllData_.show()
        Model_("cflux*tbabs*tbvarabs*powerlaw")
        # galactic absorption in the direction of MC in units of 10^22cm^-2
        AllModels_(1).TBabs.nH.values = absorption
        AllModels_(1).TBabs.nH.frozen = True
        AllModels_(1).powerlaw.PhoIndex.values = Model_tbvarabs['PhoIndex']
        for i in range(1, AllData_.nGroups + 1):
            AllModels_(i).TBvarabs.C.values = Z
            AllModels_(i).TBvarabs.N.values = Z
            AllModels_(i).TBvarabs.O.values = Z
            AllModels_(i).TBvarabs.Ne.values = Z
            AllModels_(i).TBvarabs.Na.values = Z
            AllModels_(i).TBvarabs.Mg.values = Z
            AllModels_(i).TBvarabs.Al.values = Z
            AllModels_(i).TBvarabs.Si.values = Z
            AllModels_(i).TBvarabs.S.values = Z
            AllModels_(i).TBvarabs.Cl.values = Z
            AllModels_(i).TBvarabs.Ar.values = Z
            AllModels_(i).TBvarabs.Ca.values = Z
            AllModels_(i).TBvarabs.Cr.values = Z
            AllModels_(i).TBvarabs.Fe.values = Z
            AllModels_(i).TBvarabs.Co.values = Z
            AllModels_(i).TBvarabs.Ni.values = Z
            AllModels_(i).TBvarabs.nH.values = varabs_starting_pars[1]
            AllModels_(i).powerlaw.PhoIndex.values = varabs_starting_pars[0]
        AllModels_(1).TBvarabs.nH.values = Model_tbvarabs['var_nH']
        AllModels_(1).powerlaw.norm.values = Model_tbvarabs['m1_norm']
        AllModels_(1).powerlaw.norm.frozen = True
        if AllData_.nGroups > 1:
            for i in range(2, AllData_.nGroups + 1):
                AllModels_(i).powerlaw.norm.untie()
                AllModels_(
                    i).powerlaw.norm.values = Model_tbvarabs[f'm{i}_norm']
                AllModels_(i).powerlaw.norm.frozen = True
                AllModels_(i).cflux.lg10Flux.untie()
        AllModels_(1).cflux.Emin.values = Emin
        AllModels_(1).cflux.Emax.values = Emax

        # Fit_.nIterations = 100 #in combination with Fit_.query = "no"
        Fit_.query = "yes"
        Fit_.statMethod = "cstat"
        # Fit_.renorm()
        Fit_.perform()
        print("Uncertainty Flux")
        Fit_.error("3")
        if AllData_.nGroups > 1:
            for i in range(2, AllData_.nGroups + 1):
                index = (i - 1) * 48 + 3
                Fit_.error(f"{index}")
        print("Uncertainty Variable Absorption (LMC)")
        Fit_.error("5")
        print("Uncertainty Powerlaw Index")
        # Fit_.error("47")
        Xset_.closeLog()

        # collecting table entries
        powerlaw_index = [AllModels_(1).powerlaw.PhoIndex.values[0], AllModels_(
            1).powerlaw.PhoIndex.error[0], AllModels_(1).powerlaw.PhoIndex.error[1]]
        nH_varab = [AllModels_(1).TBvarabs.nH.values[0], AllModels_(
            1).TBvarabs.nH.error[0], AllModels_(1).TBvarabs.nH.error[1]]
        flux_absorbed = []
        for i in range(1, AllData_.nGroups + 1):
            flux_absorbed.append([AllModels_(i).cflux.lg10Flux.values[0], AllModels_(
                i).cflux.lg10Flux.error[0], AllModels_(i).cflux.lg10Flux.error[1]])

        part = "part3_3"
        spec_Plotting(Plot_, AllData_, product_dir, part,
                      rebin, rescale_params, separate, plot_command, title,
                      colors, markers, visible)
        # Plot_ting("part3.3")

    #######################################
    #Part 3.4: unabsorbed Flux with varabs#
    #######################################
    if not skip_varabs:
        Xset_.openLog(f'{product_dir}/working/xspec_part3_4.log')
        AllData_.show()
        Model_("tbabs*tbvarabs*cflux*powerlaw")
        # galactic absorption in the direction of MC in units of 10^22cm^-2
        AllModels_(1).TBabs.nH.values = absorption
        AllModels_(1).TBabs.nH.frozen = True
        AllModels_(1).powerlaw.PhoIndex.values = Model_tbvarabs['PhoIndex']
        for i in range(1, AllData_.nGroups + 1):
            AllModels_(i).TBvarabs.C.values = Z
            AllModels_(i).TBvarabs.N.values = Z
            AllModels_(i).TBvarabs.O.values = Z
            AllModels_(i).TBvarabs.Ne.values = Z
            AllModels_(i).TBvarabs.Na.values = Z
            AllModels_(i).TBvarabs.Mg.values = Z
            AllModels_(i).TBvarabs.Al.values = Z
            AllModels_(i).TBvarabs.Si.values = Z
            AllModels_(i).TBvarabs.S.values = Z
            AllModels_(i).TBvarabs.Cl.values = Z
            AllModels_(i).TBvarabs.Ar.values = Z
            AllModels_(i).TBvarabs.Ca.values = Z
            AllModels_(i).TBvarabs.Cr.values = Z
            AllModels_(i).TBvarabs.Fe.values = Z
            AllModels_(i).TBvarabs.Co.values = Z
            AllModels_(i).TBvarabs.Ni.values = Z
            AllModels_(i).TBvarabs.nH.values = varabs_starting_pars[1]
            AllModels_(i).powerlaw.PhoIndex.values = varabs_starting_pars[0]
        AllModels_(1).TBvarabs.nH.values = Model_tbvarabs['var_nH']
        AllModels_(1).powerlaw.norm.values = Model_tbvarabs['m1_norm']
        AllModels_(1).powerlaw.norm.frozen = True
        if AllData_.nGroups > 1:
            for i in range(2, AllData_.nGroups + 1):
                AllModels_(i).powerlaw.norm.untie()
                AllModels_(
                    i).powerlaw.norm.values = Model_tbvarabs[f'm{i}_norm']
                AllModels_(i).powerlaw.norm.frozen = True
                AllModels_(i).cflux.lg10Flux.untie()
        AllModels_(1).cflux.Emin.values = Emin
        AllModels_(1).cflux.Emax.values = Emax

        # Fit_.nIterations = 100 #in combination with Fit_.query = "no"
        Fit_.query = "yes"
        Fit_.statMethod = "cstat"
        # Fit_.renorm()
        Fit_.perform()
        print("Uncertainty Variable Absorption (LMC)")
        Fit_.error("2")
        print("Uncertainty Flux")
        Fit_.error("46")
        if AllData_.nGroups > 1:
            for i in range(2, AllData_.nGroups + 1):
                index = (i - 1) * 48 + 46
                Fit_.error(f"{index}")
        print("Uncertainty Powerlaw Index")
        # Fit_.error("47")
        Xset_.closeLog()

        # collecting table entries
        flux_unabsorbed = []
        for i in range(1, AllData_.nGroups + 1):
            flux_unabsorbed.append([AllModels_(i).cflux.lg10Flux.values[0], AllModels_(
                i).cflux.lg10Flux.error[0], AllModels_(i).cflux.lg10Flux.error[1]])

        part = "part3_4"
        spec_Plotting(Plot_, AllData_, product_dir, part,
                      rebin, rescale_params, separate, plot_command, title,
                      colors, markers, visible)
        # Plot_ting("part3.4")

        # writing table
        for i in range(AllData_.nGroups):
            table.write(f'{epoch} & {i + 1} & {round(powerlaw_index[0], 3)}$^{{+{round(powerlaw_index[2] - powerlaw_index[0], 3)}}}_{{-{round(powerlaw_index[0] - powerlaw_index[1], 3)}}}$ & {100. * nH_varab[0]}$^{{+{100 * (nH_varab[2] - nH_varab[0])}}}_{{-{100 * (nH_varab[0] - nH_varab[1])}}}$ & {10 ** flux_absorbed[i][0]}$^{{+{10 ** flux_absorbed[i][2] - 10 ** flux_absorbed[i][0]}}}_{{-{10 ** flux_absorbed[i][0] - 10 ** flux_absorbed[i][1]}}}$ & {lum(flux_unabsorbed[i][0], distance)}$^{{+{lum(flux_unabsorbed[i][2], distance) - lum(flux_unabsorbed[i][0], distance)}}}_{{-{lum(flux_unabsorbed[i][0], distance) - lum(flux_unabsorbed[i][1], distance)}}}$\\\\\n')
            table.write('& & & & & \\\\ \n')
