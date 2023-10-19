def apl(Model, AllModels, bxa, galnh, Z, n):
    # frozen parameters
    transformations = []
    srcmod = Model('tbabs*tbvarabs*pow', modName='srcmod', sourceNum=1)
    srcmod.TBabs.nH.values = [galnh, -1]
    srcmod.TBvarabs.C.values = Z
    srcmod.TBvarabs.N.values = Z
    srcmod.TBvarabs.O.values = Z
    srcmod.TBvarabs.Ne.values = Z
    srcmod.TBvarabs.Na.values = Z
    srcmod.TBvarabs.Mg.values = Z
    srcmod.TBvarabs.Al.values = Z
    srcmod.TBvarabs.Si.values = Z
    srcmod.TBvarabs.S.values = Z
    srcmod.TBvarabs.Cl.values = Z
    srcmod.TBvarabs.Ar.values = Z
    srcmod.TBvarabs.Ca.values = Z
    srcmod.TBvarabs.Cr.values = Z
    srcmod.TBvarabs.Fe.values = Z
    srcmod.TBvarabs.Co.values = Z
    srcmod.TBvarabs.Ni.values = Z

    # fit parameters
    loc_nh = srcmod.TBvarabs.nH
    loc_nh.values = [galnh, 0.1, 1e-5, 1e-5, 10, 10]
    p_loc_nh = bxa.create_jeffreys_prior_for(srcmod, loc_nh)
    transformations.append(p_loc_nh)
    gamma = srcmod.powerlaw.PhoIndex
    gamma.values = [1, 0.1, -2, -2, 4, 4]
    p_gamma = bxa.create_uniform_prior_for(srcmod, gamma)
    transformations.append(p_gamma)
    for groupid in range(n):
        model = AllModels(groupNum=2*groupid+1, modName='srcmod')
        norm = model.powerlaw.norm
        norm.values = [1e-4, 0.1, 1e-20, 1e-20, 1e20, 1e20]
        p_norm = bxa.create_jeffreys_prior_for(model, norm)
        p_norm['name'] = f'log(norm{groupid+1})'
        transformations.append(p_norm)
        model_bkg = AllModels(groupNum=2*groupid+2, modName='srcmod')
        model_bkg.powerlaw.norm.values = [0, -1]  # this needs to be tested
    nH = srcmod.TBabs.nH
    nHs_frozen = [nH]
    nHs_modelled = [loc_nh]
    return transformations, nHs_frozen, nHs_modelled, 'apl'


def apl_simple(Model, AllModels, bxa, galnh, Z, n):
    # frozen parameters
    transformations = []
    srcmod = Model('tbabs*pow', modName='srcmod', sourceNum=1)
    srcmod.TBabs.nH.values = [galnh, -1]

    # fit parameters
    gamma = srcmod.powerlaw.PhoIndex
    gamma.values = [1, 0.1, -2, -2, 4, 4]
    p_gamma = bxa.create_uniform_prior_for(srcmod, gamma)
    transformations.append(p_gamma)
    for groupid in range(n):
        model = AllModels(groupNum=2*groupid+1, modName='srcmod')
        norm = model.powerlaw.norm
        norm.values = [1e-4, 0.1, 1e-20, 1e-20, 1e20, 1e20]
        p_norm = bxa.create_jeffreys_prior_for(model, norm)
        p_norm['name'] = f'log(norm{groupid+1})'
        transformations.append(p_norm)
        model_bkg = AllModels(groupNum=2*groupid+2, modName='srcmod')
        model_bkg.powerlaw.norm.values = [0, -1]  # this needs to be tested
    nH = srcmod.TBabs.nH
    nHs_frozen = [nH]
    nHs_modelled = []
    return transformations, nHs_frozen, nHs_modelled, 'apl_simple'
