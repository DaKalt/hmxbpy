from HiMaXBipy.spectral_analysis.modded_function import\
    modded_create_uniform_prior_for, modded_create_jeffreys_prior_for

def apl(Model, AllModels, bxa, galnh, Z, n):
    # frozen parameters
    transformations = []
    srcmod = Model('tbabs*tbvarabs*pow', modName='srcmod', sourceNum=1)
    srcmod.TBabs.nH.values = [galnh, -1, 0, 1e-5, 1e2, 1e2]
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
    loc_nh.values = [galnh, 0.01, 0, 1e-5, 1e2, 1e2]
    p_loc_nh = modded_create_jeffreys_prior_for(srcmod, loc_nh)
    transformations.append(p_loc_nh)
    gamma = srcmod.powerlaw.PhoIndex
    gamma.values = [1, 0.01, -2, -2, 4, 4]
    p_gamma = modded_create_uniform_prior_for(srcmod, gamma)
    p_gamma['name'] = '$\\alpha$' #was Gamma before but apparently xspec
    # uses alpha
    transformations.append(p_gamma)
    for groupid in range(n):
        model = AllModels(groupNum=2*groupid+1, modName='srcmod')
        norm = model.powerlaw.norm
        norm.values = [1e-4, 0.01, 0, 1e-8, 1e2, 1e20]
        p_norm = modded_create_jeffreys_prior_for(model, norm)
        p_norm['name'] = 'log(norm$_{PL,%s}$)' % (groupid+1)
        transformations.append(p_norm)
        model_bkg = AllModels(groupNum=2*groupid+2, modName='srcmod')
        model_bkg.powerlaw.norm.values = [0, -1]  # this needs to be tested
    nH = srcmod.TBabs.nH
    nHs_frozen = [nH]
    nHs_modelled = [0]
    return transformations, nHs_frozen, nHs_modelled, 'apl'


def apl_simple(Model, AllModels, bxa, galnh, Z, n):
    # frozen parameters
    transformations = []
    srcmod = Model('tbabs*pow', modName='srcmod', sourceNum=1)
    srcmod.TBabs.nH.values = [galnh, -1, 0, 1e-5, 1e2, 1e2]

    # fit parameters
    gamma = srcmod.powerlaw.PhoIndex
    gamma.values = [1, 0.01, -2, -2, 4, 4]
    p_gamma = modded_create_uniform_prior_for(srcmod, gamma)
    p_gamma['name'] = '$\\alpha$'
    transformations.append(p_gamma)
    for groupid in range(n):
        model = AllModels(groupNum=2*groupid+1, modName='srcmod')
        norm = model.powerlaw.norm
        norm.values = [1e-4, 0.01, 0, 1e-8, 1e2, 1e20]
        p_norm = modded_create_jeffreys_prior_for(model, norm)
        p_norm['name'] = 'log(norm$_{PL,%s}$)' % (groupid+1)
        transformations.append(p_norm)
        model_bkg = AllModels(groupNum=2*groupid+2, modName='srcmod')
        model_bkg.powerlaw.norm.values = [0, -1]  # this needs to be tested
    nH = srcmod.TBabs.nH
    nHs_frozen = [nH]
    nHs_modelled = []
    return transformations, nHs_frozen, nHs_modelled, 'apl_simple'


def abb_simple(Model, AllModels, bxa, galnh, Z, n):
    # frozen parameters
    transformations = []
    srcmod = Model('tbabs*bbodyrad', modName='srcmod', sourceNum=1)
    srcmod.TBabs.nH.values = [galnh, -1, 0, 1e-5, 1e2, 1e2]

    # fit parameters
    kT = srcmod.bbodyrad.kT
    kT.values = [0.1, 0.005, 0.001, 0.001, 10, 10]
    p_kT = modded_create_uniform_prior_for(srcmod, kT)
    transformations.append(p_kT)
    for groupid in range(n):
        model = AllModels(groupNum=2*groupid+1, modName='srcmod')
        norm = model.bbodyrad.norm
        norm.values = [1e-4, 0.01, 0, 1e-3, 1e7, 1e20]
        p_norm = modded_create_jeffreys_prior_for(model, norm)
        p_norm['name'] = 'log(norm$_{BB,%s}$)' % (groupid+1)
        transformations.append(p_norm)
        model_bkg = AllModels(groupNum=2*groupid+2, modName='srcmod')
        model_bkg.bbodyrad.norm.values = [0, -1]  # this needs to be tested
    nH = srcmod.TBabs.nH
    nHs_frozen = [nH]
    nHs_modelled = []
    return transformations, nHs_frozen, nHs_modelled, 'abb_simple'


def abb(Model, AllModels, bxa, galnh, Z, n):
    # frozen parameters
    transformations = []
    srcmod = Model('tbabs*tbvarabs*bbodyrad', modName='srcmod', sourceNum=1)
    srcmod.TBabs.nH.values = [galnh, -1, 0, 1e-5, 1e2, 1e2]
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
    loc_nh.values = [galnh, 0.01, 0, 1e-5, 1e2, 1e2]
    p_loc_nh = modded_create_jeffreys_prior_for(srcmod, loc_nh)
    transformations.append(p_loc_nh)
    kT = srcmod.bbodyrad.kT
    kT.values = [0.1, 0.005, 0.001, 0.001, 10, 10]
    p_kT = modded_create_uniform_prior_for(srcmod, kT)
    transformations.append(p_kT)
    for groupid in range(n):
        model = AllModels(groupNum=2*groupid+1, modName='srcmod')
        norm = model.bbodyrad.norm
        norm.values = [1e-4, 0.01, 0, 1e-3, 1e7, 1e20]
        p_norm = modded_create_jeffreys_prior_for(model, norm)
        p_norm['name'] = 'log(norm$_{BB,%s}$)' % (groupid+1)
        transformations.append(p_norm)
        model_bkg = AllModels(groupNum=2*groupid+2, modName='srcmod')
        model_bkg.bbodyrad.norm.values = [0, -1]  # this needs to be tested
    nH = srcmod.TBabs.nH
    nHs_frozen = [nH]
    nHs_modelled = [0]
    return transformations, nHs_frozen, nHs_modelled, 'abb'

def apl_diskbb_simple(Model, AllModels, bxa, galnh, Z, n):
    # frozen parameters
    transformations = []
    srcmod = Model('tbabs*(pow+diskbb)', modName='srcmod', sourceNum=1)
    srcmod.TBabs.nH.values = [galnh, -1, 0, 1e-5, 1e2, 1e2]

    # fit parameters
    gamma = srcmod.powerlaw.PhoIndex
    gamma.values = [1, 0.01, -2, -2, 4, 4]
    p_gamma = modded_create_uniform_prior_for(srcmod, gamma)
    p_gamma['name'] = '$\\alpha$'
    transformations.append(p_gamma)
    Tin = srcmod.diskbb.Tin
    Tin.values = [0.1, 0.005, 0.001, 0.001, 0.3, 10]
    p_Tin = modded_create_uniform_prior_for(srcmod, Tin)
    p_Tin['name'] = '$T_{in}$'
    transformations.append(p_Tin)
    for groupid in range(n):
        model = AllModels(groupNum=2*groupid+1, modName='srcmod')
        norm = model.powerlaw.norm
        norm.values = [1e-4, 0.01, 0, 1e-8, 1e2, 1e20]
        p_norm = modded_create_jeffreys_prior_for(model, norm)
        p_norm['name'] = 'log(norm$_{PL,%s}$)' % (groupid+1)
        transformations.append(p_norm)
        norm_disk = model.diskbb.norm
        norm_disk.values = [1e-4, 0.01, 0, 1e-6, 1e7, 1e20]
        p_norm_disk = modded_create_jeffreys_prior_for(model, norm_disk)
        p_norm_disk['name'] = 'log(norm$_{disk,%s}$)' % (groupid+1)
        transformations.append(p_norm_disk)
        model_bkg = AllModels(groupNum=2*groupid+2, modName='srcmod')
        model_bkg.powerlaw.norm.values = [0, -1]  # this needs to be tested
        model_bkg.diskbb.norm.values = [0, -1]
    nH = srcmod.TBabs.nH
    nHs_frozen = [nH]
    nHs_modelled = []
    return transformations, nHs_frozen, nHs_modelled, 'apl_diskbb_simple'

def apl_diskbb(Model, AllModels, bxa, galnh, Z, n):
    # frozen parameters
    transformations = []
    srcmod = Model('tbabs*tbvarabs*(pow+diskbb)', modName='srcmod', sourceNum=1)
    srcmod.TBabs.nH.values = [galnh, -1, 0, 1e-5, 1e2, 1e2]
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
    loc_nh.values = [galnh, 0.01, 0, 1e-5, 1e2, 1e2]
    p_loc_nh = modded_create_jeffreys_prior_for(srcmod, loc_nh)
    transformations.append(p_loc_nh)
    gamma = srcmod.powerlaw.PhoIndex
    gamma.values = [1, 0.01, -2, -2, 4, 4]
    p_gamma = modded_create_uniform_prior_for(srcmod, gamma)
    p_gamma['name'] = '$\\alpha$' #was Gamma before but apparently xspec
    # uses alpha
    transformations.append(p_gamma)
    Tin = srcmod.diskbb.Tin
    Tin.values = [0.1, 0.005, 0.001, 0.001, 10, 10]
    p_Tin = modded_create_uniform_prior_for(srcmod, Tin)
    p_Tin['name'] = '$T_{in}$'
    transformations.append(p_Tin)
    for groupid in range(n):
        model = AllModels(groupNum=2*groupid+1, modName='srcmod')
        norm = model.powerlaw.norm
        norm.values = [1e-4, 0.01, 0, 1e-8, 1e2, 1e20]
        p_norm = modded_create_jeffreys_prior_for(model, norm)
        p_norm['name'] = 'log(norm$_{PL,%s}$)' % (groupid+1)
        transformations.append(p_norm)
        norm_disk = model.diskbb.norm
        norm_disk.values = [1e-4, 0.01, 0, 1e-6, 1e7, 1e20]
        p_norm_disk = modded_create_jeffreys_prior_for(model, norm_disk)
        p_norm_disk['name'] = 'log(norm$_{disk,%s}$)' % (groupid+1)
        transformations.append(p_norm_disk)
        model_bkg = AllModels(groupNum=2*groupid+2, modName='srcmod')
        model_bkg.powerlaw.norm.values = [0, -1]  # this needs to be tested
        model_bkg.diskbb.norm.values = [0, -1]  # this needs to be tested
    nH = srcmod.TBabs.nH
    nHs_frozen = [nH]
    nHs_modelled = [0]
    return transformations, nHs_frozen, nHs_modelled, 'apl_diskbb'

def apl_hmxb(Model, AllModels, bxa, galnh, Z, n):
    # frozen parameters
    transformations = []
    srcmod = Model('tbabs*tbvarabs*pow', modName='srcmod', sourceNum=1)
    srcmod.TBabs.nH.values = [galnh, -1, 0, 1e-5, 1e2, 1e2]
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
    loc_nh.values = [galnh, 0.01, 0, 1e-5, 1e2, 1e2]
    p_loc_nh = modded_create_jeffreys_prior_for(srcmod, loc_nh)
    transformations.append(p_loc_nh)
    gamma = srcmod.powerlaw.PhoIndex
    gamma.values = [1, -1]
    for groupid in range(n):
        model = AllModels(groupNum=2*groupid+1, modName='srcmod')
        norm = model.powerlaw.norm
        norm.values = [1e-4, 0.01, 0, 1e-8, 1e2, 1e20]
        p_norm = modded_create_jeffreys_prior_for(model, norm)
        p_norm['name'] = 'log(norm$_{PL,%s}$)' % (groupid+1)
        transformations.append(p_norm)
        model_bkg = AllModels(groupNum=2*groupid+2, modName='srcmod')
        model_bkg.powerlaw.norm.values = [0, -1]  # this needs to be tested
    nH = srcmod.TBabs.nH
    nHs_frozen = [nH]
    nHs_modelled = [0]
    return transformations, nHs_frozen, nHs_modelled, 'apl_hmxb'

def apl_bb(Model, AllModels, bxa, galnh, Z, n):
    # frozen parameters
    transformations = []
    srcmod = Model('tbabs*tbvarabs*(pow+bbodyrad)', modName='srcmod', sourceNum=1)
    srcmod.TBabs.nH.values = [galnh, -1, 0, 1e-5, 1e2, 1e2]
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
    loc_nh.values = [galnh, 0.01, 0, 1e-5, 1e2, 1e2]
    p_loc_nh = modded_create_jeffreys_prior_for(srcmod, loc_nh)
    transformations.append(p_loc_nh)
    gamma = srcmod.powerlaw.PhoIndex
    gamma.values = [1, 0.01, -2, -2, 4, 4]
    p_gamma = modded_create_uniform_prior_for(srcmod, gamma)
    p_gamma['name'] = '$\\alpha$' #was Gamma before but apparently xspec
    # uses alpha
    transformations.append(p_gamma)
    kT = srcmod.bbodyrad.kT
    kT.values = [0.1, 0.005, 0.001, 0.001, 0.3, 10]
    p_kT = modded_create_uniform_prior_for(srcmod, kT)
    transformations.append(p_kT)
    for groupid in range(n):
        model = AllModels(groupNum=2*groupid+1, modName='srcmod')
        norm = model.powerlaw.norm
        norm.values = [1e-4, 0.01, 0, 1e-8, 1e2, 1e20]
        p_norm = modded_create_jeffreys_prior_for(model, norm)
        p_norm['name'] = 'log(norm$_{PL,%s}$)' % (groupid+1)
        transformations.append(p_norm)
        norm_bb = model.bbodyrad.norm
        norm_bb.values = [1e-4, 0.01, 0, 1e-3, 1e7, 1e20]
        p_norm_bb = modded_create_jeffreys_prior_for(model, norm_bb)
        p_norm_bb['name'] = 'log(norm$_{BB,%s}$)' % (groupid+1)
        transformations.append(p_norm_bb)
        model_bkg = AllModels(groupNum=2*groupid+2, modName='srcmod')
        model_bkg.powerlaw.norm.values = [0, -1]  # this needs to be tested
        model_bkg.bbodyrad.norm.values = [0, -1]  # this needs to be tested
    nH = srcmod.TBabs.nH
    nHs_frozen = [nH]
    nHs_modelled = [0]
    return transformations, nHs_frozen, nHs_modelled, 'apl_bb'

def apl_bb_simple(Model, AllModels, bxa, galnh, Z, n):
    # frozen parameters
    transformations = []
    srcmod = Model('tbabs*(pow+bbodyrad)', modName='srcmod', sourceNum=1)
    srcmod.TBabs.nH.values = [galnh, -1, 0, 1e-5, 1e2, 1e2]

    # fit parameters
    gamma = srcmod.powerlaw.PhoIndex
    gamma.values = [1, 0.01, -2, -2, 4, 4]
    p_gamma = modded_create_uniform_prior_for(srcmod, gamma)
    p_gamma['name'] = '$\\alpha$' #was Gamma before but apparently xspec
    # uses alpha
    transformations.append(p_gamma)
    kT = srcmod.bbodyrad.kT
    kT.values = [0.1, 0.005, 0.001, 0.001, 0.3, 10]
    p_kT = modded_create_uniform_prior_for(srcmod, kT)
    transformations.append(p_kT)
    for groupid in range(n):
        model = AllModels(groupNum=2*groupid+1, modName='srcmod')
        norm = model.powerlaw.norm
        norm.values = [1e-4, 0.01, 0, 1e-8, 1e2, 1e20]
        p_norm = modded_create_jeffreys_prior_for(model, norm)
        p_norm['name'] = 'log(norm$_{PL,%s}$)' % (groupid+1)
        transformations.append(p_norm)
        norm_bb = model.bbodyrad.norm
        norm_bb.values = [1e-4, 0.01, 0, 1e-3, 1e7, 1e20]
        p_norm_bb = modded_create_jeffreys_prior_for(model, norm_bb)
        p_norm_bb['name'] = 'log(norm$_{BB,%s}$)' % (groupid+1)
        transformations.append(p_norm_bb)
        model_bkg = AllModels(groupNum=2*groupid+2, modName='srcmod')
        model_bkg.powerlaw.norm.values = [0, -1]  # this needs to be tested
        model_bkg.bbodyrad.norm.values = [0, -1]  # this needs to be tested
    nH = srcmod.TBabs.nH
    nHs_frozen = [nH]
    nHs_modelled = []
    return transformations, nHs_frozen, nHs_modelled, 'apl_bb_simple'

def adiskbb_simple(Model, AllModels, bxa, galnh, Z, n):
    # frozen parameters
    transformations = []
    srcmod = Model('tbabs*(diskbb)', modName='srcmod', sourceNum=1)
    srcmod.TBabs.nH.values = [galnh, -1, 0, 1e-5, 1e2, 1e2]

    # fit parameters
    Tin = srcmod.diskbb.Tin
    Tin.values = [0.1, 0.005, 0.001, 0.001, 10, 10]
    p_Tin = modded_create_uniform_prior_for(srcmod, Tin)
    p_Tin['name'] = '$T_{in}$'
    transformations.append(p_Tin)
    for groupid in range(n):
        model = AllModels(groupNum=2*groupid+1, modName='srcmod')
        norm_disk = model.diskbb.norm
        norm_disk.values = [1e-4, 0.01, 0, 1e-6, 1e7, 1e20]
        p_norm_disk = modded_create_jeffreys_prior_for(model, norm_disk)
        p_norm_disk['name'] = 'log(norm$_{disk,%s}$)' % (groupid+1)
        transformations.append(p_norm_disk)
        model_bkg = AllModels(groupNum=2*groupid+2, modName='srcmod')
        model_bkg.diskbb.norm.values = [0, -1]  # this needs to be tested
    nH = srcmod.TBabs.nH
    nHs_frozen = [nH]
    nHs_modelled = []
    return transformations, nHs_frozen, nHs_modelled, 'adiskbb_simple'

def adiskbb(Model, AllModels, bxa, galnh, Z, n):
    # frozen parameters
    transformations = []
    srcmod = Model('tbabs*tbvarabs*(diskbb)', modName='srcmod', sourceNum=1)
    srcmod.TBabs.nH.values = [galnh, -1, 0, 1e-5, 1e2, 1e2]
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
    loc_nh.values = [galnh, 0.01, 0, 1e-5, 1e2, 1e2]
    p_loc_nh = modded_create_jeffreys_prior_for(srcmod, loc_nh)
    transformations.append(p_loc_nh)
    Tin = srcmod.diskbb.Tin
    Tin.values = [0.1, 0.005, 0.001, 0.001, 10, 10]
    p_Tin = modded_create_uniform_prior_for(srcmod, Tin)
    p_Tin['name'] = '$T_{in}$'
    transformations.append(p_Tin)
    for groupid in range(n):
        model = AllModels(groupNum=2*groupid+1, modName='srcmod')
        norm_disk = model.diskbb.norm
        norm_disk.values = [1e-4, 0.01, 0, 1e-6, 1e7, 1e20]
        p_norm_disk = modded_create_jeffreys_prior_for(model, norm_disk)
        p_norm_disk['name'] = 'log(norm$_{disk,%s}$)' % (groupid+1)
        transformations.append(p_norm_disk)
        model_bkg = AllModels(groupNum=2*groupid+2, modName='srcmod')
        model_bkg.diskbb.norm.values = [0, -1]  # this needs to be tested
    nH = srcmod.TBabs.nH
    nHs_frozen = [nH]
    nHs_modelled = [0]
    return transformations, nHs_frozen, nHs_modelled, 'adiskbb'

def apl_pl_bb_simple(Model, AllModels, bxa, galnh, Z, n):
    # frozen parameters
    transformations = []
    srcmod = Model('tbabs*(pow+pow+bbodyrad)', modName='srcmod', sourceNum=1)
    srcmod.TBabs.nH.values = [galnh, -1, 0, 1e-5, 1e2, 1e2]

    # fit parameters
    gamma = srcmod.powerlaw.PhoIndex
    gamma.values = [1, 0.01, -2, -2, 4, 4]
    p_gamma = modded_create_uniform_prior_for(srcmod, gamma)
    p_gamma['name'] = '$\\alpha$' #was Gamma before but apparently xspec
    # uses alpha
    transformations.append(p_gamma)
    gamma_hard = srcmod.powerlaw_3.PhoIndex
    gamma_hard.values = [1, 0.01, -3, -3, 1, 1]
    p_gamma_hard = modded_create_uniform_prior_for(srcmod, gamma_hard)
    p_gamma_hard['name'] = '$\\alpha$ hard' #was Gamma before but apparently xspec
    # uses alpha
    transformations.append(p_gamma_hard)
    kT = srcmod.bbodyrad.kT
    kT.values = [0.1, 0.005, 0.001, 0.001, 0.3, 10]
    p_kT = modded_create_uniform_prior_for(srcmod, kT)
    transformations.append(p_kT)
    for groupid in range(n):
        model = AllModels(groupNum=2*groupid+1, modName='srcmod')
        norm = model.powerlaw.norm
        norm.values = [1e-4, 0.01, 0, 1e-8, 1e2, 1e20]
        p_norm = modded_create_jeffreys_prior_for(model, norm)
        p_norm['name'] = 'log(norm$_{PL,%s}$)' % (groupid+1)
        transformations.append(p_norm)
        norm_hard = model.powerlaw_3.norm
        norm_hard.values = [1e-4, 0.01, 0, 1e-8, 1e2, 1e20]
        p_norm_hard = modded_create_jeffreys_prior_for(model, norm_hard)
        p_norm_hard['name'] = 'log(norm$_{PL hard,%s}$)' % (groupid+1)
        transformations.append(p_norm_hard)
        norm_bb = model.bbodyrad.norm
        norm_bb.values = [1e-4, 0.01, 0, 1e-3, 1e7, 1e20]
        p_norm_bb = modded_create_jeffreys_prior_for(model, norm_bb)
        p_norm_bb['name'] = 'log(norm$_{BB,%s}$)' % (groupid+1)
        transformations.append(p_norm_bb)
        model_bkg = AllModels(groupNum=2*groupid+2, modName='srcmod')
        model_bkg.powerlaw.norm.values = [0, -1]  # this needs to be tested
        model_bkg.bbodyrad.norm.values = [0, -1]  # this needs to be tested
    nH = srcmod.TBabs.nH
    nHs_frozen = [nH]
    nHs_modelled = []
    return transformations, nHs_frozen, nHs_modelled, 'apl_pl_bb_simple'

def apl_bh(Model, AllModels, bxa, galnh, Z, n):
    # frozen parameters
    transformations = []
    srcmod = Model('tbabs*tbvarabs*pow', modName='srcmod', sourceNum=1)
    srcmod.TBabs.nH.values = [galnh, -1, 0, 1e-5, 1e2, 1e2]
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
    loc_nh.values = [galnh, 0.01, 0, 1e-5, 1e2, 1e2]
    p_loc_nh = modded_create_jeffreys_prior_for(srcmod, loc_nh)
    transformations.append(p_loc_nh)
    gamma = srcmod.powerlaw.PhoIndex
    gamma.values = [1, 0.01, 1, 1, 3, 3]
    p_gamma = modded_create_uniform_prior_for(srcmod, gamma)
    p_gamma['name'] = '$\\alpha$' #was Gamma before but apparently xspec
    # uses alpha
    transformations.append(p_gamma)
    for groupid in range(n):
        model = AllModels(groupNum=2*groupid+1, modName='srcmod')
        norm = model.powerlaw.norm
        norm.values = [1e-4, 0.01, 0, 1e-8, 1e2, 1e20]
        p_norm = modded_create_jeffreys_prior_for(model, norm)
        p_norm['name'] = 'log(norm$_{PL,%s}$)' % (groupid+1)
        transformations.append(p_norm)
        model_bkg = AllModels(groupNum=2*groupid+2, modName='srcmod')
        model_bkg.powerlaw.norm.values = [0, -1]  # this needs to be tested
    nH = srcmod.TBabs.nH
    nHs_frozen = [nH]
    nHs_modelled = [0]
    return transformations, nHs_frozen, nHs_modelled, 'apl_bh'


def apl_bh_simple(Model, AllModels, bxa, galnh, Z, n):
    # frozen parameters
    transformations = []
    srcmod = Model('tbabs*pow', modName='srcmod', sourceNum=1)
    srcmod.TBabs.nH.values = [galnh, -1, 0, 1e-5, 1e2, 1e2]

    # fit parameters
    gamma = srcmod.powerlaw.PhoIndex
    gamma.values = [1, 0.01, 1, 1, 3, 3]
    p_gamma = modded_create_uniform_prior_for(srcmod, gamma)
    p_gamma['name'] = '$\\alpha$'
    transformations.append(p_gamma)
    for groupid in range(n):
        model = AllModels(groupNum=2*groupid+1, modName='srcmod')
        norm = model.powerlaw.norm
        norm.values = [1e-4, 0.01, 0, 1e-8, 1e2, 1e20]
        p_norm = modded_create_jeffreys_prior_for(model, norm)
        p_norm['name'] = 'log(norm$_{PL,%s}$)' % (groupid+1)
        transformations.append(p_norm)
        model_bkg = AllModels(groupNum=2*groupid+2, modName='srcmod')
        model_bkg.powerlaw.norm.values = [0, -1]  # this needs to be tested
    nH = srcmod.TBabs.nH
    nHs_frozen = [nH]
    nHs_modelled = []
    return transformations, nHs_frozen, nHs_modelled, 'apl_bh_simple'

def apl_bh_wide(Model, AllModels, bxa, galnh, Z, n):
    # frozen parameters
    transformations = []
    srcmod = Model('tbabs*tbvarabs*pow', modName='srcmod', sourceNum=1)
    srcmod.TBabs.nH.values = [galnh, -1, 0, 1e-5, 1e2, 1e2]
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
    loc_nh.values = [galnh, 0.01, 0, 1e-5, 1e2, 1e2]
    p_loc_nh = modded_create_jeffreys_prior_for(srcmod, loc_nh)
    transformations.append(p_loc_nh)
    gamma = srcmod.powerlaw.PhoIndex
    gamma.values = [1, 0.01, 1, 1, 4, 4]
    p_gamma = modded_create_uniform_prior_for(srcmod, gamma)
    p_gamma['name'] = '$\\alpha$' #was Gamma before but apparently xspec
    # uses alpha
    transformations.append(p_gamma)
    for groupid in range(n):
        model = AllModels(groupNum=2*groupid+1, modName='srcmod')
        norm = model.powerlaw.norm
        norm.values = [1e-4, 0.01, 0, 1e-8, 1e2, 1e20]
        p_norm = modded_create_jeffreys_prior_for(model, norm)
        p_norm['name'] = 'log(norm$_{PL,%s}$)' % (groupid+1)
        transformations.append(p_norm)
        model_bkg = AllModels(groupNum=2*groupid+2, modName='srcmod')
        model_bkg.powerlaw.norm.values = [0, -1]  # this needs to be tested
    nH = srcmod.TBabs.nH
    nHs_frozen = [nH]
    nHs_modelled = [0]
    return transformations, nHs_frozen, nHs_modelled, 'apl_bh_wide'

def apl_bh_wide_simple(Model, AllModels, bxa, galnh, Z, n):
    # frozen parameters
    transformations = []
    srcmod = Model('tbabs*pow', modName='srcmod', sourceNum=1)
    srcmod.TBabs.nH.values = [galnh, -1, 0, 1e-5, 1e2, 1e2]

    # fit parameters
    gamma = srcmod.powerlaw.PhoIndex
    gamma.values = [1, 0.01, 1, 1, 4, 4]
    p_gamma = modded_create_uniform_prior_for(srcmod, gamma)
    p_gamma['name'] = '$\\alpha$'
    transformations.append(p_gamma)
    for groupid in range(n):
        model = AllModels(groupNum=2*groupid+1, modName='srcmod')
        norm = model.powerlaw.norm
        norm.values = [1e-4, 0.01, 0, 1e-8, 1e2, 1e20]
        p_norm = modded_create_jeffreys_prior_for(model, norm)
        p_norm['name'] = 'log(norm$_{PL,%s}$)' % (groupid+1)
        transformations.append(p_norm)
        model_bkg = AllModels(groupNum=2*groupid+2, modName='srcmod')
        model_bkg.powerlaw.norm.values = [0, -1]  # this needs to be tested
    nH = srcmod.TBabs.nH
    nHs_frozen = [nH]
    nHs_modelled = []
    return transformations, nHs_frozen, nHs_modelled, 'apl_bh_wide_simple'

def apl_pcf(Model, AllModels, bxa, galnh, Z, n):
    # frozen parameters
    transformations = []
    srcmod = Model('tbabs*(partcov*tbvarabs)*pow', modName='srcmod', sourceNum=1)
    srcmod.TBabs.nH.values = [galnh, -1, 0, 1e-5, 1e2, 1e2]
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
    part_cov = srcmod.partcov.CvrFract
    part_cov.values = [0.5, 0.01, 0, 0, 1, 1]
    p_part_cov = modded_create_uniform_prior_for(srcmod, part_cov)
    p_part_cov['name'] = 'cf'
    transformations.append(p_part_cov)
    loc_nh = srcmod.TBvarabs.nH
    loc_nh.values = [galnh, 0.01, 0, 1e-5, 1e2, 1e2]
    p_loc_nh = modded_create_jeffreys_prior_for(srcmod, loc_nh)
    transformations.append(p_loc_nh)
    gamma = srcmod.powerlaw.PhoIndex
    gamma.values = [1, 0.01, -2, -2, 4, 4]
    p_gamma = modded_create_uniform_prior_for(srcmod, gamma)
    p_gamma['name'] = '$\\alpha$' #was Gamma before but apparently xspec
    # uses alpha
    transformations.append(p_gamma)
    for groupid in range(n):
        model = AllModels(groupNum=2*groupid+1, modName='srcmod')
        norm = model.powerlaw.norm
        norm.values = [1e-4, 0.01, 0, 1e-8, 1e2, 1e20]
        p_norm = modded_create_jeffreys_prior_for(model, norm)
        p_norm['name'] = 'log(norm$_{PL,%s}$)' % (groupid+1)
        transformations.append(p_norm)
        model_bkg = AllModels(groupNum=2*groupid+2, modName='srcmod')
        model_bkg.powerlaw.norm.values = [0, -1]  # this needs to be tested
    nH = srcmod.TBabs.nH
    nHs_frozen = [nH]
    nHs_modelled = [1]
    return transformations, nHs_frozen, nHs_modelled, 'apl_pcf'

def adiskbb_bh_simple(Model, AllModels, bxa, galnh, Z, n):
    # frozen parameters
    transformations = []
    srcmod = Model('tbabs*(diskbb)', modName='srcmod', sourceNum=1)
    srcmod.TBabs.nH.values = [galnh, -1, 0, 1e-5, 1e2, 1e2]

    # fit parameters
    Tin = srcmod.diskbb.Tin
    Tin.values = [0.615, 0.005, 0.001, 0.35, 10, 10]
    p_Tin = modded_create_uniform_prior_for(srcmod, Tin)
    p_Tin['name'] = '$T_{in}$'
    transformations.append(p_Tin)
    for groupid in range(n):
        model = AllModels(groupNum=2*groupid+1, modName='srcmod')
        norm_disk = model.diskbb.norm
        norm_disk.values = [324, 0.01, 0, 36, 2916, 1e20]
        p_norm_disk = modded_create_jeffreys_prior_for(model, norm_disk)
        p_norm_disk['name'] = 'log(norm$_{disk,%s}$)' % (groupid+1)
        transformations.append(p_norm_disk)
        model_bkg = AllModels(groupNum=2*groupid+2, modName='srcmod')
        model_bkg.diskbb.norm.values = [0, -1]  # this needs to be tested
    nH = srcmod.TBabs.nH
    nHs_frozen = [nH]
    nHs_modelled = []
    return transformations, nHs_frozen, nHs_modelled, 'adiskbb_bh_simple'

def adiskbb_bh(Model, AllModels, bxa, galnh, Z, n):
    # frozen parameters
    transformations = []
    srcmod = Model('tbabs*tbvarabs*(diskbb)', modName='srcmod', sourceNum=1)
    srcmod.TBabs.nH.values = [galnh, -1, 0, 1e-5, 1e2, 1e2]
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
    loc_nh.values = [galnh, 0.01, 0, 1e-5, 1e2, 1e2]
    p_loc_nh = modded_create_jeffreys_prior_for(srcmod, loc_nh)
    transformations.append(p_loc_nh)
    Tin = srcmod.diskbb.Tin
    Tin.values = [0.615, 0.005, 0.001, 0.35, 10, 10]
    p_Tin = modded_create_uniform_prior_for(srcmod, Tin)
    p_Tin['name'] = '$T_{in}$'
    transformations.append(p_Tin)
    for groupid in range(n):
        model = AllModels(groupNum=2*groupid+1, modName='srcmod')
        norm_disk = model.diskbb.norm
        norm_disk.values = [324, 0.01, 0, 36, 2916, 1e20]
        p_norm_disk = modded_create_jeffreys_prior_for(model, norm_disk)
        p_norm_disk['name'] = 'log(norm$_{disk,%s}$)' % (groupid+1)
        transformations.append(p_norm_disk)
        model_bkg = AllModels(groupNum=2*groupid+2, modName='srcmod')
        model_bkg.diskbb.norm.values = [0, -1]  # this needs to be tested
    nH = srcmod.TBabs.nH
    nHs_frozen = [nH]
    nHs_modelled = [0]
    return transformations, nHs_frozen, nHs_modelled, 'adiskbb_bh'