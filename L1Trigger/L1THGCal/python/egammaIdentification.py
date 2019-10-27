import FWCore.ParameterSet.Config as cms


class Category:
    def __init__(self, eta_min, eta_max, pt_min, pt_max):
        self.eta_min = eta_min
        self.eta_max = eta_max
        self.pt_min = pt_min
        self.pt_max = pt_max


categories = [
        # Low eta
        Category(eta_min=1.5, eta_max=2.7, pt_min=0., pt_max=1e6),
        # High eta
        Category(eta_min=2.7, eta_max=3.0, pt_min=0., pt_max=1e6),
        ]

# Identification for dRNN 2D clustering and cone 3D clustering
bdt_weights_drnn_cone = [
        # Low eta
        'L1Trigger/L1THGCal/data/egamma_id_drnn_cone_loweta_v0.xml',
        # High eta
        'L1Trigger/L1THGCal/data/egamma_id_drnn_cone_higheta_v0.xml',
        ]

working_points_drnn_cone = [
        # Low eta
        {
        '900':0.14057436,
        '950':0.05661769,
        '975':-0.01481255,
        '995':-0.19656579,
        },
        # High eta
        {
        '900':0.05995301,
        '950':-0.02947988,
        '975':-0.10577436,
        '995':-0.26401181,
        },
        ]


# Identification for dRNN 2D clustering and DBSCAN 3D clustering
bdt_weights_drnn_dbscan = [
        # Low eta
        'L1Trigger/L1THGCal/data/egamma_id_drnn_dbscan_loweta_v0.xml',
        # High eta
        'L1Trigger/L1THGCal/data/egamma_id_drnn_dbscan_higheta_v0.xml',
        ]

working_points_drnn_dbscan = [
        # Low eta
        {
        '900':0.08421164,
        '950':0.06436077,
        '975':-0.04547527,
        '995':-0.23717142,
        },
        # High eta
        {
        '900':0.05559443,
        '950':-0.0171725,
        '975':-0.10630798,
        '995':-0.27290947,
        },
        ]

# Identification for 3D HistoMax clustering (as in 3.5.2)
bdt_weights_histomax = [
        # Low eta
        'L1Trigger/L1THGCal/data/egid_histomax_loweta.xml',
        # High eta
        'L1Trigger/L1THGCal/data/egid_histomax_higheta.xml',
        ]

working_points_histomax = [
        # Low eta
        {
        '995': -0.9635030, # epsilon_b = 0.1135
        '975': -0.2973018, # epsilon_b = 0.0584
        '950': 0.4912865, # epsilon_b = 0.0439
        '900': 0.8831878, # epsilon_b = 0.0305
        },
        # High eta
        {
        '995': -0.9964944, # epsilon_b = 0.4657
        '975': -0.8169104, # epsilon_b = 0.1600
        '950': -0.0011038, # epsilon_b = 0.0800
        '900': 0.6431891, # epsilon_b = 0.0514
        },
        ]

tight_wp = ['975', '900']
loose_wp = ['995', '950']


egamma_identification_drnn_cone = cms.PSet(
        Inputs=cms.vstring('cl3d_firstlayer', 'cl3d_coreshowerlength', 'cl3d_maxlayer', 'cl3d_srrmean'),
        CategoriesEtaMin=cms.vdouble([cat.eta_min for cat in categories]),
        CategoriesEtaMax=cms.vdouble([cat.eta_max for cat in categories]),
        CategoriesPtMin=cms.vdouble([cat.pt_min for cat in categories]),
        CategoriesPtMax=cms.vdouble([cat.pt_max for cat in categories]),
        Weights=cms.vstring(bdt_weights_drnn_cone),
        WorkingPoints=cms.vdouble([wps[eff] for wps,eff in zip(working_points_drnn_cone,tight_wp)]),
        )

egamma_identification_drnn_dbscan = cms.PSet(
        Inputs=cms.vstring('cl3d_firstlayer', 'cl3d_coreshowerlength', 'cl3d_maxlayer', 'cl3d_srrmean'),
        CategoriesEtaMin=cms.vdouble([cat.eta_min for cat in categories]),
        CategoriesEtaMax=cms.vdouble([cat.eta_max for cat in categories]),
        CategoriesPtMin=cms.vdouble([cat.pt_min for cat in categories]),
        CategoriesPtMax=cms.vdouble([cat.pt_max for cat in categories]),
        Weights=cms.vstring(bdt_weights_drnn_dbscan),
        WorkingPoints=cms.vdouble([wps[eff] for wps,eff in zip(working_points_drnn_dbscan,tight_wp)]),
        )

egamma_identification_histomax = cms.PSet(
        Inputs=cms.vstring('cl3d_coreshowerlength', 'cl3d_showerlength', 'cl3d_firstlayer', 'cl3d_maxlayer', 'cl3d_szz', 'cl3d_srrmean', 'cl3d_srrtot', 'cl3d_seetot', 'cl3d_spptot'),
        CategoriesEtaMin=cms.vdouble([cat.eta_min for cat in categories]),
        CategoriesEtaMax=cms.vdouble([cat.eta_max for cat in categories]),
        CategoriesPtMin=cms.vdouble([cat.pt_min for cat in categories]),
        CategoriesPtMax=cms.vdouble([cat.pt_max for cat in categories]),
        Weights=cms.vstring(bdt_weights_histomax),
        WorkingPoints=cms.vdouble([wps[eff] for wps,eff in zip(working_points_histomax,tight_wp)]),
        )
