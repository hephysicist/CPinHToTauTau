sample_list = [

  'DYto2L_M_10to50_amcatnloFXFX',
  'DYto2L_M_50_0J_amcatnloFXFX',
  'DYto2L_M_50_1J_amcatnloFXFX',
  'DYto2L_M_50_2J_amcatnloFXFX',
  'DYto2L_M_50_amcatnloFXFX',
  'DYto2Tau_MLL_50_0J_amcatnloFXFX',
  'DYto2Tau_MLL_50_1J_amcatnloFXFX',
  'DYto2Tau_MLL_50_2J_amcatnloFXFX',
  'GluGluHto2Tau_UncorrelatedDecay_CPodd_Filtered_ProdAndDecay',
  'GluGluHto2Tau_UncorrelatedDecay_CPodd_UnFiltered_ProdAndDecay',
  'GluGluHto2Tau_UncorrelatedDecay_MM_Filtered_ProdAndDecay',
  'GluGluHto2Tau_UncorrelatedDecay_MM_UnFiltered_ProdAndDecay',
  'GluGluHto2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay',
  'GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay',
  'TbarWplusto4Q',
  'TWminusto4Q',
  'TbarWplusto2L2Nu',
  'TbarWplustoLNu2Q',
  'TWminusto2L2Nu',
  'TWminustoLNu2Q',
  'TTto2L2Nu',
  'TTto4Q',
  'TTtoLNu2Q',
  'VBFHto2Tau_UncorrelatedDecay_Filtered',
  'VBFHto2Tau_UncorrelatedDecay_UnFiltered',
  'WminusHto2Tau_UncorrelatedDecay_Filtered',
  'WminusHto2Tau_UncorrelatedDecay_UnFiltered',
  'WplusHto2Tau_UncorrelatedDecay_Filtered',
  'WplusHto2Tau_UncorrelatedDecay_UnFiltered',
  'WtoLNu_1J_madgraphMLM',
  'WtoLNu_2J_madgraphMLM',
  'WtoLNu_3J_madgraphMLM',
  'WtoLNu_4J_madgraphMLM',
  'WtoLNu_madgraphMLM',
  'WW',
  'WWW_4F',
  'WWZ_4F',
  'WZ',
  'WZZ',
  'ZHto2Tau_UncorrelatedDecay_Filtered',
  'ZHto2Tau_UncorrelatedDecay_UnFiltered',
  'ZZ',
  'ZZZ',
]

sample_to_process = {

  'DYto2L_M_10to50_amcatnloFXFX': 'dy_m10to50',
  'DYto2L_M_50_0J_amcatnloFXFX': 'dy_m50toinf_0j',
  'DYto2L_M_50_1J_amcatnloFXFX': 'dy_m50toinf_1j',
  'DYto2L_M_50_2J_amcatnloFXFX': 'dy_m50toinf_2j',
  'DYto2L_M_50_amcatnloFXFX': 'dy_m50toinf',
  'DYto2Tau_MLL_50_0J_amcatnloFXFX': 'dy_tautau_m50toinf_0j',
  'DYto2Tau_MLL_50_1J_amcatnloFXFX': 'dy_tautau_m50toinf_1j',
  'DYto2Tau_MLL_50_2J_amcatnloFXFX': 'dy_tautau_m50toinf_2j',

  'GluGluHto2Tau_UncorrelatedDecay_CPodd_Filtered_ProdAndDecay': 'h_ggf_htt',
  'GluGluHto2Tau_UncorrelatedDecay_CPodd_UnFiltered_ProdAndDecay': 'h_ggf_htt',
  'GluGluHto2Tau_UncorrelatedDecay_MM_Filtered_ProdAndDecay': 'h_ggf_htt',
  'GluGluHto2Tau_UncorrelatedDecay_MM_UnFiltered_ProdAndDecay': 'h_ggf_htt',
  'GluGluHto2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay': 'h_ggf_htt',
  'GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay': 'h_ggf_htt',

  'TbarWplusto4Q': 'st_twchannel_tbar_fh', #TbarWplusto4Q_ext1
  'TWminusto4Q': 'st_twchannel_t_fh', #TWminusto4Q_ext1
  'TbarWplusto2L2Nu': 'st_twchannel_tbar_dl', #TbarWplusto2L2Nu_ext1
  'TbarWplustoLNu2Q': 'st_twchannel_tbar_sl', #TbarWplustoLNu2Q_ext1
  'TWminusto2L2Nu': 'st_twchannel_t_dl', #TWminusto2L2Nu_ext1
  'TWminustoLNu2Q': 'st_twchannel_t_sl', #TWminustoLNu2Q_ext1
  
  'TTto2L2Nu': 'tt_dl',
  'TTto4Q': 'tt_fh',
  'TTtoLNu2Q': 'tt_sl',

  'VBFHto2Tau_UncorrelatedDecay_Filtered': 'h_vbf_htt',
  'VBFHto2Tau_UncorrelatedDecay_UnFiltered': 'h_vbf_htt',

  'WminusHto2Tau_UncorrelatedDecay_Filtered': 'wh_htt',
  'WminusHto2Tau_UncorrelatedDecay_UnFiltered': 'wh_htt',
  'WplusHto2Tau_UncorrelatedDecay_Filtered': 'wh_htt',
  'WplusHto2Tau_UncorrelatedDecay_UnFiltered': 'wh_htt',

  'WtoLNu_1J_madgraphMLM': 'w_lnu_1j',
  'WtoLNu_2J_madgraphMLM': 'w_lnu_2j',
  'WtoLNu_3J_madgraphMLM': 'w_lnu_ge3j',
  'WtoLNu_4J_madgraphMLM': 'w_lnu_ge3j',
  'WtoLNu_madgraphMLM': 'w_lnu',

  'WW': 'ww',
  'WWW_4F': 'www',
  'WWZ_4F': 'wwz',
  'WZ': 'wz',
  'WZZ': 'wzz',

  'ZHto2Tau_UncorrelatedDecay_Filtered': 'zh_htt',
  'ZHto2Tau_UncorrelatedDecay_UnFiltered': 'zh_htt',

  'ZZ': 'zz',
  'ZZZ': 'zzz',
}
