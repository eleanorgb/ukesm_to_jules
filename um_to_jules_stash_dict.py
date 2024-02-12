def get_jules_stash_dict():
    DICT_STASH = {
        "m01s00i627": {"name": "pop_den", "units": "1/km2"},
        "m01s00i252": {"name": "co2_mmr", "units": None},
        "m01s00i447": {"name": "deposition_n", "units": "kg/m2/s"},
        "m01s19i111": {"name": "deposition_n", "units": "kg/m2/s"},
        "m01s21i097": {"name": "flash_rate", "units": "1/km2"},
        "m01s50i082": {"name": "flash_rate", "units": "1/km2"},
        "m01s00i448": {"name": "frac_agr", "units": None},
        "m01s00i458": {"name": "frac_past", "units": None},
        "m01s00i449": {"name": "frac_agr_prev", "units": None},
        "m01s00i459": {"name": "frac_past_prev", "units": None},
        "m01s00i505": {"name": "landfrac", "units": None},
        "m01s00i157": {"name": "gridarea", "units": "m2"},
        "m01s00i030": {"name": "land_sea_mask", "units": None},
        "m01s00i151": {"name": "river_seq", "units": None},
        "m01s00i152": {"name": "river_dir", "units": None},
        "m01s00i153": {"name": "rivers_sto_rp", "units": "m2"},
        "m01s00i287": {"name": "wood_prod_fast", "units": "kg/m2"},
        "m01s00i288": {"name": "wood_prod_med", "units": "kg/m2"},
        "m01s00i289": {"name": "wood_prod_slow", "units": "kg/m2"},
        "m01s00i004": {"name": "theta", "units": "K"},
        "m01s00i033": {"name": "orography", "units": "m"},
        "m01s00i010": {"name": "spec_hum", "units": ""},
        "m01s30i111": {"name": "temp", "units": ""},
        "m01s30i417": {"name": "pstar", "units": ""},
        "m01s02i207": {"name": "lw_down", "units": ""},
        "m01s01i235": {"name": "sw_down", "units": ""},
        "m01s04i203": {"name": "ls_rain", "units": "kg/m2/s"},
        "m01s04i204": {"name": "ls_snow", "units": "kg/m2/s"},
        "m01s05i205": {"name": "conv_rain", "units": "kg/m2/s"},
        "m01s05i206": {"name": "conv_snow", "units": "kg/m2/s"},
        "m01s05i215": {"name": "snowfall", "units": "kg/m2/s"},
        "m01s05i216": {"name": "precip", "units": "kg/m2/s"},
        "m01s03i230": {"name": "wind", "units": "m/s"},
        "m01s01i231": {"name": "diffuse_sw_down", "units": ""},
        "m01s00i207": {"name": "b", "units": ""},
        "m01s00i274": {"name": "topidx", "units": ""},
        "m01s00i275": {"name": "sd_topidx", "units": ""},
        "m01s00i276": {"name": "fexp", "units": ""},
        "m01s00i040": {"name": "sm_wilt", "units": ""},
        "m01s00i041": {"name": "sm_crit", "units": ""},
        "m01s00i043": {"name": "sm_sat", "units": ""},
        "m01s00i044": {"name": "satcon", "units": ""},
        "m01s00i046": {"name": "hcap", "units": ""},
        "m01s00i047": {"name": "hcon", "units": ""},
        "m01s00i048": {"name": "sathh", "units": ""},
        "m01s00i220": {"name": "soil_albedo", "units": ""},
        "m01s00i418": {"name": "clay", "units": ""},
        "m01s00i442": {"name": "ns_dpm", "units": "kg/m2"},
        "m01s00i443": {"name": "ns_rpm", "units": "kg/m2"},
        "m01s00i444": {"name": "ns_bio", "units": "kg/m2"},
        "m01s00i445": {"name": "ns_hum", "units": "kg/m2"},
        "m01s00i446": {"name": "n_inorg", "units": "kg/m2"},
        "m01s00i466": {"name": "cs_dpm", "units": "kg/m2"},
        "m01s00i467": {"name": "cs_rpm", "units": "kg/m2"},
        "m01s00i468": {"name": "cs_bio", "units": "kg/m2"},
        "m01s00i469": {"name": "cs_hum", "units": "kg/m2"},
        "m01s00i281": {"name": "sthzw", "units": ""},
        "m01s00i278": {"name": "zw", "units": "m"},
        "m01s00i229": {"name": "canopy", "units": ""},
        "m01s00i233": {"name": "tstar_tile", "units": "K"},
        "m01s00i214": {"name": "sthu", "units": ""},
        "m01s00i215": {"name": "sthf", "units": ""},
        "m01s00i020": {"name": "t_soil", "units": "K"},
        "m01s00i240": {"name": "snow_tile", "units": "kg/m2"},
        "m01s00i216": {"name": "frac", "units": ""},
        "m01s00i217": {"name": "lai", "units": ""},
        "m01s00i218": {"name": "canht", "units": ""},
        "m01s03i462": {"name": "gs", "units": ""},
        "m01s00i213": {"name": "gs", "units": ""},
        "m01s00i231": {"name": "rgrain", "units": ""},
        "m01s00i377": {"name": "rho_snow", "units": "kg/m3"},
        "m01s00i376": {"name": "snow_depth", "units": "m"},
        "m01s00i242": {"name": "snow_grnd", "units": "kg/m2"},
        "m01s00i380": {"name": "nsnow", "units": ""},
        "m01s00i381": {"name": "snow_ds", "units": "m"},
        "m01s00i382": {"name": "snow_ice", "units": "kg/m2"},
        "m01s00i383": {"name": "snow_liq", "units": "kg/m2"},
        "m01s00i384": {"name": "tsnow", "units": "K"},
        "m01s00i386": {"name": "rgrainl", "units": ""},
        "m01s00i576": {"name": "tsurf_elev_surft", "units": "K"},
    }
    return DICT_STASH
