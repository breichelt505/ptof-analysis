# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 18:47:50 2023

@author: blr
"""
import h5py
import numpy as np
import pandas as pd
import pathlib

def get_setup_info(dir, DIM):
    dir = pathlib.Path(dir)
    for f in dir.iterdir():
        if DIM in str(f) and "Traveler" in str(f) and "~$" not in str(f):
            traveler = f
    data = pd.read_excel(traveler,converters={"Unnamed: 4":str,"Unnamed: 5":str})

    data_dict = {
        "snout drawing #" : data["Unnamed: 5"][8],
        "snout config" : data["Unnamed: 5"][9],
        "PTOF detector SN" : data["Unnamed: 5"][12],
        "Snout cable SN" : data["Unnamed: 5"][13],
        "Snout cable delay" : get_delay_from_SN(data["Unnamed: 5"][13])
    }

    # if np.isnan(data_dict["snout drawing #"]):
    #     data_dict = {
    #         "snout drawing #" : data["Unnamed: 4"][8],
    #         "snout config" : data["Unnamed: 4"][9],
    #         "PTOF detector SN" : data["Unnamed: 4"][12],
    #         "Snout cable SN" : data["Unnamed: 4"][13],
    #         "Snout cable delay" : get_delay_from_SN(data["Unnamed: 4"][13])
    #     }


    for f in dir.iterdir():
        if "INSTALLED-PART-CONFIG" in str(f) and "~$" not in str(f):
            installed_parts = f
    data = pd.read_csv(installed_parts)
    DLP_SN = data[data["Location"] == f"TC{DIM}-1-DLP"][data["Internal Loc"]=="PCD1CABLE"]["Serial"].values[0]
    data_dict.update({
        "DLP cable SN" : DLP_SN,
        "DLP cable delay" : get_delay_from_SN(DLP_SN.split("-")[0])
    })
    for i in range(1,4):
        try:
            atten_SN = data[data["Location"] == f"TC{DIM}-1-DIM-PCDMEZ"][data["Internal Loc"]==f"BIAST2NAT{i}"]["Serial"].values[0]
            data_dict.update({
            f"attenuator {i} SN" : atten_SN
            })
        except:
            continue
    return data_dict
    

def get_delay_from_SN(serial_num):
    table_dir = pathlib.Path("./tables/cable_response/component_timing_info")
    delay=np.nan
    for f in table_dir.iterdir():
        temp_delay = np.nan
        if "~$" not in str(f):
            data = pd.read_excel(f,sheet_name="PTOf cables",converters={"NIF Serial #":str,"Time (ns)":str})
            try:
                temp_delay = float(str(data[data["NIF Serial #"]==serial_num]["Time (ns)"].iloc[0]).replace("ns",""))
            except IndexError:
                pass
        if not np.isnan(temp_delay):
            delay = temp_delay
    return delay

def read_tek_file(fname):
    with h5py.File(fname) as f:
        channels = f["DATA"]["SHOT"].keys()
        ts = []
        Vs = []
        sat_flags = []
        vclip_lbs = []
        vclip_ubs = []

        for channel in channels:
            ts.append(np.array(f["DATA"]["SHOT"][channel]["SCALED"]["X_AXIS"])*1e9)
            Vs.append(np.array(f["DATA"]["SHOT"][channel]["SCALED"]["DATA"]))
            sat_flags.append(np.array(f["DATA"]["SHOT"][channel]["SCALED"]["DATA"]))

            range = float(f["SETTINGS"][channel]["VERTICAL_RANGE"][:])
            pos = float(f["SETTINGS"][channel]["VERTICAL_POSITION"][:])
            vclip_lbs.append(-(range/2 + pos))
            vclip_ubs.append((range/2 - pos))


    return ts, Vs, sat_flags, np.r_[vclip_lbs], np.r_[vclip_ubs]
