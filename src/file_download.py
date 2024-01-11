import webbrowser
import os

webdav_url = "https://nifit.llnl.gov/ArchiveWebDav/export/shotdata/shots/"

def download_DIM_config_info(shotnum, DIM="090-078", record_sys=2):
    shotnum_url = f"20{shotnum[1:3]}/{shotnum[3:5]}/{shotnum}/installed_configuration/csv/"
    shotnum_url += f"TC{DIM}-1-DIM-PCDMEZ/"
    url = webdav_url + shotnum_url

    # Recording system files
    fnames = [
        f"SCOPE{record_sys:d}/{shotnum}_TC{DIM}-1-DIM-PCDMEZ_SCOPE{record_sys:d}.csv"
    ]
    if record_sys == 1:
        fnames += [
            f"CRT1FIDUDELAY/{shotnum}_TC{DIM}-1-DIM-PCDMEZ_CRT1FIDUDELAY.csv",
            f"CRT1OUTNAT1/{shotnum}_TC{DIM}-1-DIM-PCDMEZ_CRT1OUTNAT1.csv",
            f"CRT1OUTNAT2/{shotnum}_TC{DIM}-1-DIM-PCDMEZ_CRT1OUTNAT2.csv",
            f"CRTSCOPE1/{shotnum}_TC{DIM}-1-DIM-PCDMEZ_CRTSCOPE1.csv",
                ] 
    
    for fname in fnames:
        print(url+fname)
        file = webbrowser.get().open(url+fname,autoraise=False)


    # Snout config files
    shotnum_url = f"20{shotnum[1:3]}/{shotnum[3:5]}/{shotnum}/installed_configuration/csv/"
    shotnum_url += f"TC{DIM}-1-SNOUT/AUXBRACKET/"
    url = webdav_url + shotnum_url

    fnames = [
        f"{shotnum}_TC{DIM}-1-SNOUT_AUXBRACKET.csv"
    ]
    for fname in fnames:
        print(url+fname)
        file = webbrowser.get().open(url+fname,autoraise=False)


def download_PCD_data(shotnum, DIM="090-078", record_sys=2):
    shotnum_url = f"20{shotnum[1:3]}/{shotnum[3:5]}/{shotnum}/raw_files/TD/PCD/"
    url = webdav_url + shotnum_url
    fnames = [
        f"TD_TC{DIM}_PCD_HVPS-{record_sys:02d}-DB_{shotnum}.h5",
        f"TD_TC{DIM}_PCD_SCOPE-{record_sys:02d}-DB_{shotnum}.h5",
        f"TD_TC{DIM}_PCD_CRT-SCOPE-{record_sys:02d}-DB_{shotnum}.h5"
            ] 
                    # Note that the FTD (CRT-SCOPE) is not present on record_sys 2  
                                                                            
    for fname in fnames:
        print(url+fname)
        file = webbrowser.get().open(url+fname,autoraise=False)


def download_report_files(shotnum):
    shotnum_url = f"20{shotnum[1:3]}/{shotnum[3:5]}/{shotnum}/reports/"
    url = webdav_url + shotnum_url
    fnames = [
        f"SHOT-INPUTS_{shotnum}.csv",
        f"INSTALLED-PART-CONFIG_{shotnum}.csv",
        f"NEUT_ALLDIAG_VALUES_{shotnum}.csv",
        f"BT-BW-ALLDIAG-VALUES_{shotnum}.csv",
        f"SHOT_OVERVIEW_{shotnum}.csv"
        ]
    
    for fname in fnames:
        print(url+fname)
        file = webbrowser.get().open(url+fname,autoraise=False)

