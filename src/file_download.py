import webbrowser
import os
import pathlib
import pandas as pd
import shutil

webdav_url = "https://nifit.llnl.gov/ArchiveWebDav/export/shotdata/shots/"
lpom_url = "https://lpomprod.llnl.gov/expdir/"
nif_server_path = pathlib.Path(r"\\nif-prgms-srvr\NIFROOT\NIF\Target Diagnostic Systems\TD Factory_Alignment & Build Data\Completed Shots")

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
        f"NEUT-ALLDIAG-VALUES_{shotnum}.csv",
        f"BT-BW-ALLDIAG-VALUES_{shotnum}.csv",
        f"SHOT_OVERVIEW_{shotnum}.csv",
        f"AUTHORIZED-VALUES_{shotnum}.csv",
        f"AUTH-INSTR-BANGTIME_{shotnum}.csv",
        f"AUTH-INSTR-BURNWIDTH_{shotnum}.csv",
        f"AUTH-INSTR-TION_{shotnum}.csv",
        f"AUTH-INSTR-YN_{shotnum}.csv",
        f"AUTHORIZED-VALUES_{shotnum}.csv",
        ]
    
    for fname in fnames:
        print(url+fname)
        file = webbrowser.get().open(url+fname,autoraise=False)

def download_LPOM_files(shotnum):
    shotnum_url = f"{shotnum}/data/"
    url = lpom_url + shotnum_url
    fnames = [
        f"{shotnum}_Global_Arrival_Summary.tsv",
        f"{shotnum}_Timing_Summary.tsv",
        f"Target_Powers/{shotnum}_Total_Power.tsv"
        ]
    
    for fname in fnames:
        print(url+fname)
        file = webbrowser.get().open(url+fname,autoraise=False)

def get_exp_ID(shotnum, fpath):
    file = fpath / pathlib.Path(f"SHOT_OVERVIEW_{shotnum}.csv")
    data = pd.read_csv(file)
    return str(data["CMT Experiment"][0])


def download_nif_server_files(shotnum, exp_id, download_loc, DIM="090-078"):
    shotnum_path = nif_server_path / f"20{shotnum[1:3]}\\"
    download_loc = pathlib.Path(download_loc)

    shot_dir = ""
    for d in shotnum_path.iterdir():
        if exp_id[:-1] in str(d):
            shot_dir = d

    shot_dir = shot_dir / f"{int(DIM[0:3]):d}-{int(DIM[4:]):d} SNOUT"
    shot_dir = shot_dir / "pre-shot"
    for f in shot_dir.iterdir():
        if f.suffix == ".pdf":
            shutil.copy(f,download_loc/f"{f.stem}_{DIM}{f.suffix}")
        elif "Traveler" in str(f):
            shutil.copy(f,download_loc/f"{f.stem}_{DIM}{f.suffix}")