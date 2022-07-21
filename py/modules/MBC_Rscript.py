#!/usr/bin/env python3

import subprocess
import numpy as np

def MBC_Rscript(ref_file,hst_file,sim_file,bc_file,time,path_to_Rscript,
    method="QDM",return_results=True):

    subprocess_call = [
        "Rscript",
        path_to_Rscript,
        ref_file,
        hst_file,
        sim_file,
        bc_file,
        time,
        method,
        ]

    subprocess.call(subprocess_call)

    if return_results:

        arr_import = np.load(bc_file)

        hst_bc = np.array(arr_import[0,...])
        sim_bc = np.array(arr_import[1,...])

        return {
            "hst_bc": hst_bc,
            "sim_bc": sim_bc,
            }