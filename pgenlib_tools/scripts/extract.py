#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:       :
@Date     :2024/03/05 15:19:30
@Author      :Tingfeng Xu
@version      :1.0
'''

import numpy as np 
import argparse
import sys
import warnings
import textwrap
from signal import SIG_DFL, SIGPIPE, signal
import pandas as pd 
from pgenlib_tools import PgenReaderFull
from pathlib import Path
warnings.filterwarnings("ignore")
signal(
    SIGPIPE, SIG_DFL
)  # prevent IOError: [Errno 32] Broken pipe. If pipe closed by 'head'.


def getParser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            """
        %prog Not work for big data! carefully when use large data 
        @Author: xutingfeng@big.ac.cn

        GWAS SSF: https://www.biorxiv.org/content/10.1101/2022.07.15.500230v1
        GWAS Standard

        Version: 1.0
    
        """
        ),
    )
    parser.add_argument("-p",
                        "--pgen",
                        dest="pfile",
                        required=True,
    )
    parser.add_argument("-s","--snplist",
                        dest="sfile",
                        required=True,
    )
    parser.add_argument("-o","--output",
                        dest="ofile",
                        required=True,
    )

    
    return parser

if __name__ == "__main__":
    parser = getParser()
    args = parser.parse_args()

    pfile = args.pfile
    sfile = args.sfile
    ofile = args.ofile


    pgen = PgenReaderFull(pfile_path=pfile)
    variant_ids = pd.read_csv(sfile).iloc[:, 0].tolist()


    df = pgen.extract(variant_ids =variant_ids,
                    asFrame=True, na_rep=np.nan)
    df_carrier = df.sum(axis=1).to_frame().reset_index(drop=False)
    df_carrier.columns = ["eid",'carrier']
    if ofile.endswith(".csv"):
        df_carrier.to_csv(ofile, index=False)
    elif ofile.endswith(".tsv"):
        df_carrier.to_csv(ofile, index=False, sep="\t")

