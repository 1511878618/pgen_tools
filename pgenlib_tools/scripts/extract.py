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
                        required=False,
    )
    # Must for -g with one more chr 
    parser.add_argument("--pgen-folder",
                        dest="pfile_folder",
                        required=False,
                        help="pgen file folder")
    parser.add_argument("--pgen-format",
                        dest="pfile_format",
                        required=False,
                        default=None,
                        help="pgen file name format with chr"
                        )
    parser.add_argument("-s","--snpfile",
                        dest="sfile",
                        required=True,
                        help="first col should be snp "
    )
    # chr only work when not all snp on the same chr and with --pgen-folder 
    parser.add_argument("-c", "--chrom",
                        dest="chrom",
                        required=False,
                        default=None,
                        help="chromosome number, default None"
                        )
    parser.add_argument("-g", "--group",
                        dest="group",
                        required=False,
                        default=None,
                        help="groupby this column, default None"
                        )
    
    parser.add_argument("-o","--output",
                        dest="ofile",
                        required=True,
    )
    parser.add_argument("-m","--method",default="carrier",
                        dest="method",
                        required=False,
    )

    
    return parser

def load_data(x):
    # if isinstance(x, Path.)
    x = str(x)
    if ".csv" in x:
        return pd.read_csv(x)
    elif x.endswith(".feather"):
        return pd.read_feather(x)
    elif x.endswith(".pkl"):
        return pd.read_pickle(x)
    elif ".tsv" in x:
        return pd.read_csv(x, sep="\t")
    else:
        raise ValueError(f"File format: {x} not supported")



def check_pfile_exists(pfile):
    pfile = str(pfile)
    pgen = pfile + ".pgen"
    pvar = pfile + ".pvar"
    psam = pfile + ".psam"
    if not Path(pgen).exists():
        print(f"{pgen} not exists")
        return False

    if not Path(pvar).exists():
        print(f"{pvar} not exists")
        return False

    if not Path(psam).exists():
        print(f"{psam} not exists")
        return False
    return True
def burdenSet(df, method = "carrier"):
    """
    method == "carrier" means if one of them are carrier then set carrier as 1 
    """
    if method == "carrier":
        df_carrier = (df.sum(axis=1) >= 1).astype(int) .to_frame().reset_index(drop=False)
        df_carrier.columns = ["eid",'carrier']
        df_carrier['carrier']
        
        return df_carrier
    else:
        raise NotImplementedError(f"Not implemented for {method}")
if __name__ == "__main__":
    parser = getParser()
    args = parser.parse_args()

    pfile = args.pfile
    sfile = args.sfile
    ofile = args.ofile
    Path(ofile).parent.mkdir(parents=True, exist_ok=True)

    method = args.method
    group = args.group
    chrom = args.chrom
    pfile_folder = args.pfile_folder
    pfile_format = args.pfile_format # ukb_wes_{} 

    if pfile is None and pfile_folder is None:
        print("Please provide pgen file or pgen folder")
        sys.exit(1)
    if pfile_folder is not None:
        if chrom is None:
            print("Please provide chromosome number")
            sys.exit(1)
        
    if pfile is not None and pfile_folder is not None:
        print("Please do not provide pgen file and pgen folder at the same time")
        sys.exit(1)

    
    if pfile is not None:  # only use one pfile to extract 
        pgen = PgenReaderFull(pfile_path=pfile)
        variant_ids = pd.read_csv(sfile).iloc[:, 0].tolist()
        df = pgen.extract(variant_ids =variant_ids,
                        asFrame=True, na_rep=np.nan)
        df_carrier = burdentSet(df, method)
        if ofile.endswith(".csv"):
            df_carrier.to_csv(ofile, index=False)
        elif ofile.endswith(".tsv"):
            df_carrier.to_csv(ofile, index=False, sep="\t")
    elif pfile_folder is not None:
        if chrom is None and group is None:
            print("Please provide chromosome number and group column")
            sys.exit(1)
        # check all chrom in the pfile_folder 
        pfile_folder = Path(pfile_folder)
        pfile_format = str(pfile_folder/pfile_format)
        get_chr_file = lambda x: pfile_format.format(x)
        sfile_df = load_data(sfile)
        chrs = sfile_df[chrom].unique().tolist()

        for chr in chrs:
            if not check_pfile_exists(get_chr_file(chr)):
                raise ValueError(f"Check {get_chr_file(chr)} exists plz")
            
        res = []
        for groupName, group_df in sfile_df.groupby(group):
            print(f"Extract {groupName}")
            assert group_df['chr'].nunique() == 1, "More than 1 chr in" + groupName
            chrfile = get_chr_file(group_df['chr'].iloc[0])
            chrpgen = PgenReaderFull(pfile_path=chrfile)
            variant_ids = group_df['ID'].tolist()
            df = chrpgen.extract(variant_ids =variant_ids,
                                asFrame=True, na_rep=np.nan)
            
            # only for carrier 
            print(df.sum(axis=0))
            df_burden = burdenSet(df)
            df_burden.rename(columns = {
                "carrier": f"carrier_{groupName}"}, inplace=True
            )
            print(df_burden[f'carrier_{groupName}'].value_counts())
            if df_burden[f'carrier_{groupName}'].sum() >1:
                res.append(df_burden)
                # df_burden.to_csv(saveDir/f"{groupName}.csv", index=False)

            else:
                print(f"No carrier of {groupName}")

            print("********************************************************")
        from functools import reduce 
        res_df = reduce(lambda x, y : x.merge(y), res)
        if ofile.endswith(".csv"):
            res_df.to_csv(ofile, index=False)
        elif ofile.endswith(".tsv"):
            res_df.to_csv(ofile, index=False, sep="\t")
        else:
            res_df.to_csv(Path(ofile).with_suffix(".csv"), index=False)
            print(f"Not supported file format {Path(ofile).suffix}, save as csv format")
    print("Done")
