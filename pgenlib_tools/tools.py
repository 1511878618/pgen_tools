from typing import Union, List, overload
import numpy as np
import pandas as pd
import pgenlib as pg



def read_pvar(pvar):
    comment_lines = check_comment_lines(pvar)
    return pd.read_csv(pvar, sep="\t", skiprows=comment_lines)


def check_comment_lines(file, comment="##"):
    with open (file, "r") as f:
        idx = 0
        while True:
            line = f.readline()
            if not line:
                break
            if not line.startswith(comment):
                break
            idx += 1
    return idx


class PgenReaderFull():
    """
    Modify from  Pgenlib 0.90.1 
    Url:https://github.com/chrchang/plink-ng/blob/master/2.0/Python/python_api.txt

    """
    def __init__(self, pfile_path:str, pgen_path:str =None, pvar_path:str=None, sample_path:str=None ):

        self.pgen_path:str = pgen_path if pgen_path is not None else pfile_path + '.pgen'
        self.pvar_path:str = pvar_path if pvar_path is not None else pfile_path + '.pvar'
        self.sample_path:str = sample_path if sample_path is not None else pfile_path + '.psam'

        self.pvar:pg.PvarReader = pg.PvarReader(bytes(self.pvar_path, 'utf-8'))
        self._init_var() # this may be too large to load
        
        self._init_sample()

        self.pgen:pg.PgenReader = pg.PgenReader(bytes(self.pgen_path, 'utf-8'), pvar=self.pvar, 
                                                sample_subset=np.arange(len(self.sample_list), dtype=np.uint32)
                                                )


        # basic information
        self.sample_ct:int = self.get_raw_sample_ct() # number of samples as same as get_raw_sample_ct, but only call once
        self.variant_ct:int = self.get_variant_ct() # number of variants as same as get_variant_ct, but only call once


    def _init_sample(self):
        """
        load from psam and to list[[FID, IID]]
        sample_list = [[FID, IID], ...]
        """
        self.sample_list = pd.read_csv(self.sample_path, sep='\t',low_memory=False).iloc[:, :2].values.tolist()
        self.IID_list = [i[1] for i in self.sample_list]
    def _init_var(self):
        """
        load from pvar and to list[variant_id]
        """
        self.var_list = read_pvar(self.pvar_path).iloc[:, 2].tolist()


    @overload
    def get_variant_ids(self, variant_idx:List[int])->List[str]:...
    @overload
    def get_variant_ids(self, variant_idx:int)->str:...

    def get_variant_ids(self, variant_idx:Union[int, List[int]])->Union[str, List[str]]:
        """
        from PvarReader.get_variant_id 
        Modified version for interface        
        """
        if isinstance(variant_idx, int):
            return self.pvar.get_variant_id(variant_idx).decode('utf-8')
        else:
            return [self.get_variant_ids(i) for i in variant_idx]
        
    @overload
    def get_variant_idx(self, variant_id:str)->int:...
    @overload
    def get_variant_idx(self, variant_id:List[str])->List[int]:...
    def get_variant_idx(self, variant_id:Union[str, List[str]])->Union[int, List[int]]:
        """
        Modified version for interface, if not found return na
        """
        if isinstance(variant_id, str):
            try:
                return self.var_list.index(variant_id)
            except:
                return None
        else:
            return [self.get_variant_idx(i) for i in variant_id]
        
    def get_sample_ids(self, sample_idx:Union[int, List[int]])->Union[str, List[str]]:
        """
        TF version for get sample ids, if not found return na
        """
        if isinstance(sample_idx, int):
            return self.sample_list[sample_idx]
        else:
            return [self.get_sample_ids(i) for i in sample_idx]

    def get_sample_idx_by_IID(self, sample_id:Union[str, List[str]])->Union[int, List[int]]:
        """
        TF version for get sample index by IID
        """
        if isinstance(sample_id, (str, int)):
            try:
                return self.IID_list.index(sample_id)
            except:
                print(f'{sample_id} not found in IID_list')
                return None
        else:
            return [self.get_sample_idx_by_IID(i) for i in sample_id]

    def get_raw_sample_ct(self) -> int:
        """
        from PgenReader.get_raw_sample_ct
        """
        return self.pgen.get_raw_sample_ct()

    def get_variant_ct(self) -> int:
        """
        from PgenReader.get_variant_ct
        """
        return self.pgen.get_variant_ct()
    
    def hardcall_phase_present(self) -> bool:
        """
        from PgenReader.hardcall_phase_present
        """
        return self.pgen.hardcall_phase_present()
    def read(self, variant_idx, sample_idx=None, allele_idx=1)->np.ndarray:
        """
        from PgenReader.read with modified 

        Parameters
            - variant_idx:int32, index of the variant to read
            - allele_idx (int, optional): 等位基因索引，默认为 1。当 allele_idx 为 1 时，返回的基因型信息是基于ALT等位基因的计数；当 allele_idx 为 0 时，返回的基因型信息是基于REF等位基因的计数。

        """
        
        sample_idx = np.array(sample_idx) if sample_idx is not None else np.array(range(self.sample_ct))

        # use buf to store the result; the author of pgenlib think this is faster 
        buf = np.empty(self.sample_ct, dtype=np.int8)
        self.pgen.read(variant_idx = variant_idx, 
                       geno_int_out = buf,
                       allele_idx = allele_idx,
                       
                       )
        return buf[sample_idx]
    def read_range(self, variant_idx_start, variant_idx_end, allele_idx=1, sample_idx=None)->np.ndarray:
        """
        from PgenReader.read_range
        Parameters
            - variant_idx_start: int32, index of the first variant to read
            - variant_idx_end: int32, index of the last variant to read
            - allele_idx: int32, index of the allele to read, defualt 1 for ALT count, while 0 for REF count
            - sample_idx: int32, index of the sample to read
        """
        sample_idx = np.array(sample_idx) if sample_idx is not None else np.array(range(self.sample_ct))

        buf = np.empty((self.sample_ct, variant_idx_end-variant_idx_start), dtype=np.int8)
        self.pgen.read_range(variant_idx_start = variant_idx_start, 
                             variant_idx_end = variant_idx_end,
                             geno_int_out = buf,
                             allele_idx = allele_idx,
                             )
        
        return buf[:, sample_idx]
    def read_list(self, variant_idxs:list,  allele_idx=1, sample_idx=None) -> np.ndarray:
        """
        from PgenReader.read_list;

        if Sample is too large ,this will cost lots of memory; This situation, recommend to use plink2 split files
        """
        sample_idx = np.array(sample_idx) if sample_idx is not None else np.array(range(self.sample_ct))

        buf = np.empty((len(variant_idxs), self.sample_ct), dtype=np.int32) # shape (variant_nums, sample_nums)

        variant_idxs_array = np.array(variant_idxs, dtype=np.uint32)
        self.pgen.read_list(variant_idxs = variant_idxs_array, 
                            geno_int_out = buf,
                            allele_idx = allele_idx,
                            )
        return buf[:, sample_idx]

    # def read_list(self, variant_idxs:list,  allele_idx=1, sample_idx=None) -> np.ndarray:
    #     """
    #     from PgenReader.read_list
    #     """
    #     sample_idx = np.array(sample_idx) if sample_idx is not None else np.array(range(self.sample_ct))

    #     buf = np.empty((len(variant_idxs), len(sample_idx)), dtype=np.int32) # shape (variant_nums, sample_nums)
    #     print(buf.shape)
    #     variant_idxs_array = np.array(variant_idxs, dtype=np.uint32)
    #     self.pgen.read_list(variant_idxs = variant_idxs_array, 
    #                         geno_int_out = buf,
    #                         allele_idx = allele_idx,
    #                         )
    #     # return buf[:, sample_idx]
    #     return buf 

    @overload
    def extract(self, variant_idx:List[int],allele_idx=1, sample_idx=None, asFrame = False)-> np.ndarray:...
    @overload
    def extract(self, variant_idx:int,  allele_idx=1, sample_idx=None, asFrame = False)-> np.ndarray:...
    @overload
    def extract(self, variant_ids:List[str], allele_idx=1, sample_idx=None, asFrame = False)-> np.ndarray:...
    @overload
    def extract(self, variant_ids:str, allele_idx=1, sample_idx=None, asFrame = False)-> np.ndarray:...
    
    def extract(self, variant_idx=None, variant_ids=None, allele_idx=1, sample_idx=None, asFrame = False, na_rep=np.nan)->np.ndarray:
        """
        TF version for extract data from pgen file 
        Usage:
        pgfull = PgenReaderFull(pfile_path) # suffix of pgen,pvar and psam like plink2 usage do 

        # extract the first 3 variants, the second allele, and the first 3 samples
        snp_idxs = [1,2,3]
        geno, variants, samples = pgfull.extract(variant_idx=[1,2,3])

        # or extract by variant_ids
        geno, variants, samples = pgfull.extract(variant_ids=['rs1','rs2','rs3'])

        # or only part of the samples
        geno, variants, samples = pgfull.extract(variant_idx=[1,2,3], sample_idx=[1,2,3])

        # dataframe also can be returned, but not as default, so you should call like this
        geno_df = pd.DataFrame(geno_array.astype(np.int8), index=var_id, columns = [iid for fid, iid in sample_id]).T
        """
        if variant_idx is None and variant_ids is None:
            print(f'will extract all data from pgen file, numpy array shape ({self.get_variant_ct()}, {self.get_raw_sample_ct()}) ,this may cause memory error')
            variant_idx = range(self.get_variant_ct())

        if variant_idx and variant_ids:
            raise ValueError('variant_idx and variant_ids cannot be assigned at the same time')
        if variant_idx:
            if isinstance(variant_idx, int):
                variant_idx = [variant_idx]

            extracted_geno = self.read_list(variant_idx, allele_idx, sample_idx)
            samples = self.get_sample_ids(sample_idx) if sample_idx else self.sample_list
            variants = self.get_variant_ids(variant_idx)

        if variant_ids:
            if isinstance(variant_ids, str):
                variant_ids = [variant_ids]

            variant_idx = self.get_variant_idx(variant_ids)

            # if any variant_idx is None means this variants not in dataset
            variant_idx_not_found = [index for index in range(len(variant_idx)) if variant_idx[index] is None]
            if len(variant_idx_not_found) >0:
                not_found_variants = [variant_ids[i] for i in variant_idx_not_found]
            
                print(f'Totally {len(not_found_variants)} not founded, part of them are {",".join(not_found_variants[:10])}.')
                variant_idx = [i for i in variant_idx if i is not None]
                variant_ids = [i for i in variant_ids if i not in not_found_variants]

            extracted_geno = self.read_list(variant_idx, allele_idx, sample_idx)
            samples = self.get_sample_ids(sample_idx) if sample_idx else self.sample_list
            variants = variant_ids

        # return pd.DataFrame(extracted_geno, index=variants, columns=samples)

        if na_rep is not None:
            extracted_geno = np.where(extracted_geno == -9, na_rep, extracted_geno)

        if not asFrame:
            return extracted_geno, variants, samples
        else:
            return pd.DataFrame(extracted_geno, index=variants, columns=[iid for fid, iid in samples]).T

    def read_snp_list(self, snp_list, allele_idx=1, sample_idx=None)->np.ndarray:
        pass 


    def read_dosage(self):
        """
        soon modify from PgenReader.read_dosage
        """
        raise NotImplementedError
    def read_alleles(self):
        """
        soon modify from PgenReader.read_alleles
        """
        raise NotImplementedError
    
    def read_alleles_and_phasepresent(self):
        """
        soon modify from PgenReader.read_alleles_and_phasepresent
        """
        raise NotImplementedError
    

    
    def __enter__(self):
        return self
    def __exit__(self):
        self.close()

    def close(self):
        """
        from PgenReader.close
        """
        self.pgen.close()
        self.pvar.close()
