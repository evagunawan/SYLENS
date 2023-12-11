import random

def subsample_single(Read_ID_List: list, Subsample_Level: int, seed: int):
    random.seed(seed)
    
    # TODO add functionality to randomly select reads to keep in subsample
    SubsampleList = []

    return SubsampleList #the function will return a list of ids to keep


def subsample_paired(Read1_ID_List: list, Read2_ID_List: list, Subsample_Level: int, seed: int):
    random.seed(seed)

    # TODO add functionality to randomly select reads to keep in subsample
    SubsampleR1List = []
    SubsampleR2List = []

    return SubsampleR1List, SubsampleR2List #the function will return two lists of ids to keep
