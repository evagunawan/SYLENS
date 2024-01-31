import random

def subsample_single(Read_ID_List: list, Subsample_Level: int, seed: int):
    random.seed(seed)
    
    SubsampleList = []

    SubsampleList = random.sample(Read_ID_List, Subsample_Level)

    return SubsampleList #the function will return a list of ids to keep


def subsample_paired(Read1_ID_List: list, Read2_ID_List: list, Subsample_Level: int, seed: int):
    random.seed(seed)

    SubsampleR1List = []
    SubsampleR2List = []

    Zipped_IDs = list(zip(Read1_ID_List,Read2_ID_List))
    Random_IDs = random.sample(Zipped_IDs, Subsample_Level)

    SubsampleR1List, SubsampleR2List = zip(*Random_IDs)

    return SubsampleR1List, SubsampleR2List #the function will return two lists of ids to keep