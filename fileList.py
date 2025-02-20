def fileList():
    filelist = []
    for r in range(1,9):
        filelist.append(f"run_{r}/example.pdb")
    return filelist