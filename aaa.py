import MGSimHelps as MG

def readerLammpsFile(file):
    system = MG.AtomsSystem()
    command = None
    counter_frames = 1
    for line in file:
        line = ""  # TODO: clear
        keyvalue = line.strip().partition(": ")
        if keyvalue[0] == "ITEM":
            if command == "atoms":
                yield system
                system = MG.AtomsSystem()
                counter_frames += 1

