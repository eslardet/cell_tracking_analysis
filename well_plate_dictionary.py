all_wells = ["B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", 
             "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11",
             "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11",
             "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11",
             "F2", "F3", "F4", "F5", "F6", "F7",
             "G2", "G3", "G4", "G5", "G6", "G7"]

wells_clean_green = ["B2", "B3", "B4", "B7", "B8", "B9", "B10", "B11", 
              "C5", "C6", "C7", "C8", "C9", "C10", "C11",
             "D2", "D3", "D4", "D5", "D6", "D7",
             "E3", "E4", "E5", "E8", "E9", "E10", "E11",
             "F2", "F3", "F4", "F5", "F6", "F7",
             "G2", "G3", "G4", "G5", "G6", "G7"]

all_groups = ["Only glial cells", "Uninfected healthy control", "AC lPVL", "AC hPVL", "HAM"]
all_stim = ["With stimulation", "No stimulation"]
all_t_cell = ["Non-specific CD4", "Specific CD4", "Non-specific CD8", "Specific CD8"]

patient_group_dict = {"Only glial cells": ["F6", "F7", "G6", "G7"],
                      "Uninfected healthy control": ["F2", "F3", "F4", "F5","G2", "G3", "G4", "G5"],
                      "AC lPVL": ["B2", "B3", "C2", "C3", "D2", "D3", "E2", "E3"],
                      "AC hPVL": ["B4", "B5", "B6", "B7", "C4", "C5", "C6", "C7", "D4", "D5", "D6", "D7", "E4", "E5", "E6", "E7"],
                      "HAM": ["B8", "B9", "B10", "B11",  "C8", "C9", "C10", "C11", "D8", "D9", "D10", "D11", "E7", "E8", "E9", "E10", "E11"]}

stimulation_dict = {"With stimulation": ["B2", "B4", "B6", "B8", "B10", "C2", "C4", "C6", "C8", "C10", "D2", "D4", "D6", "D8", "D10", 
                                         "E2", "E4", "E6", "E8", "E10", "F2", "F4", "G2", "G4"],
                    "No stimulation": ["B3", "B5", "B7", "B9", "B11", "C3", "C5", "C7", "C9", "C11", "D3", "D5", "D7", "D9", "D11", 
                                         "E3", "E5", "E7", "E9", "E11", "F3", "F5", "G3", "G5"]}

## Check whether F4, F5, G4, G5 are CD4 or CD8????
t_cell_dict = {"Non-specific CD4": ["B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "F2", "F3"],
               "Specific CD4": ["C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "G2", "G3"],
               "Non-specific CD8": ["D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "F4", "F5"],
               "Specific CD8": ["E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "G4", "G5"]}

plate_dict = {"Astrocytes & T-cells": ["1722", "1727", "1736"],
              "Microglia & T-cells": ["1723", "1728", "1737"],
              "T-cells only": ["1724", "1729", "1738"]}



def get_wells(group, stim, t_cell, clean_green=False):
    # Example:
    # inputs: group = "HAM", stim = "With stimulation", t_cell = "Specific CD4"
    # output: ['C8', 'C10']
    if clean_green == True:
        set_wells = set(wells_clean_green)
    else:
        set_wells = set(all_wells)
    return list(set.intersection(set(patient_group_dict[group]), set(stimulation_dict[stim]), set(t_cell_dict[t_cell]), set_wells))

def get_well_info(well):
    well_group = ""
    well_stim = ""
    well_t_cell = ""
    for group, w in patient_group_dict.items():
        if well in w:
            well_group = group
    for stim, w in stimulation_dict.items():
        if well in w:
            well_stim = stim
    for t_cell, w in t_cell_dict.items():
        if well in w:
            well_t_cell = t_cell
    
    return well_group, well_t_cell, well_stim

def get_plates(cell_types):
    return plate_dict[cell_types]

def get_plate_info(plate):
    cell_types = ""
    for group, p in plate_dict.items():
        if plate in p:
            cell_types = group
    return cell_types

# print(get_well_info("F6"))

# print(get_wells("HAM", "With stimulation", "Specific CD4"))
# print(get_wells_clean_green("HAM", "With stimulation", "Specific CD4"))