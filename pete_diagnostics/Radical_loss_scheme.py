# Open globchem
def main():
    MASTER_DIR ="/mnt/lustre/users/mjr583/GC/old_versions/12.9.3/KPP/Tropchem/Tropchem.eqn"
    LINES = open(MASTER_DIR).readlines()
    
    ROx_LIST = ["A3O2",
                "ATO2",
                "B3O2",
                "BRO2",
                "C4HVP1",
                "C4HVP2",
                "ETO2",
                "HPALD1OO",
                "HPALD2OO",
                "ICHOO",
                "ICNOO",
                "IDHNBOO",
                "IDHNDOO1",
                "IDHNDOO2",
                "IDNOO",
                "IEPOXAOO",
                "IEPOXBOO",
                "IHPNBOO",
                "IHPNDOO",
                "IHOO1",
                "IHOO4",
                "IHPOO1",
                "IHPOO2",
                "IHPOO3",
                "INO2B",
                "INO2D",
                "ISOPNOO1",
                "ISOPNOO2",
                "KO2",
                "LIMO2",
                "MACR1OO",
                "MCROHOO",
                "MO2",
                "MCO3",
                "MVKOHOO",
                "OTHRO2",
                "PIO2",
                "PRN1",
                "R4N1",
                "PO2",
                "R4O2",
                "RCO3",
                "OH",
                "HO2"]

    num_list = "0123456789."

    #First loop cleans up file
    BEGIN = False
    ADD_TO_NEXT = False
    ALL_REACTS = []
    LOSS_REACTS = []
    NON_ROX_LOSS = []
    for x in range(1,len(LINES)):
        # Strip \n
        LINE = LINES[x].strip()
        #Replace white space with a single space
        LINE = " ".join(LINE.split())
       
        # Wait until header ends
        if LINE == "#EQUATIONS":
            BEGIN = True
            continue
        if not BEGIN:
            continue  
        print(LINE) 
        # If first character is a forward slash then skip
        if LINE[0] == "/":
            continue
        
        # Add the previous line if no end of equation character
        if ADD_TO_NEXT:
            LINE = f"{ADD_LINE} {LINE}"
            #Reset switch
            ADD_TO_NEXT = False    
        
        # Remove heterogenous reactions
        if "HET" in LINE:
            continue
        
        # Remove anything after ";" on a line to remove refrences
        line_end = LINE.find(";")
        if line_end > 0 : 
            LINE = LINE[:line_end]
        
        # Check for end of equation character
        if LINE.count(":") == 0:
            ADD_TO_NEXT = True
            ADD_LINE = LINE
            continue
        
        # Radicals on the left and right of reaction.
        LOSS_RADICALS = 0
        PROD_RADCIALS = 0
        
        # Split each side of equation
        REACTS = LINE [:LINE.find("=")]
        PRODS  = LINE [LINE.find("="):LINE.find(":")]
        
        # Count radical reactants
        for R in REACTS.split(" "):
            R = R.strip()
            #Cycle through each RO2 species
            for S in ROx_LIST:
                if R == S:
                    LOSS_RADICALS += 1
        
        # Count radical products
        for P in PRODS.split(" "):
            P = P.strip()
            #Remove stoichiometry 
            P_NO_STOIC = P.lstrip(num_list)
            # Cycle through each RO2 species
            for S in ROx_LIST:
                 if S == P_NO_STOIC:
                    # If after the species is remove the string is empty then 1 whole radical is present
                    STOIC = P.replace(S, "")
                    if STOIC == "":
                        PROD_RADCIALS +=1
                    else:
                        PROD_RADCIALS += float(STOIC)
        
        RADICAL_CHANGE = PROD_RADCIALS - LOSS_RADICALS
     
        #! provides a unique character to split on
        LINE = f"{LINE} ! RADICAL_CHANGE({RADICAL_CHANGE})"
        
        ALL_REACTS.append(LINE)
        #If radical change is a loss and at least 2 radical reactants save as ROx + ROx reaction
        if (RADICAL_CHANGE < 0) & (LOSS_RADICALS >= 2):
            LOSS_REACTS.append(LINE)

        if (RADICAL_CHANGE < 0) & (LOSS_RADICALS < 2):
            NON_ROX_LOSS.append(LINE)

    # Write all reactions with calculated radical loss
    with open("TROPCHEM_All_reactions_radical_loss.txt", "w") as file:
        for L in ALL_REACTS:
            file.write(L+"\n")
            
    with open("TROPCHEM_ROx_self_reactions.txt", "w") as file:
        for L in LOSS_REACTS:
            file.write(L+"\n")       

    with open("TROPCHEM_NON_ROx_self_reactions.txt", "w") as file:
        for L in NON_ROX_LOSS:
            file.write(L+"\n")

if __name__ == "__main__":
    main()
